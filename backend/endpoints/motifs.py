import math
from typing import Optional

from fastapi import APIRouter, UploadFile, File, Form, HTTPException

from core.bg_freqs import BackgroundFrequencies
from core.fasta_reader import FastaReader, FastaRecord
from core.fimo_runner import FimoRunner
from core.p_value_calc import PValueCalculator
from core.pwm import PositionWeightMatrix
from core.seq_scanner import SequenceScanner

router = APIRouter()

ENGINE_VALUES = {"custom", "fimo", "both"}
CUSTOM_CALIBRATION_VALUES = {"legacy", "fimo_like"}

# legacy: mantem o comportamento atual do projeto
# fimo_like: aproxima a construcao de PWM do fluxo MEME/FIMO para motivo IUPAC com nsites=20
CUSTOM_CALIBRATION_PRESETS = {
    "legacy": {
        "pseudo_fraction": 0.0001,
        "lambda_scale": 100.0,
        "score_mode": "integer",
    },
    "fimo_like": {
        # FIMO usa pseudocount total padrao 0.1; para nsites=20 isso equivale a 0.1/20 = 0.005
        "pseudo_fraction": 0.005,
        "lambda_scale": 100.0,
        "score_mode": "log_odds",
    },
}


def _resolve_custom_calibration(mode: str) -> dict:
    preset = CUSTOM_CALIBRATION_PRESETS.get(mode)
    if not preset:
        raise ValueError(f"custom_calibration_mode invalido: {mode}")
    return preset


def _attach_bh_q_values(matches: list[dict]) -> None:
    if not matches:
        return

    m = len(matches)
    ranked_indices = sorted(
        range(m),
        key=lambda idx: (matches[idx]["p_value_numeric"], -matches[idx]["score_log_odds"]),
    )

    q_values = [1.0] * m
    running_min = 1.0

    for rank in range(m, 0, -1):
        idx = ranked_indices[rank - 1]
        p_value = float(matches[idx]["p_value_numeric"])
        adjusted = min(1.0, (p_value * m) / rank)
        if adjusted < running_min:
            running_min = adjusted
        q_values[idx] = running_min

    for idx, q_value in enumerate(q_values):
        matches[idx]["q_value_numeric"] = q_value
        matches[idx]["q_value"] = f"{q_value:.2e}"
        matches[idx]["-log10_qval"] = round(-math.log10(q_value), 2) if q_value > 0 else None


def _run_custom_scan(
        records: list[FastaRecord],
        motif_pattern: str,
        bg_freqs: dict[str, float],
        p_value_threshold: float,
        score_threshold: Optional[int],
        scan_reverse_complement: bool,
        custom_calibration_mode: str,
) -> tuple[dict, list[dict]]:
    calibration = _resolve_custom_calibration(custom_calibration_mode)
    pseudo_fraction = float(calibration["pseudo_fraction"])
    lambda_scale = float(calibration["lambda_scale"])
    score_mode = str(calibration["score_mode"])

    pfm = PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern)
    pwm_log_fwd = PositionWeightMatrix.build_log_odds_pwm(
        pfm,
        bg_freqs,
        pseudo_fraction=pseudo_fraction,
    )
    pwm_log_rev = PositionWeightMatrix.get_reverse_complement(pwm_log_fwd)
    pwm_fwd = PositionWeightMatrix.build_integer_pwm(
        pfm,
        bg_freqs,
        pseudo_fraction=pseudo_fraction,
        lambda_scale=lambda_scale,
    )
    pwm_rev = PositionWeightMatrix.get_reverse_complement(pwm_fwd)
    pcalc_fwd = PValueCalculator(pwm_fwd, bg_freqs)
    pcalc_rev = PValueCalculator(pwm_rev, bg_freqs)
    motif_length = pwm_fwd.shape[1]

    if score_threshold is None:
        cutoff_fwd = pcalc_fwd.get_score_threshold_for_pvalue(p_value_threshold)
        cutoff_rev = pcalc_rev.get_score_threshold_for_pvalue(p_value_threshold)
        cutoff_mode = "pvalue"
        resolved_score_threshold = None
    else:
        cutoff_fwd = score_threshold
        cutoff_rev = score_threshold
        cutoff_mode = "score"
        resolved_score_threshold = score_threshold

    scan_configs = [("forward", "+", pwm_fwd, pwm_log_fwd, pcalc_fwd, cutoff_fwd)]
    if scan_reverse_complement:
        scan_configs.append(("reverse", "-", pwm_rev, pwm_log_rev, pcalc_rev, cutoff_rev))

    matches: list[dict] = []
    for record in records:
        record_length = len(record.sequence)
        record_orientation = "reverse_complement" if record.is_reverse else "forward"

        for strand_name, strand_symbol, pwm, pwm_log, pcalc, cutoff in scan_configs:
            hit_positions, hit_scores = SequenceScanner.scan_integer_matrix(
                sequence=record.sequence,
                pwm=pwm,
                threshold=cutoff,
            )

            for pos_idx, score_value in zip(hit_positions, hit_scores):
                score_int = int(score_value)
                p_value = float(pcalc.get_pvalue(score_int))

                if score_threshold is None and p_value > p_value_threshold:
                    continue

                start = int(pos_idx) + 1
                end = start + motif_length - 1
                sequence_window = record.sequence[pos_idx:pos_idx + motif_length]
                matched_sequence = (
                    PositionWeightMatrix.reverse_complement_sequence(sequence_window)
                    if strand_symbol == "-"
                    else sequence_window
                )

                if record.is_reverse:
                    source_start = record_length - end + 1
                    source_end = record_length - start + 1
                    source_strand_symbol = "-" if strand_symbol == "+" else "+"
                else:
                    source_start = start
                    source_end = end
                    source_strand_symbol = strand_symbol

                source_strand = "reverse" if source_strand_symbol == "-" else "forward"
                neg_log10_pval = round(-math.log10(p_value), 2) if p_value > 0 else None

                window_bases = [PositionWeightMatrix.BASE_TO_ID[base] for base in sequence_window]
                score_log_odds = float(sum(pwm_log[base_id, col] for col, base_id in enumerate(window_bases)))
                score_output = score_int if score_mode == "integer" else round(score_log_odds, 4)

                matches.append({
                    "gene": record.gene_id,
                    "record_orientation": record_orientation,
                    "strand": strand_name,
                    "strand_symbol": strand_symbol,
                    "source_strand": source_strand,
                    "source_strand_symbol": source_strand_symbol,
                    "matched_sequence": matched_sequence,
                    "position": start,
                    "end_position": end,
                    "source_position": source_start,
                    "source_end_position": source_end,
                    "score": score_output,
                    "score_integer": score_int,
                    "score_log_odds": round(score_log_odds, 4),
                    "p_value": f"{p_value:.2e}",
                    "p_value_numeric": p_value,
                    "-log10_pval": neg_log10_pval,
                })

    matches.sort(key=lambda item: (item["p_value_numeric"], -item["score_log_odds"]))
    _attach_bh_q_values(matches)

    metadata = {
        "motif_pattern": motif_pattern,
        "motif_length": motif_length,
        "score_threshold": resolved_score_threshold,
        "cutoff_mode": cutoff_mode,
        "background_frequencies": bg_freqs,
        "scan_reverse_complement": scan_reverse_complement,
        "custom_calibration_mode": custom_calibration_mode,
        "custom_pseudo_fraction": pseudo_fraction,
        "custom_lambda_scale": lambda_scale,
        "custom_score_mode": score_mode,
        "qvalue_method": "BH",
    }
    return metadata, matches


def _build_comparison(custom_matches: list[dict], fimo_matches: list[dict]) -> dict:
    custom_keys = {
        (
            item["gene"],
            int(item["source_position"]),
            int(item["source_end_position"]),
            item["source_strand_symbol"],
        )
        for item in custom_matches
    }
    fimo_keys = {
        (
            item["gene"],
            int(item["source_position"]),
            int(item["source_end_position"]),
            item["source_strand_symbol"],
        )
        for item in fimo_matches
    }

    overlap = custom_keys & fimo_keys
    only_custom = custom_keys - fimo_keys
    only_fimo = fimo_keys - custom_keys

    overlap_ratio_custom = (len(overlap) / len(custom_keys)) if custom_keys else 0.0
    overlap_ratio_fimo = (len(overlap) / len(fimo_keys)) if fimo_keys else 0.0

    return {
        "custom_total": len(custom_keys),
        "fimo_total": len(fimo_keys),
        "overlap_total": len(overlap),
        "only_custom_total": len(only_custom),
        "only_fimo_total": len(only_fimo),
        "overlap_ratio_custom": overlap_ratio_custom,
        "overlap_ratio_fimo": overlap_ratio_fimo,
    }


@router.post("/", status_code=200)
async def find_motifs(
        file: UploadFile = File(...),
        motif_pattern: str = Form(..., description="IUPAC motif pattern with optional bracket groups"),
        p_value_threshold: float = Form(1e-4),
        score_threshold: Optional[int] = Form(None),
        scan_reverse_complement: bool = Form(True),
        engine: str = Form("custom"),
        fimo_p_value_threshold: float = Form(1e-4),
        custom_calibration_mode: str = Form("legacy"),
):
    filename = file.filename or ""
    if not filename.endswith((".fasta", ".fa", ".fna")):
        raise HTTPException(status_code=400, detail="O formato deve ser FASTA.")

    if engine not in ENGINE_VALUES:
        raise HTTPException(status_code=400, detail=f"engine invalido: {engine}")

    if custom_calibration_mode not in CUSTOM_CALIBRATION_VALUES:
        raise HTTPException(
            status_code=400,
            detail=f"custom_calibration_mode invalido: {custom_calibration_mode}",
        )

    if score_threshold is None and not 0 < p_value_threshold <= 1:
        raise HTTPException(
            status_code=400,
            detail="p_value_threshold deve estar no intervalo (0, 1].",
        )

    if not 0 < fimo_p_value_threshold <= 1:
        raise HTTPException(
            status_code=400,
            detail="fimo_p_value_threshold deve estar no intervalo (0, 1].",
        )

    fasta_bytes = await file.read()
    try:
        fasta_text = fasta_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HTTPException(status_code=400, detail=f"FASTA invalido (utf-8): {exc}") from exc

    records = FastaReader.parse_text(fasta_text)
    if not records:
        raise HTTPException(status_code=400, detail="Nenhum registro FASTA valido encontrado.")

    bg_freqs = BackgroundFrequencies.calculate_from_text(fasta_text)

    custom_metadata: Optional[dict] = None
    custom_matches: Optional[list[dict]] = None
    if engine in {"custom", "both"}:
        try:
            custom_metadata, custom_matches = _run_custom_scan(
                records=records,
                motif_pattern=motif_pattern,
                bg_freqs=bg_freqs,
                p_value_threshold=p_value_threshold,
                score_threshold=score_threshold,
                scan_reverse_complement=scan_reverse_complement,
                custom_calibration_mode=custom_calibration_mode,
            )
        except Exception as exc:
            raise HTTPException(status_code=400, detail=f"Erro no scanner custom: {exc}") from exc

    fimo_matches: Optional[list[dict]] = None
    if engine in {"fimo", "both"}:
        try:
            fimo_matches = FimoRunner.scan(
                motif_pattern=motif_pattern,
                records=records,
                bg_freqs=bg_freqs,
                p_value_threshold=fimo_p_value_threshold,
                scan_reverse_complement=scan_reverse_complement,
            )
        except RuntimeError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        except Exception as exc:
            raise HTTPException(status_code=400, detail=f"Erro na execucao do FIMO: {exc}") from exc

    if engine == "custom":
        return {
            "engine": "custom",
            "metadata": {
                "filename": filename,
                **(custom_metadata or {}),
            },
            "matches": custom_matches or [],
        }

    if engine == "fimo":
        return {
            "engine": "fimo",
            "metadata": {
                "filename": filename,
                "motif_pattern": motif_pattern,
                "motif_length": len(PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern).T),
                "score_threshold": None,
                "cutoff_mode": "pvalue",
                "background_frequencies": bg_freqs,
                "scan_reverse_complement": scan_reverse_complement,
                "fimo_p_value_threshold": fimo_p_value_threshold,
                "fimo_motif_pseudocount": 0.1,
            },
            "matches": fimo_matches or [],
        }

    comparison = _build_comparison(custom_matches or [], fimo_matches or [])
    return {
        "engine": "both",
        "custom": {
            "metadata": {
                "filename": filename,
                **(custom_metadata or {}),
            },
            "matches": custom_matches or [],
        },
        "fimo": {
            "metadata": {
                "filename": filename,
                "motif_pattern": motif_pattern,
                "motif_length": len(PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern).T),
                "score_threshold": None,
                "cutoff_mode": "pvalue",
                "background_frequencies": bg_freqs,
                "scan_reverse_complement": scan_reverse_complement,
                "fimo_p_value_threshold": fimo_p_value_threshold,
                "fimo_motif_pseudocount": 0.1,
            },
            "matches": fimo_matches or [],
        },
        "comparison": comparison,
    }
