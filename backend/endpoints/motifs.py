import math
from typing import Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from core.bg_freqs import BackgroundFrequencies
from core.fasta_reader import FastaReader, FastaRecord
from core.p_value_calc import PValueCalculator
from core.pwm import PositionWeightMatrix
from core.seq_scanner import SequenceScanner

router = APIRouter()

DEFAULT_PSEUDO_FRACTION = 0.0001
DEFAULT_INITIAL_GRANULARITY = 0.1
DEFAULT_DECREASING_STEP = 10.0
DEFAULT_MAX_REFINEMENT_STEPS = 6
DEFAULT_MAX_SCORE_STATES = 800_000


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
    score_threshold: Optional[float],
    scan_reverse_complement: bool,
) -> tuple[dict, list[dict]]:
    pfm = PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern)
    pwm_log_fwd = PositionWeightMatrix.build_log_odds_pwm(
        pfm,
        bg_freqs,
        pseudo_fraction=DEFAULT_PSEUDO_FRACTION,
    )
    pwm_log_rev = PositionWeightMatrix.get_reverse_complement(pwm_log_fwd)

    pcalc_fwd = PValueCalculator(
        pwm_log_odds=pwm_log_fwd,
        bg_freqs=bg_freqs,
        initial_granularity=DEFAULT_INITIAL_GRANULARITY,
        decreasing_step=DEFAULT_DECREASING_STEP,
        max_refinement_steps=DEFAULT_MAX_REFINEMENT_STEPS,
        max_score_states=DEFAULT_MAX_SCORE_STATES,
    )
    pcalc_rev = PValueCalculator(
        pwm_log_odds=pwm_log_rev,
        bg_freqs=bg_freqs,
        initial_granularity=DEFAULT_INITIAL_GRANULARITY,
        decreasing_step=DEFAULT_DECREASING_STEP,
        max_refinement_steps=DEFAULT_MAX_REFINEMENT_STEPS,
        max_score_states=DEFAULT_MAX_SCORE_STATES,
    )
    motif_length = pwm_log_fwd.shape[1]

    if score_threshold is None:
        cutoff_fwd = float(pcalc_fwd.get_score_threshold_for_pvalue(p_value_threshold))
        cutoff_rev = float(pcalc_rev.get_score_threshold_for_pvalue(p_value_threshold))
        cutoff_mode = "pvalue"
        resolved_score_threshold = None
    else:
        cutoff_fwd = float(score_threshold)
        cutoff_rev = float(score_threshold)
        cutoff_mode = "score"
        resolved_score_threshold = float(score_threshold)

    scan_configs = [("forward", "+", pwm_log_fwd, pcalc_fwd, cutoff_fwd)]
    if scan_reverse_complement:
        scan_configs.append(("reverse", "-", pwm_log_rev, pcalc_rev, cutoff_rev))

    matches: list[dict] = []
    for record in records:
        record_length = len(record.sequence)
        record_orientation = "reverse_complement" if record.is_reverse else "forward"

        for strand_name, strand_symbol, pwm_log, pcalc, cutoff in scan_configs:
            hit_positions, hit_scores = SequenceScanner.scan_log_odds_matrix(
                sequence=record.sequence,
                pwm=pwm_log,
                threshold=cutoff,
            )

            for pos_idx, score_value in zip(hit_positions, hit_scores):
                score_log_odds = float(score_value)
                p_value = float(pcalc.get_pvalue(score_log_odds))

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

                matches.append(
                    {
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
                        "score": round(score_log_odds, 4),
                        "score_log_odds": round(score_log_odds, 4),
                        "p_value": f"{p_value:.2e}",
                        "p_value_numeric": p_value,
                        "-log10_pval": neg_log10_pval,
                    }
                )

    matches.sort(key=lambda item: (item["p_value_numeric"], -item["score_log_odds"]))
    _attach_bh_q_values(matches)

    metadata = {
        "motif_pattern": motif_pattern,
        "motif_length": motif_length,
        "score_threshold": resolved_score_threshold,
        "cutoff_mode": cutoff_mode,
        "background_frequencies": bg_freqs,
        "scan_reverse_complement": scan_reverse_complement,
        "custom_pseudo_fraction": DEFAULT_PSEUDO_FRACTION,
        "qvalue_method": "BH",
        "pvalue_engine": pcalc_fwd.get_diagnostics(),
    }
    return metadata, matches


@router.post("/", status_code=200)
async def find_motifs(
    file: UploadFile = File(...),
    motif_pattern: str = Form(..., description="IUPAC motif pattern with optional bracket groups"),
    p_value_threshold: float = Form(1e-4),
    score_threshold: Optional[float] = Form(None),
    scan_reverse_complement: bool = Form(True),
):
    filename = file.filename or ""
    if not filename.endswith((".fasta", ".fa", ".fna")):
        raise HTTPException(status_code=400, detail="O formato deve ser FASTA.")

    if score_threshold is None and not 0 < p_value_threshold <= 1:
        raise HTTPException(
            status_code=400,
            detail="p_value_threshold deve estar no intervalo (0, 1].",
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

    try:
        metadata, matches = _run_custom_scan(
            records=records,
            motif_pattern=motif_pattern,
            bg_freqs=bg_freqs,
            p_value_threshold=p_value_threshold,
            score_threshold=score_threshold,
            scan_reverse_complement=scan_reverse_complement,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Erro no scanner custom: {exc}") from exc

    return {
        "engine": "custom",
        "metadata": {
            "filename": filename,
            **metadata,
        },
        "matches": matches,
    }
