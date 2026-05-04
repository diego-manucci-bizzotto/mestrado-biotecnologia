import math
from dataclasses import dataclass
from typing import Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from core.bg_freqs import BackgroundFrequencies
from core.fasta_reader import FastaReader, FastaRecord
from core.integer_p_value import IntegerPValueCalculator
from core.p_value_calc import PValueCalculator
from core.pwm import PositionWeightMatrix
from core.seq_scanner import SequenceScanner

router = APIRouter()

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")

PSEUDO_FRACTION = 0.0001
INITIAL_GRANULARITY = 0.1
DECREASING_STEP = 10.0
MAX_REFINEMENT_STEPS = 6
MAX_SCORE_STATES = 800_000
INTEGER_LAMBDA_SCALE = 100.0

FIELD_DEFINITIONS = {
    "record_orientation": "Orientation of the uploaded FASTA record: forward or reverse_complement.",
    "strand": "Scanner strand relative to the scanned sequence: forward or reverse.",
    "source_strand": "Biological strand after mapping the hit back to source coordinates.",
    "score": "PWM log2-odds score. Larger means a stronger motif match.",
    "p_value": "Probability of observing this score or higher under the FASTA background model.",
    "bed_start_0based": "BED-compatible 0-based start.",
    "bed_end_0based": "BED-compatible end (half-open).",
}

STATISTICAL_NOTES = {
    "background_model": "A/C/G/T frequencies are estimated from the uploaded FASTA.",
    "thresholding": "Either threshold by p-value (derived score cutoff) or by direct score.",
    "article_connection": "P-values use discretized score distributions with progressive refinement.",
    "boundary_reconciliation": "Boundary hits are reconciled against an internal integer-DP acceptance check.",
}


@dataclass(frozen=True)
class ModelBundle:
    pfm: object
    motif_length: int
    pwm_forward_float: object
    pwm_reverse_float: object
    pwm_forward_int: object
    pwm_reverse_int: object
    pcalc_forward_float: PValueCalculator
    pcalc_reverse_float: PValueCalculator
    pcalc_forward_int: IntegerPValueCalculator
    pcalc_reverse_int: IntegerPValueCalculator


@dataclass(frozen=True)
class ScanConfig:
    strand_name: str
    strand_symbol: str
    pwm_float: object
    cutoff_float: float
    pwm_int: object
    pcalc_int: IntegerPValueCalculator


@router.post("/", status_code=200)
async def find_motifs(
    file: UploadFile = File(...),
    motif_pattern: str = Form(..., description="IUPAC motif pattern with optional bracket groups"),
    p_value_threshold: float = Form(
        1e-4,
        description="Maximum p-value to report when score_threshold is not supplied.",
    ),
    score_threshold: Optional[float] = Form(
        None,
        description="Optional direct PWM log2-odds score cutoff. If supplied, p_value_threshold is not used as the cutoff.",
    ),
    scan_reverse_complement: bool = Form(
        True,
        description="Scan both forward and reverse-complement motif orientations.",
    ),
):
    filename = file.filename or ""
    _validate_request(filename, p_value_threshold, score_threshold)

    fasta_text = await _read_fasta_text(file)
    records = FastaReader.parse_text(fasta_text)
    if not records:
        raise HTTPException(status_code=400, detail="Nenhum registro FASTA valido encontrado.")

    background_frequencies = BackgroundFrequencies.calculate_from_text(fasta_text)

    try:
        metadata, matches = _scan_records(
            records=records,
            motif_pattern=motif_pattern.strip().upper(),
            background_frequencies=background_frequencies,
            p_value_threshold=p_value_threshold,
            score_threshold=score_threshold,
            scan_reverse_complement=scan_reverse_complement,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Erro no scanner de motivos: {exc}") from exc

    return {
        "engine": "custom",
        "metadata": {
            "filename": filename,
            **metadata,
        },
        "matches": matches,
    }


def _validate_request(
    filename: str,
    p_value_threshold: float,
    score_threshold: Optional[float],
) -> None:
    if not filename.lower().endswith(FASTA_EXTENSIONS):
        raise HTTPException(status_code=400, detail="O formato deve ser FASTA (.fasta, .fa ou .fna).")

    if score_threshold is None and not 0 < p_value_threshold <= 1:
        raise HTTPException(status_code=400, detail="p_value_threshold deve estar no intervalo (0, 1].")


async def _read_fasta_text(file: UploadFile) -> str:
    fasta_bytes = await file.read()
    try:
        return fasta_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HTTPException(status_code=400, detail=f"FASTA invalido (utf-8): {exc}") from exc


def _scan_records(
    *,
    records: list[FastaRecord],
    motif_pattern: str,
    background_frequencies: dict[str, float],
    p_value_threshold: float,
    score_threshold: Optional[float],
    scan_reverse_complement: bool,
) -> tuple[dict, list[dict]]:
    return _scan_records_custom(
        records=records,
        motif_pattern=motif_pattern,
        background_frequencies=background_frequencies,
        p_value_threshold=p_value_threshold,
        score_threshold=score_threshold,
        scan_reverse_complement=scan_reverse_complement,
    )


def _scan_records_custom(
    *,
    records: list[FastaRecord],
    motif_pattern: str,
    background_frequencies: dict[str, float],
    p_value_threshold: float,
    score_threshold: Optional[float],
    scan_reverse_complement: bool,
) -> tuple[dict, list[dict]]:
    # Step 1: Build motif models and p-value engines in float + integer spaces.
    models = _build_models(motif_pattern=motif_pattern, background_frequencies=background_frequencies)

    # Step 2: Resolve strand-specific cutoffs and scan configuration.
    scan_configs, cutoff_mode, reported_score_threshold = _build_scan_configs(
        models=models,
        p_value_threshold=p_value_threshold,
        score_threshold=score_threshold,
        scan_reverse_complement=scan_reverse_complement,
    )

    # Step 3: In p-value mode, build an integer-DP acceptance set for boundary reconciliation.
    boundary_accept_keys: Optional[set[tuple[str, str, int, int]]] = None
    if score_threshold is None:
        boundary_accept_keys = _collect_boundary_accept_keys(
            records=records,
            models=models,
            p_value_threshold=p_value_threshold,
            scan_reverse_complement=scan_reverse_complement,
        )

    # Step 4: Scan all records and project hits back to source-space coordinates.
    matches = _scan_all_records(
        records=records,
        motif_length=models.motif_length,
        scan_configs=scan_configs,
        p_value_threshold=p_value_threshold,
        score_threshold=score_threshold,
        boundary_accept_keys=boundary_accept_keys,
    )

    # Step 5: Sort and build response metadata.
    matches.sort(key=lambda item: (item["p_value_numeric"], -_score_for_rank(item)))
    metadata = _build_metadata(
        motif_pattern=motif_pattern,
        motif_length=models.motif_length,
        cutoff_mode=cutoff_mode,
        reported_score_threshold=reported_score_threshold,
        background_frequencies=background_frequencies,
        scan_reverse_complement=scan_reverse_complement,
        pvalue_engine=models.pcalc_forward_int.get_diagnostics(),
    )
    return metadata, matches


def _build_models(
    *,
    motif_pattern: str,
    background_frequencies: dict[str, float],
) -> ModelBundle:
    pfm = PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern)
    pwm_forward_float = PositionWeightMatrix.build_log_odds_pwm(
        pfm,
        background_frequencies,
        pseudo_fraction=PSEUDO_FRACTION,
    )
    pwm_reverse_float = PositionWeightMatrix.get_reverse_complement(pwm_forward_float)

    pwm_forward_int = PositionWeightMatrix.build_integer_pwm(
        pfm,
        background_frequencies,
        pseudo_fraction=PSEUDO_FRACTION,
        lambda_scale=INTEGER_LAMBDA_SCALE,
    )
    pwm_reverse_int = PositionWeightMatrix.get_reverse_complement(pwm_forward_int)

    pcalc_forward_float = _new_pvalue_calculator(pwm_forward_float, background_frequencies)
    pcalc_reverse_float = _new_pvalue_calculator(pwm_reverse_float, background_frequencies)
    pcalc_forward_int = _new_integer_pvalue_calculator(pwm_forward_int, background_frequencies)
    pcalc_reverse_int = _new_integer_pvalue_calculator(pwm_reverse_int, background_frequencies)

    return ModelBundle(
        pfm=pfm,
        motif_length=int(pwm_forward_float.shape[1]),
        pwm_forward_float=pwm_forward_float,
        pwm_reverse_float=pwm_reverse_float,
        pwm_forward_int=pwm_forward_int,
        pwm_reverse_int=pwm_reverse_int,
        pcalc_forward_float=pcalc_forward_float,
        pcalc_reverse_float=pcalc_reverse_float,
        pcalc_forward_int=pcalc_forward_int,
        pcalc_reverse_int=pcalc_reverse_int,
    )


def _build_scan_configs(
    *,
    models: ModelBundle,
    p_value_threshold: float,
    score_threshold: Optional[float],
    scan_reverse_complement: bool,
) -> tuple[list[ScanConfig], str, Optional[float]]:
    if score_threshold is None:
        cutoff_forward = float(models.pcalc_forward_float.get_score_threshold_for_pvalue(p_value_threshold))
        cutoff_reverse = float(models.pcalc_reverse_float.get_score_threshold_for_pvalue(p_value_threshold))
        cutoff_mode = "pvalue"
        reported_score_threshold = None
    else:
        cutoff_forward = float(score_threshold)
        cutoff_reverse = float(score_threshold)
        cutoff_mode = "score"
        reported_score_threshold = float(score_threshold)

    scan_configs = [
        ScanConfig(
            strand_name="forward",
            strand_symbol="+",
            pwm_float=models.pwm_forward_float,
            cutoff_float=cutoff_forward,
            pwm_int=models.pwm_forward_int,
            pcalc_int=models.pcalc_forward_int,
        )
    ]
    if scan_reverse_complement:
        scan_configs.append(
            ScanConfig(
                strand_name="reverse",
                strand_symbol="-",
                pwm_float=models.pwm_reverse_float,
                cutoff_float=cutoff_reverse,
                pwm_int=models.pwm_reverse_int,
                pcalc_int=models.pcalc_reverse_int,
            )
        )

    return scan_configs, cutoff_mode, reported_score_threshold


def _collect_boundary_accept_keys(
    *,
    records: list[FastaRecord],
    models: ModelBundle,
    p_value_threshold: float,
    scan_reverse_complement: bool,
) -> set[tuple[str, str, int, int]]:
    cutoff_forward = int(models.pcalc_forward_int.get_score_threshold_for_pvalue(p_value_threshold))
    cutoff_reverse = int(models.pcalc_reverse_int.get_score_threshold_for_pvalue(p_value_threshold))

    int_scan_configs = [
        ("+", models.pwm_forward_int, models.pcalc_forward_int, cutoff_forward),
    ]
    if scan_reverse_complement:
        int_scan_configs.append(("-", models.pwm_reverse_int, models.pcalc_reverse_int, cutoff_reverse))

    accept_keys: set[tuple[str, str, int, int]] = set()
    for record in records:
        for strand_symbol, pwm_int, pcalc_int, cutoff in int_scan_configs:
            hit_positions, hit_scores = SequenceScanner.scan_integer_matrix(
                sequence=record.sequence,
                pwm=pwm_int,
                threshold=cutoff,
            )
            for position_idx, score_value in zip(hit_positions, hit_scores):
                p_value = float(pcalc_int.get_pvalue(int(score_value)))
                if p_value > p_value_threshold:
                    continue
                scan_start = int(position_idx) + 1
                scan_end = scan_start + models.motif_length - 1
                accept_keys.add((record.gene_id, strand_symbol, scan_start, scan_end))

    return accept_keys


def _scan_all_records(
    *,
    records: list[FastaRecord],
    motif_length: int,
    scan_configs: list[ScanConfig],
    p_value_threshold: float,
    score_threshold: Optional[float],
    boundary_accept_keys: Optional[set[tuple[str, str, int, int]]],
) -> list[dict]:
    matches: list[dict] = []
    for record in records:
        matches.extend(
            _scan_one_record_custom(
                record=record,
                motif_length=motif_length,
                scan_configs=scan_configs,
                p_value_threshold=p_value_threshold,
                score_threshold=score_threshold,
                boundary_accept_keys=boundary_accept_keys,
            )
        )
    return matches


def _build_metadata(
    *,
    motif_pattern: str,
    motif_length: int,
    cutoff_mode: str,
    reported_score_threshold: Optional[float],
    background_frequencies: dict[str, float],
    scan_reverse_complement: bool,
    pvalue_engine: dict,
) -> dict:
    return {
        "motif_pattern": motif_pattern,
        "motif_length": motif_length,
        "score_threshold": reported_score_threshold,
        "cutoff_mode": cutoff_mode,
        "background_frequencies": background_frequencies,
        "scan_reverse_complement": scan_reverse_complement,
        "custom_pseudo_fraction": PSEUDO_FRACTION,
        "pvalue_engine": pvalue_engine,
        "field_definitions": FIELD_DEFINITIONS,
        "statistical_notes": STATISTICAL_NOTES,
    }


def _new_pvalue_calculator(
    pwm_log_odds,
    background_frequencies: dict[str, float],
) -> PValueCalculator:
    return PValueCalculator(
        pwm_log_odds=pwm_log_odds,
        bg_freqs=background_frequencies,
        initial_granularity=INITIAL_GRANULARITY,
        decreasing_step=DECREASING_STEP,
        max_refinement_steps=MAX_REFINEMENT_STEPS,
        max_score_states=MAX_SCORE_STATES,
    )


def _new_integer_pvalue_calculator(
    pwm_integer,
    background_frequencies: dict[str, float],
) -> IntegerPValueCalculator:
    return IntegerPValueCalculator(
        pwm_integer=pwm_integer,
        bg_freqs=background_frequencies,
    )


def _scan_one_record_custom(
    *,
    record: FastaRecord,
    motif_length: int,
    scan_configs: list[ScanConfig],
    p_value_threshold: float,
    score_threshold: Optional[float],
    boundary_accept_keys: Optional[set[tuple[str, str, int, int]]],
) -> list[dict]:
    record_orientation = "reverse_complement" if record.is_reverse else "forward"
    record_length = len(record.sequence)
    record_matches: list[dict] = []

    for cfg in scan_configs:
        hit_positions, hit_scores = SequenceScanner.scan_log_odds_matrix(
            sequence=record.sequence,
            pwm=cfg.pwm_float,
            threshold=cfg.cutoff_float,
        )

        for position_idx, score_value in zip(hit_positions, hit_scores):
            score_log_odds = float(score_value)
            scan_start = int(position_idx) + 1
            scan_end = scan_start + motif_length - 1
            candidate_key = (record.gene_id, cfg.strand_symbol, scan_start, scan_end)
            sequence_window = record.sequence[position_idx:position_idx + motif_length]

            # Report p-values from integer DP (PWMScan-style computation space).
            score_integer = _score_sequence_with_integer_pwm(sequence_window, cfg.pwm_int)
            p_value = float(cfg.pcalc_int.get_pvalue(score_integer))

            if score_threshold is None and p_value > p_value_threshold:
                if boundary_accept_keys is None or candidate_key not in boundary_accept_keys:
                    continue

            matched_sequence = (
                PositionWeightMatrix.reverse_complement_sequence(sequence_window)
                if cfg.strand_symbol == "-"
                else sequence_window
            )

            source_start, source_end, source_strand_symbol = _to_source_coordinates(
                record_length=record_length,
                scan_start=scan_start,
                scan_end=scan_end,
                record_is_reverse=record.is_reverse,
                strand_symbol=cfg.strand_symbol,
            )
            source_strand = "reverse" if source_strand_symbol == "-" else "forward"

            record_matches.append(
                {
                    "gene": record.gene_id,
                    "record_orientation": record_orientation,
                    "strand": cfg.strand_name,
                    "strand_symbol": cfg.strand_symbol,
                    "source_strand": source_strand,
                    "source_strand_symbol": source_strand_symbol,
                    "matched_sequence": matched_sequence,
                    "position": scan_start,
                    "end_position": scan_end,
                    "bed_start_0based": scan_start - 1,
                    "bed_end_0based": scan_end,
                    "source_position": source_start,
                    "source_end_position": source_end,
                    "score": round(score_log_odds, 4),
                    "score_log_odds": round(score_log_odds, 4),
                    "score_integer": score_integer,
                    "p_value": f"{p_value:.2e}",
                    "p_value_numeric": p_value,
                    "-log10_pval": round(-math.log10(p_value), 2) if p_value > 0 else None,
                }
            )

    return record_matches


def _to_source_coordinates(
    *,
    record_length: int,
    scan_start: int,
    scan_end: int,
    record_is_reverse: bool,
    strand_symbol: str,
) -> tuple[int, int, str]:
    if not record_is_reverse:
        return scan_start, scan_end, strand_symbol

    source_start = record_length - scan_end + 1
    source_end = record_length - scan_start + 1
    source_strand_symbol = "-" if strand_symbol == "+" else "+"
    return source_start, source_end, source_strand_symbol


def _score_sequence_with_integer_pwm(sequence: str, pwm_integer) -> int:
    score = 0
    for idx, base in enumerate(sequence):
        base_idx = PositionWeightMatrix.BASE_TO_ID.get(base)
        if base_idx is None:
            raise ValueError(f"Base invalida '{base}' ao calcular score inteiro.")
        score += int(pwm_integer[base_idx, idx])
    return score


def _score_for_rank(match: dict) -> float:
    return float(match["score_log_odds"])
