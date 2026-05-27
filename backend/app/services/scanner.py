from __future__ import annotations

import uuid

from app.schemas.analysis import AnalysisDetail, AnalysisSummary, BackgroundFrequencies, Hit
from app.schemas.motif import Motif
from app.services.fasta import FastaRecord, background_frequencies, parse_fasta
from app.services.iupac import reverse_complement
from app.services.motifs import motif_pwm
from app.services.pwm import build_score_distribution, distribution_p_value, score_window
from app.services.statistics import (
    benjamini_hochberg,
    gene_presence_p_value,
    poisson_cdf,
    poisson_survival,
    windows_for_length,
)

ANALYSES: dict[str, AnalysisDetail] = {}
DEFAULT_PWM_P_THRESHOLD = 0.001


def run_analysis(
    fasta_text: str,
    motif: Motif,
    pwm_p_value_cutoff: float = DEFAULT_PWM_P_THRESHOLD,
) -> AnalysisDetail:
    if not (0.0 < pwm_p_value_cutoff <= 1.0):
        raise ValueError("p_value_cutoff must be > 0 and <= 1.")

    records = parse_fasta(fasta_text)
    background, total_bases, _ = background_frequencies(records)
    matrix = motif_pwm(motif)
    motif_length = len(matrix)
    total_windows = sum(windows_for_length(len(record.sequence), motif_length) for record in records) * 2

    raw_hits = _scan_pwm(records, matrix, background, pwm_p_value_cutoff)

    q_values = benjamini_hochberg([hit.window_p_value for hit in raw_hits])
    hits_with_q = [
        hit.model_copy(update={"q_value": q_values[index]})
        for index, hit in enumerate(raw_hits)
    ]
    hits = [
        hit.model_copy(update={"id": index + 1})
        for index, hit in enumerate(
            sorted(hits_with_q, key=lambda item: (item.window_p_value, item.q_value))
        )
    ]

    expected_hit_rate = pwm_p_value_cutoff
    enrichment_test_count = total_windows
    expected_hits = enrichment_test_count * expected_hit_rate
    observed_hits = len(hits)
    enrichment_p_value = poisson_survival(observed_hits, expected_hits)
    enrichment_q_value = benjamini_hochberg([enrichment_p_value])[0]
    depletion_p_value = poisson_cdf(observed_hits, expected_hits)
    two_sided_p_value = min(1.0, 2.0 * min(enrichment_p_value, depletion_p_value))
    fold_change = observed_hits / expected_hits if expected_hits > 0 else float("inf")

    analysis_id = str(uuid.uuid4())
    summary = AnalysisSummary(
        id=analysis_id,
        motif_id=motif.id,
        motif_name=motif.name,
        sequence_count=len(records),
        total_bases=total_bases,
        total_windows=total_windows,
        observed_hits=observed_hits,
        expected_hits=expected_hits,
        expected_hit_rate=expected_hit_rate,
        enrichment_test_count=enrichment_test_count,
        enrichment_model="poisson_pwm_threshold",
        enrichment_p_value=enrichment_p_value,
        enrichment_q_value=enrichment_q_value,
        depletion_p_value=depletion_p_value,
        two_sided_p_value=two_sided_p_value,
        fold_change=fold_change,
        background=BackgroundFrequencies(**background),
        status="completed",
    )
    detail = AnalysisDetail(summary=summary, hits=hits)
    ANALYSES[analysis_id] = detail
    return detail


def get_analysis(analysis_id: str) -> AnalysisDetail:
    try:
        return ANALYSES[analysis_id]
    except KeyError as exc:
        raise ValueError(f"Unknown analysis_id: {analysis_id}") from exc


def _scan_pwm(
    records: list[FastaRecord],
    matrix: list[dict[str, float]],
    background: dict[str, float],
    p_value_cutoff: float,
) -> list[Hit]:
    motif_length = len(matrix)
    distribution = build_score_distribution(matrix, background)
    hits: list[Hit] = []
    for record in records:
        for strand, sequence in (("+", record.sequence), ("-", reverse_complement(record.sequence))):
            window_count = windows_for_length(len(sequence), motif_length)
            for start in range(window_count):
                window = sequence[start : start + motif_length].upper()
                score = score_window(window, matrix, background)
                if score == float("-inf"):
                    continue
                p_value = distribution_p_value(score, distribution)
                if p_value <= p_value_cutoff:
                    hits.append(
                        Hit(
                            id=0,
                            gene_id=record.gene_id,
                            description=record.description,
                            strand=strand,
                            start=start + 1,
                            end=start + motif_length,
                            matched_sequence=window,
                            score=score,
                            window_p_value=p_value,
                            gene_p_value=gene_presence_p_value(p_value, window_count),
                            q_value=1.0,
                        )
                    )
    return hits
