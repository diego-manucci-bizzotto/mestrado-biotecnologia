from typing import Literal

from pydantic import BaseModel


class BackgroundFrequencies(BaseModel):
    A: float
    C: float
    G: float
    T: float


class AnalysisSummary(BaseModel):
    id: str
    motif_id: str
    motif_name: str
    sequence_count: int
    total_bases: int
    total_windows: int
    observed_hits: int
    expected_hits: float
    expected_hit_rate: float
    enrichment_test_count: int
    enrichment_model: Literal["poisson_pwm_threshold"]
    enrichment_p_value: float
    enrichment_q_value: float
    depletion_p_value: float
    two_sided_p_value: float
    fold_change: float
    background: BackgroundFrequencies
    status: Literal["completed"]


class Hit(BaseModel):
    id: int
    gene_id: str
    description: str
    strand: Literal["+", "-"]
    start: int
    end: int
    matched_sequence: str
    score: float | None
    window_p_value: float
    gene_p_value: float
    q_value: float


class AnalysisDetail(BaseModel):
    summary: AnalysisSummary
    hits: list[Hit]
