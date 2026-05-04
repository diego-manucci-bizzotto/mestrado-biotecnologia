from dataclasses import dataclass
from typing import Dict

import numpy as np


@dataclass(frozen=True)
class IntegerDistribution:
    score_int_min: int
    score_int_max: int
    tail_probabilities: np.ndarray


class IntegerPValueCalculator:
    """
    Exact dynamic-programming p-value calculator for integer PWMs.
    """

    def __init__(self, pwm_integer: np.ndarray, bg_freqs: Dict[str, float]):
        matrix = np.asarray(pwm_integer, dtype=np.int32)
        if matrix.ndim != 2 or matrix.shape[0] != 4:
            raise ValueError("PWM inteiro deve ter shape (4, L).")

        background = np.array(
            [bg_freqs["A"], bg_freqs["C"], bg_freqs["G"], bg_freqs["T"]],
            dtype=np.float64,
        )
        if np.any(background <= 0):
            raise ValueError("As frequencias de fundo devem ser > 0.")
        background = background / background.sum()

        self._matrix = matrix
        self._background = background
        self._distribution = self._build_distribution()

    def _build_distribution(self) -> IntegerDistribution:
        col_min = self._matrix.min(axis=0).astype(np.int64)
        offsets = (self._matrix - col_min[np.newaxis, :]).astype(np.int64)
        score_int_min = int(col_min.sum())
        max_offset = int(offsets.max(axis=0).sum())

        distribution = np.zeros(max_offset + 1, dtype=np.float64)
        distribution[0] = 1.0

        for col_idx in range(self._matrix.shape[1]):
            base_offsets = offsets[:, col_idx]
            next_distribution = np.zeros(
                distribution.size + int(base_offsets.max()),
                dtype=np.float64,
            )

            for base_idx in range(4):
                shift = int(base_offsets[base_idx])
                prob = float(self._background[base_idx])
                if prob <= 0:
                    continue
                next_distribution[shift:shift + distribution.size] += distribution * prob

            distribution = next_distribution

        tail = np.cumsum(distribution[::-1])[::-1]
        tail = np.clip(tail, 0.0, 1.0)

        return IntegerDistribution(
            score_int_min=score_int_min,
            score_int_max=score_int_min + distribution.size - 1,
            tail_probabilities=tail,
        )

    def get_pvalue(self, score: int | float) -> float:
        if not np.isfinite(score):
            raise ValueError("score invalido.")

        score_int = int(np.floor(score))
        if score_int <= self._distribution.score_int_min:
            return 1.0
        if score_int > self._distribution.score_int_max:
            return 0.0

        idx = score_int - self._distribution.score_int_min
        return float(self._distribution.tail_probabilities[idx])

    def get_score_threshold_for_pvalue(self, p_value: float) -> int:
        if not 0 < p_value <= 1:
            raise ValueError("p_value deve estar em (0, 1].")

        tail = self._distribution.tail_probabilities
        candidates = np.flatnonzero(tail <= p_value)
        if candidates.size == 0:
            return int(self._distribution.score_int_max + 1)

        idx = int(candidates[0])
        return int(self._distribution.score_int_min + idx)

    def get_diagnostics(self) -> dict:
        return {
            "method": "integer_exact_dp",
            "score_int_min": self._distribution.score_int_min,
            "score_int_max": self._distribution.score_int_max,
            "score_state_count": int(self._distribution.tail_probabilities.size),
        }
