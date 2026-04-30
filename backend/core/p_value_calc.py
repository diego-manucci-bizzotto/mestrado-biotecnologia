from dataclasses import dataclass
from typing import Dict

import numpy as np


@dataclass(frozen=True)
class RoundDistribution:
    epsilon: float
    max_error: float
    score_int_min: int
    score_int_max: int
    tail_probabilities: np.ndarray
    score_count: int


class PValueCalculator:
    """
    P-value calculator inspired by Touzet & Varre (2007):
    - progressive score discretization (round matrices with decreasing epsilon)
    - dynamic-programming score distribution under zero-order background model
    - stop criterion based on plateau condition:
      P(M_eps, alpha - E) == P(M_eps, alpha)
    """

    def __init__(
        self,
        pwm_log_odds: np.ndarray,
        bg_freqs: Dict[str, float],
        initial_granularity: float = 0.1,
        decreasing_step: float = 10.0,
        max_refinement_steps: int = 6,
        max_score_states: int = 800_000,
    ):
        matrix = np.asarray(pwm_log_odds, dtype=np.float64)
        if matrix.ndim != 2 or matrix.shape[0] != 4:
            raise ValueError("PWM deve ter shape (4, L).")

        background = np.array(
            [bg_freqs["A"], bg_freqs["C"], bg_freqs["G"], bg_freqs["T"]],
            dtype=np.float64,
        )
        if np.any(background <= 0):
            raise ValueError("As frequencias de fundo devem ser > 0.")
        background = background / background.sum()

        if initial_granularity <= 0:
            raise ValueError("initial_granularity deve ser > 0.")
        if decreasing_step <= 1:
            raise ValueError("decreasing_step deve ser > 1.")
        if max_refinement_steps < 1:
            raise ValueError("max_refinement_steps deve ser >= 1.")
        if max_score_states < 1000:
            raise ValueError("max_score_states deve ser >= 1000.")

        self._matrix = matrix
        self._background = background
        self._motif_length = matrix.shape[1]
        self._initial_granularity = float(initial_granularity)
        self._decreasing_step = float(decreasing_step)
        self._max_refinement_steps = int(max_refinement_steps)
        self._max_score_states = int(max_score_states)

        self._levels = self._build_levels()
        self._finest_level = self._levels[-1]
        self._diagnostics = self._build_diagnostics()

    def _build_levels(self) -> list[RoundDistribution]:
        levels: list[RoundDistribution] = []
        epsilon = self._initial_granularity

        for _ in range(self._max_refinement_steps):
            rounded_int = np.floor(self._matrix / epsilon).astype(np.int32)
            col_min = rounded_int.min(axis=0)
            offsets = rounded_int - col_min[np.newaxis, :]
            score_count = int(offsets.max(axis=0).sum()) + 1

            if score_count > self._max_score_states and levels:
                break
            if score_count > self._max_score_states:
                raise ValueError(
                    "Granularidade inicial gerou distribuicao muito grande; "
                    "aumente initial_granularity ou reduza o motivo."
                )

            dist = self._build_distribution_for_int_matrix(rounded_int)
            levels.append(
                RoundDistribution(
                    epsilon=epsilon,
                    max_error=float(self._motif_length * epsilon),
                    score_int_min=dist["score_int_min"],
                    score_int_max=dist["score_int_max"],
                    tail_probabilities=dist["tail_probabilities"],
                    score_count=score_count,
                )
            )

            epsilon /= self._decreasing_step

        if not levels:
            raise ValueError("Nao foi possivel construir distribuicao de score.")
        return levels

    def _build_distribution_for_int_matrix(self, int_matrix: np.ndarray) -> dict:
        col_min = int_matrix.min(axis=0)
        offsets = int_matrix - col_min[np.newaxis, :]
        score_int_min = int(col_min.sum())
        max_offset = int(offsets.max(axis=0).sum())

        distribution = np.zeros(max_offset + 1, dtype=np.float64)
        distribution[0] = 1.0

        for col_idx in range(int_matrix.shape[1]):
            base_offsets = offsets[:, col_idx]
            next_distribution = np.zeros(distribution.size + int(base_offsets.max()), dtype=np.float64)

            for base_idx in range(4):
                shift = int(base_offsets[base_idx])
                prob = float(self._background[base_idx])
                if prob <= 0:
                    continue
                next_distribution[shift:shift + distribution.size] += distribution * prob

            distribution = next_distribution

        tail = np.cumsum(distribution[::-1])[::-1]
        tail = np.clip(tail, 0.0, 1.0)

        return {
            "score_int_min": score_int_min,
            "score_int_max": score_int_min + distribution.size - 1,
            "tail_probabilities": tail,
        }

    @staticmethod
    def _tail_at_or_above(level: RoundDistribution, score: float) -> float:
        scaled = np.nextafter(score / level.epsilon, -np.inf)
        threshold_int = int(np.ceil(scaled))

        if threshold_int <= level.score_int_min:
            return 1.0
        if threshold_int > level.score_int_max:
            return 0.0

        idx = threshold_int - level.score_int_min
        return float(level.tail_probabilities[idx])

    @staticmethod
    def _threshold_with_tail_at_most(level: RoundDistribution, p_value: float) -> float:
        tail = level.tail_probabilities
        candidates = np.flatnonzero(tail <= p_value)
        if candidates.size == 0:
            return float((level.score_int_max + 1) * level.epsilon)
        score_int = level.score_int_min + int(candidates[0])
        return float(score_int * level.epsilon)

    def get_pvalue(self, score: float) -> float:
        if not np.isfinite(score):
            raise ValueError("score invalido.")

        fallback = self._tail_at_or_above(self._finest_level, score)

        for level in self._levels:
            p_at_score = self._tail_at_or_above(level, score)
            p_at_score_minus_error = self._tail_at_or_above(level, score - level.max_error)

            if np.isclose(p_at_score, p_at_score_minus_error, atol=1e-15):
                return p_at_score

            fallback = p_at_score_minus_error

        return float(np.clip(fallback, 0.0, 1.0))

    def get_score_threshold_for_pvalue(self, p_value: float) -> float:
        if not 0 < p_value <= 1:
            raise ValueError("p_value deve estar em (0, 1].")

        fallback = self._threshold_with_tail_at_most(self._finest_level, p_value)

        for level in self._levels:
            candidate = self._threshold_with_tail_at_most(level, p_value)
            p_at_candidate = self._tail_at_or_above(level, candidate)
            p_at_candidate_minus_error = self._tail_at_or_above(level, candidate - level.max_error)

            if np.isclose(p_at_candidate, p_at_candidate_minus_error, atol=1e-15):
                return candidate

            fallback = candidate

        return fallback

    def _build_diagnostics(self) -> dict:
        return {
            "method": "touzet_varre_2007_discretized_dp",
            "initial_granularity": self._initial_granularity,
            "decreasing_step": self._decreasing_step,
            "max_refinement_steps": self._max_refinement_steps,
            "levels_used": len(self._levels),
            "levels": [
                {
                    "epsilon": level.epsilon,
                    "max_error": level.max_error,
                    "score_count": level.score_count,
                }
                for level in self._levels
            ],
        }

    def get_diagnostics(self) -> dict:
        return self._diagnostics
