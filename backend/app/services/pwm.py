from __future__ import annotations

import itertools
import math
import random
from bisect import bisect_left

BASES = "ACGT"
ScoreDistribution = tuple[list[float], list[float]]


def normalize_pwm(matrix: list[dict[str, float]]) -> list[dict[str, float]]:
    normalized: list[dict[str, float]] = []
    for index, column in enumerate(matrix):
        values = {base: float(column.get(base, 0.0)) for base in BASES}
        total = sum(values.values())
        if total <= 0:
            raise ValueError(f"PWM column {index + 1} has no positive probabilities.")
        normalized.append({base: max(values[base] / total, 1e-12) for base in BASES})
    return normalized


def score_window(window: str, matrix: list[dict[str, float]], background: dict[str, float]) -> float:
    score = 0.0
    for base, column in zip(window.upper(), matrix):
        if base not in BASES:
            return float("-inf")
        score += math.log2(column[base] / max(background[base], 1e-12))
    return score


def score_p_value(score: float, matrix: list[dict[str, float]], background: dict[str, float]) -> float:
    return distribution_p_value(score, build_score_distribution(matrix, background))


def build_score_distribution(
    matrix: list[dict[str, float]],
    background: dict[str, float],
    sample_count: int = 100_000,
) -> ScoreDistribution:
    if len(matrix) <= 10:
        return _enumerated_distribution(matrix, background)
    return _sampled_distribution(matrix, background, sample_count=sample_count)


def distribution_p_value(score: float, distribution: ScoreDistribution) -> float:
    scores, suffix_probabilities = distribution
    index = bisect_left(scores, score - 1e-12)
    if index >= len(suffix_probabilities):
        # Report the strongest resolvable lower bound instead of p=0.
        return suffix_probabilities[-1] if suffix_probabilities else 1.0
    return max(0.0, min(1.0, suffix_probabilities[index]))


def _enumerated_distribution(
    matrix: list[dict[str, float]],
    background: dict[str, float],
) -> ScoreDistribution:
    pairs: list[tuple[float, float]] = []
    for kmer in itertools.product(BASES, repeat=len(matrix)):
        sequence = "".join(kmer)
        null_probability = 1.0
        for base in sequence:
            null_probability *= background[base]
        pairs.append((score_window(sequence, matrix, background), null_probability))
    return _suffix_distribution(pairs)


def _sampled_distribution(
    matrix: list[dict[str, float]],
    background: dict[str, float],
    sample_count: int,
) -> ScoreDistribution:
    rng = random.Random(42)
    bases = list(BASES)
    weights = [background[base] for base in bases]
    probability = 1.0 / sample_count
    pairs: list[tuple[float, float]] = []
    for _ in range(sample_count):
        sequence = "".join(rng.choices(bases, weights=weights, k=len(matrix)))
        pairs.append((score_window(sequence, matrix, background), probability))
    return _suffix_distribution(pairs)


def _suffix_distribution(pairs: list[tuple[float, float]]) -> ScoreDistribution:
    pairs.sort(key=lambda item: item[0])
    scores = [score for score, _ in pairs]
    suffix_probabilities = [0.0] * len(pairs)
    running = 0.0
    for index in range(len(pairs) - 1, -1, -1):
        running += pairs[index][1]
        suffix_probabilities[index] = running
    return scores, suffix_probabilities
