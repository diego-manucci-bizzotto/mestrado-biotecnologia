from __future__ import annotations

import math

_GAMMA_EPS = 1e-14
_GAMMA_MAX_ITERS = 100_000
_GAMMA_FPMIN = 1e-300


def windows_for_length(sequence_length: int, motif_length: int) -> int:
    return max(0, sequence_length - motif_length + 1)


def gene_presence_p_value(window_probability: float, window_count: int) -> float:
    if window_count <= 0:
        return 1.0
    return 1.0 - (1.0 - window_probability) ** window_count


def poisson_survival(observed: int, expected: float) -> float:
    if observed <= 0:
        return 1.0
    if expected <= 0:
        return 0.0
    # P(X >= k) for X ~ Poisson(lambda) equals regularized lower gamma P(k, lambda).
    # Using incomplete-gamma routines avoids underflow from exp(-lambda) when lambda is large.
    return max(0.0, min(1.0, _regularized_gamma_p(float(observed), expected)))


def poisson_cdf(observed: int, expected: float) -> float:
    if observed < 0:
        return 0.0
    if expected <= 0:
        return 1.0
    # P(X <= k) for X ~ Poisson(lambda) equals regularized upper gamma Q(k+1, lambda).
    return max(0.0, min(1.0, _regularized_gamma_q(float(observed + 1), expected)))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    m = len(p_values)
    if m == 0:
        return []
    ordered = sorted(enumerate(p_values), key=lambda item: item[1], reverse=True)
    q_values = [1.0] * m
    previous = 1.0
    for rank_from_end, (index, p_value) in enumerate(ordered):
        rank = m - rank_from_end
        adjusted = min(previous, p_value * m / rank)
        previous = adjusted
        q_values[index] = max(0.0, min(1.0, adjusted))
    return q_values


def _regularized_gamma_p(a: float, x: float) -> float:
    if a <= 0 or x < 0:
        raise ValueError("regularized gamma requires a > 0 and x >= 0.")
    if x == 0:
        return 0.0
    if x < a + 1.0:
        return _regularized_gamma_p_series(a, x)
    return 1.0 - _regularized_gamma_q_continued_fraction(a, x)


def _regularized_gamma_q(a: float, x: float) -> float:
    if a <= 0 or x < 0:
        raise ValueError("regularized gamma requires a > 0 and x >= 0.")
    if x == 0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _regularized_gamma_p_series(a, x)
    return _regularized_gamma_q_continued_fraction(a, x)


def _regularized_gamma_p_series(a: float, x: float) -> float:
    gln = math.lgamma(a)
    ap = a
    term = 1.0 / a
    total = term
    for _ in range(_GAMMA_MAX_ITERS):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * _GAMMA_EPS:
            prefactor = math.exp(-x + a * math.log(x) - gln)
            return total * prefactor
    raise RuntimeError("regularized gamma series did not converge.")


def _regularized_gamma_q_continued_fraction(a: float, x: float) -> float:
    gln = math.lgamma(a)
    b = x + 1.0 - a
    c = 1.0 / _GAMMA_FPMIN
    d = 1.0 / max(abs(b), _GAMMA_FPMIN)
    h = d

    for i in range(1, _GAMMA_MAX_ITERS + 1):
        an = -float(i) * (float(i) - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _GAMMA_FPMIN:
            d = _GAMMA_FPMIN
        c = b + an / c
        if abs(c) < _GAMMA_FPMIN:
            c = _GAMMA_FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _GAMMA_EPS:
            prefactor = math.exp(-x + a * math.log(x) - gln)
            return prefactor * h
    raise RuntimeError("regularized gamma continued fraction did not converge.")
