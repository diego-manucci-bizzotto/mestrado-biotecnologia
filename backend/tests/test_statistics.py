import math

from app.services.statistics import poisson_cdf, poisson_survival


def _poisson_tail_direct(observed: int, expected: float) -> float:
    # Direct definition used only for tiny-parameter regression checks.
    cumulative = 0.0
    for k in range(observed):
        cumulative += math.exp(-expected) * (expected**k) / math.factorial(k)
    return 1.0 - cumulative


def test_poisson_survival_matches_direct_tail_for_small_values() -> None:
    expected = 0.1
    assert abs(poisson_survival(1, expected) - _poisson_tail_direct(1, expected)) < 1e-12
    assert abs(poisson_survival(2, expected) - _poisson_tail_direct(2, expected)) < 1e-12


def test_poisson_cdf_and_survival_are_complementary() -> None:
    expected = 17.25
    observed = 9
    assert abs(poisson_survival(observed, expected) - (1.0 - poisson_cdf(observed - 1, expected))) < 1e-12


def test_poisson_survival_large_expected_large_observed_stays_small() -> None:
    p_value = poisson_survival(48_925, 1_712.8608)
    assert 0.0 <= p_value <= 1.0
    assert p_value < 1e-30


def test_poisson_survival_large_expected_small_observed_stays_near_one() -> None:
    p_value = poisson_survival(1, 1_000.0)
    assert 0.0 <= p_value <= 1.0
    assert p_value > 0.999999999


def test_poisson_cdf_degenerate_expected_zero() -> None:
    assert poisson_cdf(0, 0.0) == 1.0
    assert poisson_cdf(10, 0.0) == 1.0
