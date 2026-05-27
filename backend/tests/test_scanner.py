from app.services.motifs import create_iupac_motif
from app.services.scanner import run_analysis


def test_pwm_scan_finds_consensus_derived_hit() -> None:
    fasta = ">UPSTREAM1000 (GENE_002) example\nTTTGCCAGGTTT\n"
    detail = run_analysis(fasta, create_iupac_motif("GCCARG", "PacC test"))

    assert detail.summary.observed_hits >= 1
    assert all(hit.score is not None for hit in detail.hits)


def test_pwm_enrichment_expected_hits_uses_cutoff_rate() -> None:
    fasta = ">UPSTREAM1000 (GENE_003) example\nTTTGCCAGGTTT\n"
    cutoff = 1e-4
    detail = run_analysis(fasta, create_iupac_motif("GCCARG", "PacC test"), cutoff)

    assert detail.summary.expected_hit_rate == cutoff
    assert detail.summary.enrichment_test_count == detail.summary.total_windows
    assert detail.summary.enrichment_model == "poisson_pwm_threshold"
    assert detail.summary.expected_hits == detail.summary.total_windows * cutoff
    assert detail.summary.fold_change == detail.summary.observed_hits / detail.summary.expected_hits
    assert 0.0 <= detail.summary.depletion_p_value <= 1.0
    assert 0.0 <= detail.summary.two_sided_p_value <= 1.0
    assert detail.summary.two_sided_p_value == min(
        1.0,
        2.0 * min(detail.summary.enrichment_p_value, detail.summary.depletion_p_value),
    )
