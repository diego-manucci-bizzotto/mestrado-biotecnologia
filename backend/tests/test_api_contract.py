from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def _fasta_bytes() -> bytes:
    return b">UPSTREAM1000 (GENE_001) example\nTTTGCCAGGTTT\n"


def test_create_analysis_accepts_pwm_only_contract() -> None:
    response = client.post(
        "/api/analyses",
        data={"motif": "GCCARG", "p_value_cutoff": "1e-4"},
        files={"file": ("test.fasta", _fasta_bytes(), "text/plain")},
    )

    assert response.status_code == 200
    payload = response.json()
    summary = payload["summary"]
    assert "scan_mode" not in summary
    assert "expected_hit_rate" in summary
    assert "enrichment_test_count" in summary
    assert "enrichment_model" in summary
    assert "depletion_p_value" in summary
    assert "two_sided_p_value" in summary
    assert "fold_change" in summary
    assert summary["enrichment_model"] == "poisson_pwm_threshold"


def test_create_analysis_rejects_deprecated_fields() -> None:
    response = client.post(
        "/api/analyses",
        data={
            "motif": "GCCARG",
            "p_value_cutoff": "1e-4",
            "scan_mode": "pwm",
            "motif_name": "Deprecated",
        },
        files={"file": ("test.fasta", _fasta_bytes(), "text/plain")},
    )

    assert response.status_code == 400
    detail = response.json()["detail"]
    assert "Deprecated field(s) are not supported" in detail
    assert "scan_mode" in detail
    assert "motif_name" in detail


def test_create_analysis_rejects_unknown_fields() -> None:
    response = client.post(
        "/api/analyses",
        data={"motif": "GCCARG", "p_value_cutoff": "1e-4", "extra_field": "x"},
        files={"file": ("test.fasta", _fasta_bytes(), "text/plain")},
    )

    assert response.status_code == 400
    detail = response.json()["detail"]
    assert "Unknown field(s): extra_field." in detail
