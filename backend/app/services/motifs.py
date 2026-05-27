from __future__ import annotations

import re

from app.schemas.motif import Motif, MotifCreate
from app.services.iupac import iupac_to_pwm, normalize_iupac
from app.services.pwm import normalize_pwm


MOTIFS: dict[str, Motif] = {}


def list_motifs() -> list[Motif]:
    return sorted(MOTIFS.values(), key=lambda motif: motif.id)


def get_motif(motif_id: str) -> Motif:
    try:
        return MOTIFS[motif_id]
    except KeyError as exc:
        raise ValueError(f"Unknown motif_id: {motif_id}") from exc


def create_motif(payload: MotifCreate) -> Motif:
    if bool(payload.iupac_consensus) == bool(payload.pwm_matrix):
        raise ValueError("Provide exactly one of iupac_consensus or pwm_matrix.")

    source_type = "iupac" if payload.iupac_consensus else "pwm"
    motif = Motif(
        id=payload.id,
        name=payload.name,
        factor=payload.factor,
        source_type=source_type,
        iupac_consensus=normalize_iupac(payload.iupac_consensus) if payload.iupac_consensus else None,
        pwm_matrix=normalize_pwm(payload.pwm_matrix) if payload.pwm_matrix else None,
        organism_scope=payload.organism_scope,
        references=payload.references,
        limitations=payload.limitations,
    )
    MOTIFS[motif.id] = motif
    return motif


def create_iupac_motif(consensus: str, name: str | None = None) -> Motif:
    normalized = normalize_iupac(consensus)
    if not normalized:
        raise ValueError("Motif IUPAC is required.")

    # Validates every IUPAC code before the scan starts.
    iupac_to_pwm(normalized)

    clean_id = re.sub(r"[^A-Z0-9]+", "-", normalized).strip("-").lower()
    return Motif(
        id=f"custom-{clean_id}",
        name=name.strip() if name and name.strip() else f"Motif IUPAC {normalized}",
        factor=name.strip() if name and name.strip() else "Custom IUPAC motif",
        source_type="iupac",
        iupac_consensus=normalized,
        organism_scope="user-supplied FASTA",
        references=[],
        limitations=(
            "Motif informado pelo usuario. A relevancia biologica depende da fonte do consenso "
            "e deve ser defendida com p-value, q-value, conservacao e contexto funcional."
        ),
    )


def motif_pwm(motif: Motif) -> list[dict[str, float]]:
    if motif.pwm_matrix:
        return normalize_pwm(motif.pwm_matrix)
    if motif.iupac_consensus:
        return iupac_to_pwm(motif.iupac_consensus)
    raise ValueError("Motif has no usable consensus or PWM.")
