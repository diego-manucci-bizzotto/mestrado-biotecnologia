from __future__ import annotations

from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile

from app.schemas.analysis import AnalysisDetail, AnalysisSummary, Hit
from app.schemas.motif import Motif, MotifCreate
from app.services.motifs import create_iupac_motif, create_motif, list_motifs
from app.services.scanner import DEFAULT_PWM_P_THRESHOLD, get_analysis, run_analysis

router = APIRouter()
ALLOWED_ANALYSIS_FIELDS = {"file", "motif", "p_value_cutoff"}
DEPRECATED_ANALYSIS_FIELDS = {"scan_mode", "motif_name"}


def _validate_analysis_form_fields(provided_fields: set[str]) -> None:
    deprecated = sorted(provided_fields & DEPRECATED_ANALYSIS_FIELDS)
    unknown = sorted(provided_fields - ALLOWED_ANALYSIS_FIELDS - DEPRECATED_ANALYSIS_FIELDS)
    if deprecated or unknown:
        details: list[str] = []
        if deprecated:
            details.append(f"Deprecated field(s) are not supported: {', '.join(deprecated)}.")
        if unknown:
            details.append(f"Unknown field(s): {', '.join(unknown)}.")
        raise HTTPException(status_code=400, detail=" ".join(details))


@router.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@router.get("/motifs", response_model=list[Motif])
def motifs() -> list[Motif]:
    return list_motifs()


@router.post("/motifs", response_model=Motif)
def add_motif(payload: MotifCreate) -> Motif:
    try:
        return create_motif(payload)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.post("/analyses", response_model=AnalysisDetail)
async def create_analysis(
    request: Request,
    file: UploadFile = File(...),
    motif: str = Form(...),
    p_value_cutoff: float = Form(DEFAULT_PWM_P_THRESHOLD),
) -> AnalysisDetail:
    try:
        form_data = await request.form()
        _validate_analysis_form_fields(set(form_data.keys()))
        motif_model = create_iupac_motif(motif)
        fasta_text = (await file.read()).decode("utf-8")
        return run_analysis(fasta_text, motif_model, p_value_cutoff)
    except UnicodeDecodeError as exc:
        raise HTTPException(status_code=400, detail="FASTA must be UTF-8 text.") from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get("/analyses/{analysis_id}", response_model=AnalysisSummary)
def analysis_summary(analysis_id: str) -> AnalysisSummary:
    try:
        return get_analysis(analysis_id).summary
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@router.get("/analyses/{analysis_id}/results", response_model=list[Hit])
def analysis_results(analysis_id: str, limit: int = 100, offset: int = 0) -> list[Hit]:
    try:
        detail = get_analysis(analysis_id)
        return detail.hits[offset : offset + limit]
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
