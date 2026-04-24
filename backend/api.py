from fastapi import APIRouter

from endpoints import motifs

api_router = APIRouter()

api_router.include_router(motifs.router, prefix="/motifs", tags=["Motifs"])