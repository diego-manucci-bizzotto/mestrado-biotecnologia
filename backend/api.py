from fastapi import APIRouter

from endpoints import motifs, motifs2

api_router = APIRouter()

api_router.include_router(motifs2.router, prefix="/motifs", tags=["Motifs"])