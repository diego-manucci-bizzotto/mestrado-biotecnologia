from fastapi import FastAPI, APIRouter
from starlette.middleware.cors import CORSMiddleware

from endpoints import motifs

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

api_router = APIRouter()

api_router.include_router(motifs.router, prefix="/motifs", tags=["Motifs"])

app.include_router(api_router, prefix="/api")
