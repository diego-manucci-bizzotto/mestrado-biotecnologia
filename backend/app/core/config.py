from pydantic import BaseModel


class Settings(BaseModel):
    app_name: str = "Mestrado Biotecnologia Motif Scanner"
    api_prefix: str = "/api"
    cors_origins: list[str] = [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "http://localhost:5174",
        "http://127.0.0.1:5174",
    ]


settings = Settings()
