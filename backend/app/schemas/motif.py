from typing import Literal

from pydantic import BaseModel, Field


class Reference(BaseModel):
    label: str
    url: str
    note: str


class Motif(BaseModel):
    id: str
    name: str
    factor: str
    source_type: Literal["iupac", "pwm"]
    iupac_consensus: str | None = None
    pwm_matrix: list[dict[str, float]] | None = None
    organism_scope: str
    references: list[Reference] = Field(default_factory=list)
    limitations: str


class MotifCreate(BaseModel):
    id: str
    name: str
    factor: str
    iupac_consensus: str | None = None
    pwm_matrix: list[dict[str, float]] | None = None
    organism_scope: str = "custom"
    references: list[Reference] = Field(default_factory=list)
    limitations: str = "Custom motif supplied by the user; validate source before publication."
