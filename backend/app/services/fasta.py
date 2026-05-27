from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class FastaRecord:
    gene_id: str
    description: str
    sequence: str


GENE_ID_RE = re.compile(r"\(([^)]+)\)")


def parse_fasta(text: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_parts: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(_record_from_parts(header, sequence_parts))
            header = line[1:].strip()
            sequence_parts = []
        else:
            sequence_parts.append(line.upper())

    if header is not None:
        records.append(_record_from_parts(header, sequence_parts))

    if not records:
        raise ValueError("No FASTA records were found.")
    return records


def _record_from_parts(header: str, sequence_parts: list[str]) -> FastaRecord:
    sequence = "".join(sequence_parts).replace("U", "T")
    if not sequence:
        raise ValueError(f"FASTA record has no sequence: {header}")

    match = GENE_ID_RE.search(header)
    gene_id = match.group(1) if match else header.split()[0]
    description = header
    if match:
        description = header[match.end() :].strip() or header

    return FastaRecord(gene_id=gene_id, description=description, sequence=sequence)


def background_frequencies(records: list[FastaRecord]) -> tuple[dict[str, float], int, int]:
    counts = {base: 0 for base in "ACGT"}
    total_bases = 0
    for record in records:
        for base in record.sequence:
            if base in counts:
                counts[base] += 1
                total_bases += 1
    if total_bases == 0:
        raise ValueError("No A/C/G/T bases were found in the FASTA.")
    return ({base: counts[base] / total_bases for base in "ACGT"}, total_bases, sum(counts.values()))
