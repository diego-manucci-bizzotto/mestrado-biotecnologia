import json
import re
from typing import Any, AsyncGenerator
from fastapi import UploadFile


class FastaRecord:
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence

    def __str__(self):
        return f"{self.header}\n{self.sequence}"


class FastaRecordResult:
    def __init__(self, gene: str, strand_type: str, motif: str, position: int, p_value: str):
        self.gene = gene
        self.strand_type = strand_type
        self.motif = motif
        self.position = position
        self.p_value = p_value


async def process_fasta_file(file: UploadFile, forward_pattern: re.Pattern[str], reverse_pattern: re.Pattern[str]):
    lines = stream_lines(file)

    async for fasta_record in parse_fasta_records(lines):
        output = process_fasta_record(
            fasta_record,
            forward_pattern,
            reverse_pattern
        )

        if output:
            yield output


async def stream_lines(file: UploadFile, read_size: int = SIXTYFOUR_KB) -> AsyncGenerator[str, None]:
    buffer = ""
    while True:
        chunk = await file.read(read_size)
        if not chunk:
            if buffer: yield buffer.strip()
            break

        buffer += chunk.decode('utf-8')
        lines = buffer.split('\n')
        buffer = lines.pop()

        for line in lines:
            yield line.strip()


async def parse_fasta_records(lines: AsyncGenerator[str, None]) -> AsyncGenerator[FastaRecord, None]:
    current_header = None
    current_sequence = []

    async for line in lines:
        if not line: continue

        if is_new_record(line):
            if current_header:
                yield FastaRecord(current_header, "".join(current_sequence))

            current_header = line
            current_sequence = []
        else:
            current_sequence.append(line)

    if current_header:
        yield FastaRecord(current_header, "".join(current_sequence))


def process_fasta_record(fasta_record, forward_pattern: re.Pattern[str], reverse_pattern: re.Pattern[str]) -> str:
    results: list[FastaRecordResult] = []

    sequence_base_proportions = get_sequence_base_proportions(sequence)

    def map_hits(pattern, strand_type):
        for match in pattern.finditer(sequence):
            motif_sequence = match.group(1).upper()
            p_val = calculate_motif_p_value(motif_seq, bg_freqs)
            pos = match.start() + 1
            if is_reverse_strand(fasta_record.header):
                pos = 1000 - pos - len(motif_sequence) + 2

            results.append(
                FastaRecordResult(
                    get_gene_id(fasta_record.header),
                    strand_type,
                    motif_sequence,
                    pos,
                    f"{p_val:.2e}"
                )
            )

    map_hits(forward_pattern, "forward")
    map_hits(reverse_pattern, "reverse")

    if results:
        return "".join(json.dumps(res) + "\n" for res in results)
    return ""


def is_reverse_strand(header: str) -> bool:
    return "[reverse complement]" in header.lower()


def is_new_record(line: str) -> bool:
    return line.startswith(">")

def get_gene_id(header: str) -> str:
    gene_id_match = re.search(r'\((.*?)\)', header)
    return gene_id_match.group(0) if gene_id_match else "unknown"
