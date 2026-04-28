from dataclasses import dataclass

from fastapi import UploadFile

ONE_MEGABYTE = 1024 * 1024


@dataclass
class FastaRecord:
    gene_id: str
    sequence: str
    is_reverse: bool


class FastaReader:
    def __init__(self, file: UploadFile, chunk_size: int = ONE_MEGABYTE):
        self.file = file
        self.chunk_size = chunk_size

    async def parse(self):
        buffer = ""
        current_header = ""
        current_seq = []

        while True:
            chunk = await self.file.read(self.chunk_size)
            if not chunk:
                break

            buffer += chunk.decode('utf-8')
            lines = buffer.split('\n')

            buffer = lines.pop()

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_header:
                        yield self._build_record(current_header, "".join(current_seq))
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line.upper())

        if current_header:
            current_seq.append(buffer.strip().upper())
            yield self._build_record(current_header, "".join(current_seq))

    @classmethod
    def parse_text(cls, text: str) -> list[FastaRecord]:
        records: list[FastaRecord] = []
        current_header = ""
        current_seq: list[str] = []

        for raw_line in text.splitlines():
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header:
                    records.append(cls._build_record(current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_header:
            records.append(cls._build_record(current_header, "".join(current_seq)))

        return records

    @staticmethod
    def _build_record(header: str, sequence: str) -> FastaRecord:
        parts = header.replace(">", "").split()

        gene_id = ""
        for part in parts:
            if part.startswith("(") and part.endswith(")"):
                gene_id = part[1:-1]
                break

        if not gene_id and parts:
            gene_id = parts[0]

        is_reverse = "[reverse complement]" in header.lower()

        return FastaRecord(
            gene_id=gene_id,
            sequence=sequence,
            is_reverse=is_reverse
        )
