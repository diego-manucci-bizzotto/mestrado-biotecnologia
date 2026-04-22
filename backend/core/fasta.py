from typing import AsyncIterator

from core.fasta_record import FastaRecord

SIXTY_FOUR_KB = 64 * 1024


class Fasta:
    def __init__(self, file, chunk_size: int = SIXTY_FOUR_KB):
        self.file = file
        self.chunk_size = chunk_size
        self.buffer = ""

    async def __aiter__(self) -> AsyncIterator[FastaRecord]:
        current_header = None
        current_sequence = []

        async for line in self._read_lines():
            if self._is_new_record(line):
                if current_header:
                    yield FastaRecord(current_header, "".join(current_sequence))

                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header:
            yield FastaRecord(current_header, "".join(current_sequence))

    async def _read_lines(self) -> AsyncIterator[str]:
        while True:
            chunk = await self.file.read(self.chunk_size)
            if not chunk:
                if self.buffer:
                    yield self.buffer.strip()
                break

            if isinstance(chunk, bytes):
                chunk = chunk.decode('utf-8')

            self.buffer += chunk
            lines = self.buffer.split('\n')
            self.buffer = lines.pop()

            for line in lines:
                clean_line = line.strip()
                if clean_line:
                    yield clean_line

    def _is_new_record(self, line: str) -> bool:
        return line.startswith(">")

    def __repr__(self):
        return f"Fasta(file={self.file.filename}, chunk_size={self.chunk_size})"
