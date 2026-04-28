from collections import Counter
from typing import Dict

from fastapi import UploadFile

ONE_MEGABYTE = 1024 * 1024

class BackgroundFrequencies:
    @staticmethod
    def calculate_from_text(text: str) -> Dict[str, float]:
        total_counts = Counter()

        for raw_line in text.splitlines():
            line = raw_line.strip().upper()
            if not line or line.startswith(">"):
                continue
            total_counts.update(base for base in line if base in "ACGT")
        total_bases = sum(total_counts.values())

        if total_bases == 0:
            return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

        return {
            "A": total_counts.get("A", 0) / total_bases,
            "C": total_counts.get("C", 0) / total_bases,
            "G": total_counts.get("G", 0) / total_bases,
            "T": total_counts.get("T", 0) / total_bases,
        }

    @staticmethod
    async def calculate_from_file(file: UploadFile, chunk_size: int = ONE_MEGABYTE) -> Dict[str, float]:
        total_counts = Counter()
        buffer = ""
        in_header = False

        while contents := await file.read(chunk_size):
            text = buffer + contents.decode('utf-8')
            lines = text.split('\n')
            buffer = lines.pop()

            for line in lines:
                line = line.strip().upper()
                if not line:
                    continue

                if line.startswith('>'):
                    in_header = True
                    continue

                if in_header:
                    in_header = False

                total_counts.update(base for base in line if base in "ACGT")

        if buffer:
            line = buffer.strip().upper()
            if line and not line.startswith('>') and not in_header:
                total_counts.update(base for base in line if base in "ACGT")

        total_bases = sum(total_counts.values())

        if total_bases == 0:
            return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

        return {
            "A": total_counts.get("A", 0) / total_bases,
            "C": total_counts.get("C", 0) / total_bases,
            "G": total_counts.get("G", 0) / total_bases,
            "T": total_counts.get("T", 0) / total_bases
        }
