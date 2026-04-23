from collections import Counter
from typing import Dict

from fastapi import UploadFile

ONE_MEGABYTE = 1024 * 1024

class BackgroundFrequencies:
    @staticmethod
    async def calculate_from_file(file: UploadFile, chunk_size: int = ONE_MEGABYTE) -> Dict[str, float]:
        total_counts = Counter()

        while contents := await file.read(chunk_size):
            text = contents.decode('utf-8')
            lines = text.split('\n')

            clean_seq = "".join([line.upper() for line in lines if not line.startswith('>')])

            clean_seq = clean_seq.replace('N', '')

            total_counts.update(clean_seq)

        total_bases = sum(total_counts.values())

        if total_bases == 0:
            return {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

        return {
            'A': total_counts.get('A', 0) / total_bases,
            'C': total_counts.get('C', 0) / total_bases,
            'G': total_counts.get('G', 0) / total_bases,
            'T': total_counts.get('T', 0) / total_bases
        }