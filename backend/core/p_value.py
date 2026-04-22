from collections import Counter
from typing import Dict


class PValue:
    def __init__(self, motif_sequence: str, raw_sequence: str):
        self.motif_sequence = motif_sequence
        self.raw_sequence = raw_sequence

        self.background_frequencies = self._get_background_frequencies(self.raw_sequence)

        self.value = self._calculate_p_value(self.motif_sequence, self.background_frequencies)

    def _calculate_p_value(self, motif_sequence: str, background_frequencies: Dict[str, float]) -> float:
        p_value = 1.0
        for base in motif_sequence:
            if base == 'N':
                p_value *= 1.0
            else:
                p_value *= background_frequencies.get(base, 0.0)
        return p_value

    def _get_background_frequencies(self, sequence: str) -> Dict[str, float]:
        total = len(sequence)
        if total == 0:
            return {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

        counts = Counter(sequence)
        return {
            'A': counts.get('A', 0) / total,
            'C': counts.get('C', 0) / total,
            'G': counts.get('G', 0) / total,
            'T': counts.get('T', 0) / total
        }

    def __str__(self):
        return f"P-value: {self.value}"