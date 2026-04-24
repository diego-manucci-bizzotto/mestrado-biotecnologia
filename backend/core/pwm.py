import numpy as np
from typing import Dict

class PositionWeightMatrix:

    BASE_TO_ID = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }

    @classmethod
    def parse_motif_pattern_to_pfm(cls, motif_pattern: str) -> np.ndarray:
        columns = []
        i = 0
        motif_pattern = motif_pattern.upper()

        while i < len(motif_pattern):
            char = motif_pattern[i]
            col = np.zeros(4)

            if char == 'N':
                col[:] = 0.25
                columns.append(col)

            elif char == '[':
                end_id = motif_pattern.find(']', i)
                if end_id == -1:
                    raise ValueError(f"Padrão malformado: falta ']' após posição {i}")

                group = motif_pattern[i+1:end_id]
                prob = 1.0 / len(group)
                for base in group:
                    if base in cls.BASE_TO_ID:
                        col[cls.BASE_TO_ID[base]] = prob
                columns.append(col)
                i = end_id

            elif char in cls.BASE_TO_ID:
                col[cls.BASE_TO_ID[char]] = 1.0
                columns.append(col)

            i += 1

        return np.array(columns).T

    @classmethod
    def build_integer_pwm(cls, pfm: np.ndarray, bg_freqs: Dict[str, float]) -> np.ndarray:
        B = np.array([bg_freqs['A'], bg_freqs['C'], bg_freqs['G'], bg_freqs['T']])

        pseudo_fraction = 0.0001
        lambda_scale = 100.0

        # Cálculo vetorizado:
        # P_ib + c*Q_b
        numerator = pfm + (pseudo_fraction * B[:, np.newaxis])
        # Q_b + c*Q_b
        denominator = B[:, np.newaxis] + (pseudo_fraction * B[:, np.newaxis])

        # Log Odds
        log_odds = np.log2(numerator / denominator)

        # Escalonamento para Inteiros
        integer_pwm = np.round(log_odds * lambda_scale).astype(np.int32)

        return integer_pwm

    @classmethod
    def get_reverse_complement(cls, pwm: np.ndarray) -> np.ndarray:
        reversed_cols = np.flip(pwm, axis=1)
        rc_pwm = reversed_cols[[3, 2, 1, 0], :]
        return rc_pwm