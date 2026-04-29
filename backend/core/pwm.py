from typing import Dict

import numpy as np


class PositionWeightMatrix:

    BASE_TO_ID = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
    }

    IUPAC_CODES = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "AG",
        "Y": "CT",
        "S": "CG",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT",
    }

    @classmethod
    def parse_motif_pattern_to_pfm(cls, motif_pattern: str) -> np.ndarray:
        columns = []
        i = 0
        motif_pattern = motif_pattern.strip().upper()

        if not motif_pattern:
            raise ValueError("O motivo nao pode estar vazio.")

        while i < len(motif_pattern):
            char = motif_pattern[i]
            col = np.zeros(4)

            if char == "[":
                end_id = motif_pattern.find("]", i)
                if end_id == -1:
                    raise ValueError(f"Padrao malformado: falta ']' apos posicao {i}")

                group = motif_pattern[i + 1:end_id]
                if not group:
                    raise ValueError(f"Padrao malformado: grupo vazio na posicao {i}")

                invalid = [base for base in group if base not in cls.BASE_TO_ID]
                if invalid:
                    raise ValueError(
                        f"Bases invalidas no grupo [{group}]: {', '.join(sorted(set(invalid)))}"
                    )

                prob = 1.0 / len(group)
                for base in group:
                    col[cls.BASE_TO_ID[base]] = prob
                columns.append(col)
                i = end_id

            elif char in cls.IUPAC_CODES:
                bases = cls.IUPAC_CODES[char]
                prob = 1.0 / len(bases)
                for base in bases:
                    col[cls.BASE_TO_ID[base]] = prob
                columns.append(col)

            else:
                raise ValueError(f"Simbolo IUPAC invalido '{char}' na posicao {i + 1}.")

            i += 1

        if not columns:
            raise ValueError("Nenhuma coluna valida foi gerada a partir do motivo.")

        return np.array(columns).T

    @classmethod
    def build_integer_pwm(
            cls,
            pfm: np.ndarray,
            bg_freqs: Dict[str, float],
            pseudo_fraction: float = 0.0001,
            lambda_scale: float = 100.0,
    ) -> np.ndarray:
        if pseudo_fraction < 0:
            raise ValueError("pseudo_fraction deve ser >= 0.")
        if lambda_scale <= 0:
            raise ValueError("lambda_scale deve ser > 0.")

        B = np.array([bg_freqs["A"], bg_freqs["C"], bg_freqs["G"], bg_freqs["T"]])

        numerator = pfm + (pseudo_fraction * B[:, np.newaxis])
        denominator = B[:, np.newaxis] + (pseudo_fraction * B[:, np.newaxis])

        log_odds = np.log2(numerator / denominator)
        integer_pwm = np.round(log_odds * lambda_scale).astype(np.int32)

        return integer_pwm

    @classmethod
    def build_log_odds_pwm(
            cls,
            pfm: np.ndarray,
            bg_freqs: Dict[str, float],
            pseudo_fraction: float = 0.0001,
    ) -> np.ndarray:
        if pseudo_fraction < 0:
            raise ValueError("pseudo_fraction deve ser >= 0.")

        B = np.array([bg_freqs["A"], bg_freqs["C"], bg_freqs["G"], bg_freqs["T"]], dtype=np.float64)
        numerator = pfm + (pseudo_fraction * B[:, np.newaxis])
        denominator = B[:, np.newaxis] + (pseudo_fraction * B[:, np.newaxis])
        return np.log2(numerator / denominator)

    @classmethod
    def get_reverse_complement(cls, pwm: np.ndarray) -> np.ndarray:
        reversed_cols = np.flip(pwm, axis=1)
        rc_pwm = reversed_cols[[3, 2, 1, 0], :]
        return rc_pwm

    @classmethod
    def reverse_complement_sequence(cls, sequence: str) -> str:
        complement = str.maketrans(
            "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
            "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
        )
        return sequence.translate(complement)[::-1].upper()
