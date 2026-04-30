from typing import Tuple

import numpy as np


class SequenceScanner:

    _ASCII_MAP = np.full(256, -1, dtype=np.int8)

    _ASCII_MAP[ord("A")] = 0
    _ASCII_MAP[ord("C")] = 1
    _ASCII_MAP[ord("G")] = 2
    _ASCII_MAP[ord("T")] = 3

    @classmethod
    def scan_integer_matrix(cls, sequence: str, pwm: np.ndarray, threshold: int) -> Tuple[np.ndarray, np.ndarray]:
        L = len(sequence)
        W = pwm.shape[1]

        if L < W:
            return np.array([]), np.array([])

        seq_bytes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        seq_encoded = cls._ASCII_MAP[seq_bytes]

        window_count = L - W + 1
        scores = np.zeros(window_count, dtype=np.int32)
        valid_bases_per_window = np.zeros(window_count, dtype=np.int16)

        for i in range(W):
            bases_at_i = seq_encoded[i:window_count + i]
            valid = bases_at_i >= 0
            safe_bases = np.where(valid, bases_at_i, 0)

            scores += np.where(valid, pwm[safe_bases, i], 0)
            valid_bases_per_window += valid

        valid_windows = valid_bases_per_window == W
        hits = np.where(valid_windows & (scores >= threshold))[0]

        return hits, scores[hits]

    @classmethod
    def scan_log_odds_matrix(
        cls,
        sequence: str,
        pwm: np.ndarray,
        threshold: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        L = len(sequence)
        W = pwm.shape[1]

        if L < W:
            return np.array([]), np.array([])

        seq_bytes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        seq_encoded = cls._ASCII_MAP[seq_bytes]

        window_count = L - W + 1
        scores = np.zeros(window_count, dtype=np.float64)
        valid_bases_per_window = np.zeros(window_count, dtype=np.int16)

        for i in range(W):
            bases_at_i = seq_encoded[i:window_count + i]
            valid = bases_at_i >= 0
            safe_bases = np.where(valid, bases_at_i, 0)

            scores += np.where(valid, pwm[safe_bases, i], 0.0)
            valid_bases_per_window += valid

        valid_windows = valid_bases_per_window == W
        hits = np.where(valid_windows & (scores >= threshold))[0]

        return hits, scores[hits]
