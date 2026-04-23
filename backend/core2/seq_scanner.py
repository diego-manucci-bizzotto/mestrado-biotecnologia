from typing import Tuple

import numpy as np


class SequenceScanner:

    _ASCII_MAP = np.zeros(256, dtype=np.int8)

    _ASCII_MAP[ord('A')] = 0
    _ASCII_MAP[ord('C')] = 1
    _ASCII_MAP[ord('G')] = 2
    _ASCII_MAP[ord('T')] = 3

    @classmethod
    def scan_integer_matrix(cls, sequence: str, pwm: np.ndarray, threshold: int) -> Tuple[np.ndarray, np.ndarray]:
        L = len(sequence)
        W = pwm.shape[1]

        # Validação biológica e estrutural
        if L < W:
            return np.array([]), np.array([])

        # -----------------------------------------------------------
        # OTIMIZAÇÃO 1: Mapeamento em Nível C (Bypass no Python)
        # -----------------------------------------------------------
        # Transforma a string de DNA num array de bytes (ASCII)
        seq_bytes = np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)
        # Mapeia instantaneamente os bytes para 0, 1, 2 ou 3
        seq_encoded = cls._ASCII_MAP[seq_bytes]

        # -----------------------------------------------------------
        # OTIMIZAÇÃO 2: Acumulador de Janela Deslizante
        # -----------------------------------------------------------
        # Em vez de um loop de tamanho L (1000), faremos um array de tamanho L-W+1
        # Usamos int32 para evitar estouro numérico na soma dos scores escalonados
        scores = np.zeros(L - W + 1, dtype=np.int32)

        # O truque de mestre: Nós iteramos apenas o tamanho do Motivo (W=8 vezes).
        # Somamos a coluna inteira do genoma de uma só vez usando fatiamento (slicing)
        for i in range(W):
            # Quais bases estão na posição 'i' da janela em todos os pontos do genoma?
            bases_at_i = seq_encoded[i : L - W + 1 + i]
            # Extrai o peso daquela base naquela coluna da matriz e soma ao acumulador global
            scores += pwm[bases_at_i, i]

        # -----------------------------------------------------------
        # FILTRAGEM TERMODINÂMICA
        # -----------------------------------------------------------
        # np.where roda em C e retorna apenas os índices que passaram no teste
        hits = np.where(scores >= threshold)[0]

        return hits, scores[hits]
