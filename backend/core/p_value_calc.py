import numpy as np
from typing import Dict

class PValueCalculator:
    def __init__(self, integer_pwm: np.ndarray, bg_freqs: Dict[str, float]):
        self.pwm = integer_pwm
        self.W = integer_pwm.shape[1]
        self.bg = np.array([bg_freqs['A'], bg_freqs['C'], bg_freqs['G'], bg_freqs['T']])

        # 1. Mapear os limites termodinâmicos possíveis da Matriz
        self.min_score = int(np.sum(np.min(self.pwm, axis=0)))
        self.max_score = int(np.sum(np.max(self.pwm, axis=0)))

        # Array para armazenar as probabilidades cumulativas (P-value)
        # O offset lida com índices negativos no array
        self.range_size = self.max_score - self.min_score + 1
        self.offset = -self.min_score
        self.pvalue_table = np.zeros(self.range_size, dtype=np.float64)

        # Calcula a distribuição no momento da inicialização
        self._build_distribution()

    def _build_distribution(self):
        # current_pdf guarda a probabilidade de atingir um 'score' na posição atual
        current_pdf = {0: 1.0}

        for i in range(self.W):
            next_pdf = {}
            for b in range(4):
                weight = int(self.pwm[b, i])
                prob_b = self.bg[b]

                # Para cada score já alcançado, somamos o peso da nova base
                for prev_score, prev_prob in current_pdf.items():
                    new_score = prev_score + weight
                    # A probabilidade de chegar aqui é P(Anterior) * P(Base Atual)
                    next_pdf[new_score] = next_pdf.get(new_score, 0.0) + (prev_prob * prob_b)

            current_pdf = next_pdf

        # Função de Sobrevivência (Survivor Function): P(X >= S)
        # Varremos de trás para frente somando as probabilidades para gerar o P-valor
        accumulated_p = 0.0
        for s in range(self.max_score, self.min_score - 1, -1):
            if s in current_pdf:
                accumulated_p += current_pdf[s]
            # Salva o P-valor exato para aquele score no array de consulta O(1)
            self.pvalue_table[s + self.offset] = accumulated_p

    def get_pvalue(self, score: int) -> float:
        if score >= self.max_score:
            return self.pvalue_table[self.max_score + self.offset]
        if score <= self.min_score:
            return 1.0
        return self.pvalue_table[score + self.offset]

    def get_score_threshold_for_pvalue(self, pvalue_threshold: float) -> int:
        if pvalue_threshold <= 0:
            return self.max_score + 1
        if pvalue_threshold >= 1:
            return self.min_score

        for score in range(self.min_score, self.max_score + 1):
            if self.pvalue_table[score + self.offset] <= pvalue_threshold:
                return score

        return self.max_score + 1
