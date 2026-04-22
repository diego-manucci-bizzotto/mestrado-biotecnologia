from typing import Dict, List
from collections import defaultdict


class PWMPValue:
    def __init__(self, pwm: List[Dict[str, int]], background_frequencies: Dict[str, float]):
        """
        pwm: Lista de dicionários representando a Matriz de Pesos (PWM).
             Cada item da lista é uma posição do motivo.
             Ex: [{'A': 2, 'C': -1, 'G': -1, 'T': 0}, {'A': -2, 'C': 3, ...}]
             * Nota: Para seguir o artigo e evitar lentidão extrema, a matriz deve
               usar números inteiros (granularidade) e não decimais complexos.
        bg_freqs: Frequências de fundo. Ex: {'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}
        """
        self.pwm = pwm
        self.background_frequencies = background_frequencies
        self.m = len(pwm)  # Tamanho do motivo

    def calculate_p_value(self, threshold_score: int) -> float:

        q = {0: 1.0}  # Q(M[1..0], 0) = 1

        for motif_position_scores in self.pwm:
            next_q = defaultdict(float)
            for current_score, probability in q.items():
                for base, background_probability in self.background_frequencies.items():
                    score_increment = motif_position_scores.get(base, 0)
                    new_score = current_score + score_increment
                    next_q[new_score] += probability * background_probability
            q = next_q
        p_value = sum(prob for score, prob in q.items() if score >= threshold_score)

        return p_value
