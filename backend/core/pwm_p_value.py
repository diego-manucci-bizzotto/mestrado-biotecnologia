from collections import defaultdict
from typing import Dict, List


class PWMPValue:
    def __init__(self, pwm: List[Dict[str, int]], background_frequencies: Dict[str, float]):
        self.pwm = pwm
        self.background_frequencies = background_frequencies
        self.m = len(pwm)  # Tamanho do motivo (Neste caso, 9)

        # PRE-CALCULA A DISTRIBUIÇÃO DE PONTUAÇÕES (Equação 1 do artigo)
        self.q_distribution = self._build_score_distribution()

    def _build_score_distribution(self) -> Dict[int, float]:
        q = {0: 1.0}
        for motif_position_scores in self.pwm:
            next_q = defaultdict(float)
            for current_score, probability in q.items():
                for base, background_probability in self.background_frequencies.items():
                    score_increment = motif_position_scores.get(base, 0)
                    new_score = current_score + score_increment
                    next_q[new_score] += probability * background_probability
            q = next_q
        return q

    def calculate_p_value(self, threshold_score: int) -> float:
        # Apenas soma probabilidades calculadas previamente (Equação 2 do artigo)
        return sum(prob for score, prob in self.q_distribution.items() if score >= threshold_score)

