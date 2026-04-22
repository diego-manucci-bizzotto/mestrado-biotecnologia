import math

from fastapi import APIRouter, UploadFile, File, Form, HTTPException

from core.fasta import Fasta
from core.motif_regex import MotifRegex
from core.p_value import PValue
from core.pwm_p_value import PWMPValue

router = APIRouter()

raw_matrix = [
    {'A': 20, 'C': 48, 'G': 18, 'T': 14},  # Pos 1
    {'A': 28, 'C': 22, 'G': 17, 'T': 32},  # Pos 2
    {'A': 37, 'C': 30, 'G': 22, 'T': 11},  # Pos 3
    {'A': 63, 'C': 36, 'G': 1, 'T': 1},  # Pos 4
    {'A': 0, 'C': 0, 'G': 100, 'T': 0},  # Pos 5 (G)
    {'A': 0, 'C': 100, 'G': 0, 'T': 0},  # Pos 6 (C)
    {'A': 0, 'C': 99, 'G': 0, 'T': 0},  # Pos 7 (C)
    {'A': 42, 'C': 30, 'G': 7, 'T': 21},  # Pos 8
    {'A': 7, 'C': 76, 'G': 8, 'T': 8}  # Pos 9
]

bg_freqs = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
pseudocount = 0.01
epsilon = 0.1
pwm_matrix = []

# Transformando porcentagens em Log-Razão Inteira (Granularidade da Definição 3 do artigo)
for col in raw_matrix:
    col_inteira = {}
    for base in ['A', 'C', 'G', 'T']:
        prob = (col[base] / 100.0) + pseudocount
        prob_normalizada = prob / (1.0 + (4 * pseudocount))
        log_odds = math.log2(prob_normalizada / bg_freqs[base])
        score_discreto = math.floor(log_odds / epsilon)
        col_inteira[base] = score_discreto
    pwm_matrix.append(col_inteira)

pwm_calculator = PWMPValue(pwm=pwm_matrix, background_frequencies=bg_freqs)


@router.post("/", status_code=200)
async def find_motifs(
        file: UploadFile = File(...),
        regex: str = Form(...)
):
    if not file.filename.endswith(('.fasta', '.fa')):
        raise HTTPException(status_code=400, detail="Ficheiro inválido.")

    if not regex:
        raise HTTPException(status_code=400, detail="Padrão regex é obrigatório.")

    motif_regex = MotifRegex(regex)

    results = []
    fasta_stream = Fasta(file)

    async for record in fasta_stream:
        regex_pattern = motif_regex.forward_regex_pattern if not record.is_reverse else motif_regex.reverse_regex_pattern

        for match in regex_pattern.finditer(record.sequence):
            motif_sequence = match.group(1).upper()

            if len(motif_sequence) != pwm_calculator.m:
                continue  # Pula se não tiver o tamanho exato da PWM

            # 4. CALCULAR O SCORE ALVO DA SEQUENCIA ENCONTRADA (Ex: "CAGCGCCCA")
            hit_score = 0
            for i, base in enumerate(motif_sequence):
                hit_score += pwm_matrix[i].get(base, 0)

            # Cálculo de background isolado (o seu original)
            p_value = PValue(motif_sequence, record.sequence)

            # 5. GERAR O P-VALUE EXATO DA PWM PASSANDO O SCORE ALVO
            pwm_p_val = pwm_calculator.calculate_p_value(threshold_score=hit_score)

            pos = match.start() + 1
            if record.is_reverse:
                pos = 1000 - pos - len(motif_sequence) + 2

            results.append(
                {
                    "gene": record.gene_id,
                    "type": "reverse" if record.is_reverse else "forward",
                    "motif": motif_sequence,
                    "position": pos,
                    "p_value": f"{p_value.value:.2e}",
                    "log_p_value": -math.log10(p_value.value) if p_value.value > 0 else 0,
                    "pwm_p_value": f"{pwm_p_val:.2e}"
                }
            )

    return {
        "metadata": {
            "filename": file.filename,
            "patterns": {"forward": motif_regex.forward_regex_pattern, "reverse": motif_regex.reverse_regex_pattern},
            "summary": {
                "total_hits": len(results),
            }
        },
        "hits": results
    }

#   if results:
#       return "".join(json.dumps(res) + "\n" for res in results)
#   return ""
#
#   return {
#       "metadata": {"filename": file.filename, "pattern": forward_pattern},
#       "hits": results
#   }
#
# StreamingResponse(fasta_stream_generator(), media_type="application/x-ndjson")
