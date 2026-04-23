import json
import math
from core.fasta import Fasta
from core.motif_regex import MotifRegex
from core.p_value import PValue
from core.pwm_p_value import PWMPValue
from core2.bg_freqs import BackgroundFrequencies
from core2.fasta_reader import FastaReader
from core2.p_value_calc import PValueCalculator
from core2.pwm import PositionWeightMatrix
from core2.seq_scanner import SequenceScanner
from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from starlette.responses import StreamingResponse

router = APIRouter()


@router.post("/", status_code=200)
async def find_motifs(
        file: UploadFile = File(...),
        motif_pattern: str = Form(..., description="IUPAC motif pattern or regex with brackets"),
):
    if not file.filename.endswith(('.fasta', '.fa', '.fna')):
        raise HTTPException(status_code=400, detail="O formato deve ser FASTA.")

    bg_freqs = await BackgroundFrequencies.calculate_from_file(file)
    await file.seek(0)  # Volta o ponteiro do arquivo para o início para a leitura subsequente

    # 3. COMPILAÇÃO DA MATRIZ E CALCULADORA DE P-VALOR
    try:
        # Fita Forward
        pfm = PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern)
        pwm_fwd = PositionWeightMatrix.build_integer_pwm(pfm, bg_freqs)
        pcalc_fwd = PValueCalculator(pwm_fwd, bg_freqs)

        # Fita Reverse Complement
        pwm_rev = PositionWeightMatrix.get_reverse_complement(pwm_fwd)
        pcalc_rev = PValueCalculator(pwm_rev, bg_freqs)

        W = pwm_fwd.shape[1]
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Erro na construção da matriz: {str(e)}")

    # 4. THRESHOLD DE CORTE (Rigor Estatístico)
    # Exigimos no mínimo um P-valor de 1e-4 para não poluir o Front-end com falsos positivos
    P_VALUE_CUTOFF = 1e-4

    # 5. MOTOR DE VARREDURA (Streaming)
    async def process_and_stream():
        fasta_reader = FastaReader(file)

        async for record in fasta_reader.parse():
            current_pwm = pwm_rev if record.is_reverse else pwm_fwd
            current_pcalc = pcalc_rev if record.is_reverse else pcalc_fwd

            # Varredura (Só pedimos hits com score > 0 para acelerar, filtramos o rigor a seguir)
            hit_positions, hit_scores = SequenceScanner.scan_integer_matrix(
                sequence=record.sequence,
                pwm=current_pwm,
                threshold=0
            )

            for pos_idx, score in zip(hit_positions, hit_scores):
                # Conversão de Score C-Level para P-Valor Biológico
                p_value = current_pcalc.get_pvalue(int(score))

                # Só emitimos se for estatisticamente significante
                if p_value <= P_VALUE_CUTOFF:
                    motif_seq = record.sequence[pos_idx: pos_idx + W]

                    final_pos = pos_idx + 1
                    if record.is_reverse:
                        final_pos = 1000 - final_pos - W + 2

                    yield json.dumps({
                        "gene": record.gene_id,
                        "strand": "reverse" if record.is_reverse else "forward",
                        "motif": motif_seq,
                        "position": int(final_pos),
                        "p_value": f"{p_value:.2e}",
                        "-log10_pval": round(-math.log10(p_value) if p_value > 0 else 0, 2)
                    }) + "\n"

    return StreamingResponse(process_and_stream(), media_type="application/x-ndjson")
