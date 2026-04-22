import json
import math

from fastapi import APIRouter, UploadFile, File, Form, HTTPException

from core.fasta import Fasta
from core.motif_regex import MotifRegex
from core.p_value import PValue

router = APIRouter()

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

            p_value = PValue(motif_sequence, record.sequence)

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
                    "log_p_value": -math.log10(p_value.value)
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