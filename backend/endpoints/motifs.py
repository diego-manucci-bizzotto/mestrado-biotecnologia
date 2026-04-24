import json
import math
from typing import Optional

from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from starlette.responses import StreamingResponse

from core.bg_freqs import BackgroundFrequencies
from core.fasta_reader import FastaReader
from core.p_value_calc import PValueCalculator
from core.pwm import PositionWeightMatrix
from core.seq_scanner import SequenceScanner

router = APIRouter()


@router.post("/", status_code=200)
async def find_motifs(
        file: UploadFile = File(...),
        motif_pattern: str = Form(..., description="IUPAC motif pattern with optional bracket groups"),
        p_value_threshold: float = Form(1e-4),
        score_threshold: Optional[int] = Form(None),
        scan_reverse_complement: bool = Form(True),
):
    filename = file.filename or ""
    if not filename.endswith((".fasta", ".fa", ".fna")):
        raise HTTPException(status_code=400, detail="O formato deve ser FASTA.")

    if score_threshold is None and not 0 < p_value_threshold <= 1:
        raise HTTPException(
            status_code=400,
            detail="p_value_threshold deve estar no intervalo (0, 1].",
        )

    bg_freqs = await BackgroundFrequencies.calculate_from_file(file)
    await file.seek(0)

    try:
        pfm = PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern)
        pwm_fwd = PositionWeightMatrix.build_integer_pwm(pfm, bg_freqs)
        pcalc_fwd = PValueCalculator(pwm_fwd, bg_freqs)

        pwm_rev = PositionWeightMatrix.get_reverse_complement(pwm_fwd)
        pcalc_rev = PValueCalculator(pwm_rev, bg_freqs)

        motif_length = pwm_fwd.shape[1]
        if score_threshold is None:
            cutoff_fwd = pcalc_fwd.get_score_threshold_for_pvalue(p_value_threshold)
            cutoff_rev = pcalc_rev.get_score_threshold_for_pvalue(p_value_threshold)
            cutoff_mode = "pvalue"
        else:
            cutoff_fwd = score_threshold
            cutoff_rev = score_threshold
            cutoff_mode = "score"

    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Erro na construcao da matriz: {str(e)}")

    async def process_and_stream():
        fasta_reader = FastaReader(file)
        scan_configs = [
            ("forward", "+", pwm_fwd, pcalc_fwd, cutoff_fwd),
        ]

        if scan_reverse_complement:
            scan_configs.append(("reverse", "-", pwm_rev, pcalc_rev, cutoff_rev))

        matches = []

        async for record in fasta_reader.parse():
            record_length = len(record.sequence)
            record_orientation = "reverse_complement" if record.is_reverse else "forward"

            for strand_name, strand_symbol, pwm, pcalc, cutoff in scan_configs:
                hit_positions, hit_scores = SequenceScanner.scan_integer_matrix(
                    sequence=record.sequence,
                    pwm=pwm,
                    threshold=cutoff,
                )

                for pos_idx, score in zip(hit_positions, hit_scores):
                    score = int(score)
                    p_value = float(pcalc.get_pvalue(score))

                    if score_threshold is None and p_value > p_value_threshold:
                        continue

                    start = int(pos_idx) + 1
                    end = start + motif_length - 1
                    sequence_window = record.sequence[pos_idx:pos_idx + motif_length]
                    matched_sequence = (
                        PositionWeightMatrix.reverse_complement_sequence(sequence_window)
                        if strand_symbol == "-"
                        else sequence_window
                    )

                    if record.is_reverse:
                        source_start = record_length - end + 1
                        source_end = record_length - start + 1
                        source_strand_symbol = "-" if strand_symbol == "+" else "+"
                    else:
                        source_start = start
                        source_end = end
                        source_strand_symbol = strand_symbol

                    source_strand = "reverse" if source_strand_symbol == "-" else "forward"

                    neg_log10_pval = round(-math.log10(p_value), 2) if p_value > 0 else None

                    matches.append({
                        "gene": record.gene_id,
                        "record_orientation": record_orientation,
                        "strand": strand_name,
                        "strand_symbol": strand_symbol,
                        "source_strand": source_strand,
                        "source_strand_symbol": source_strand_symbol,
                        "matched_sequence": matched_sequence,
                        "position": start,
                        "end_position": end,
                        "source_position": source_start,
                        "source_end_position": source_end,
                        "score": score,
                        "p_value": f"{p_value:.2e}",
                        "p_value_numeric": p_value,
                        "-log10_pval": neg_log10_pval,
                    })

        yield json.dumps({
            "metadata": {
                "filename": filename,
                "motif_pattern": motif_pattern,
                "motif_length": motif_length,
                "score_threshold": score_threshold,
                "cutoff_mode": cutoff_mode,
                "background_frequencies": bg_freqs,
            },
            "matches": matches,
        })

    return StreamingResponse(process_and_stream(), media_type="application/x-ndjson")
