import json
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from collections import Counter
from starlette.responses import StreamingResponse
from constants import COMPLEMENT
from re import Pattern, search as re_search, compile as re_compile, IGNORECASE, error as re_error, findall as re_findall

app = FastAPI()

def get_reverse_complement_regex(regex: str) -> str:
    tokens = re_findall(r'\[.*?]|.', regex)
    reverse_complement_tokens = []
    for token in reversed(tokens):
        if token.startswith('[') and token.endswith(']'):
            letters = token[1:-1]
            reverse_complement_letters = letters.translate(COMPLEMENT)[::-1]
            reverse_complement_tokens.append(f'[{reverse_complement_letters}]')
        else:
            reverse_complement_tokens.append(token.translate(COMPLEMENT))
    return ''.join(reverse_complement_tokens)

def get_sequence_base_proportions(sequence: str) -> dict:
    counts = Counter(sequence.upper())
    base_sum = sum(counts.values())

    if base_sum == 0:
        return {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    return {base: count / base_sum for base, count in counts.items() if base in "ACGTN"}

def process_record(header: str, sequence: str, is_reverse_strand: bool, forward_pattern: Pattern[str], reverse_pattern: Pattern[str]) -> str:
    results = []

    gene_match = re_search(r'\((.*?)\)', header)
    gene_name = gene_match.group(1) if gene_match else header.split()[0].replace('>', '')

    sequence_base_proportions = get_sequence_base_proportions(sequence)

    def map_hits(pattern, strand_type):
        for match in pattern.finditer(sequence):
            motif_sequence = match.group(1).upper()
            ##p_val = calculate_motif_p_value(motif_seq, bg_freqs)
            pos = match.start() + 1
            if is_reverse_strand:
                pos = 1000 - pos - len(motif_sequence) + 2

            results.append({
                "gene": gene_name,
                "type": strand_type,
                "motif": motif_sequence,
                "position": pos,
                #"p_value": f"{p_val:.2e}" # Formato científico exigido em Bioinformática
            })

    map_hits(forward_pattern, "forward")
    map_hits(reverse_pattern, "reverse")

    if results:
        return "".join(json.dumps(res) + "\n" for res in results)
    return ""



@app.post("/analyze-promoters")
async def analyze_promoters(
        file : UploadFile = File(...),
        regex: str = Form(...)
):
    try:
        forward_pattern = re_compile(f'(?=({regex.strip()}))', IGNORECASE)
        reverse_pattern = re_compile(f'(?=({get_reverse_complement_regex(regex.strip())}))', IGNORECASE)
    except re_error:
        raise HTTPException(status_code=400, detail="Invalid regex pattern.")

    async def fasta_stream_generator():
        current_header = ""
        current_seq = []
        is_reverse_strand = False

        # Lê blocos de forma assíncrona (segurança para ficheiros pesados)
        buffer = ""
        while True:
            chunk = await file.read(64 * 1024)  # Lê blocos de 64KB
            if not chunk:
                break

            buffer += chunk.decode('utf-8')
            lines = buffer.split('\n')

            # Deixa a última linha (que pode estar cortada a meio) para o próximo ciclo
            buffer = lines.pop()

            for line in lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Processa a sequência anterior antes de iniciar a nova
                    if current_header and current_seq:
                        output = process_record(
                            current_header,
                            "".join(current_seq),
                            is_reverse_strand,
                            forward_pattern,
                            reverse_pattern
                        )
                        if output:
                            yield output

                    # Atualiza o estado
                    current_header = line
                    is_reverse_strand = "[reverse complement]" in line.lower()
                    current_seq = []
                else:
                    current_seq.append(line)

        # Garante o processamento do último registo que ficou no buffer ou array
        if current_header:
            if buffer:  # Adiciona qualquer resto de sequência
                current_seq.append(buffer.strip())
            output = process_record(
                current_header,
                "".join(current_seq),
                is_reverse_strand,
                forward_pattern,
                reverse_pattern
            )
            if output:
                yield output

    # Devolvemos Content-Type especial: application/x-ndjson
    # Isto instrui os browsers e Workers a processarem linha-a-linha em tempo real.
    return StreamingResponse(fasta_stream_generator(), media_type="application/x-ndjson")