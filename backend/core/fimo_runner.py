import csv
import os
import shlex
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable

import numpy as np

from core.fasta_reader import FastaRecord
from core.pwm import PositionWeightMatrix


@dataclass
class FimoRecord:
    sequence_name: str
    start: int
    end: int
    strand_symbol: str
    score: float
    p_value: float
    matched_sequence: str


class FimoRunner:
    @staticmethod
    def _resolve_local_binary() -> str | None:
        binary = shutil.which("fimo") or shutil.which("fimo.exe")
        return binary

    @staticmethod
    def _resolve_docker_container() -> str | None:
        container_from_env = os.getenv("FIMO_DOCKER_CONTAINER", "").strip()
        if container_from_env:
            return container_from_env

        try:
            result = subprocess.run(
                ["docker", "ps", "--format", "{{.Names}}"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=False,
            )
        except Exception:
            return None

        if result.returncode != 0:
            return None

        names = [line.strip() for line in (result.stdout or "").splitlines() if line.strip()]
        if not names:
            return None

        for candidate in names:
            lowered = candidate.lower()
            if "meme" in lowered:
                return candidate

        return None

    @staticmethod
    def _check_fimo_in_container(container_name: str) -> bool:
        result = subprocess.run(
            ["docker", "exec", container_name, "sh", "-lc", "which fimo"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        return result.returncode == 0 and bool((result.stdout or "").strip())

    @staticmethod
    def _build_meme_motif_file(
        motif_pattern: str,
        pfm: np.ndarray,
        bg_freqs: Dict[str, float],
        output_path: Path,
    ) -> None:
        width = pfm.shape[1]
        lines = [
            "MEME version 5",
            "",
            "ALPHABET= ACGT",
            "strands: + -",
            "",
            "Background letter frequencies",
            f"A {bg_freqs['A']:.6f} C {bg_freqs['C']:.6f} G {bg_freqs['G']:.6f} T {bg_freqs['T']:.6f}",
            "",
            f"MOTIF USER_MOTIF {motif_pattern}",
            f"letter-probability matrix: alength= 4 w= {width} nsites= 20 E= 0",
        ]

        for i in range(width):
            row = [pfm[0, i], pfm[1, i], pfm[2, i], pfm[3, i]]
            lines.append(" ".join(f"{value:.6f}" for value in row))

        output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    @staticmethod
    def _build_fasta_file(records: Iterable[FastaRecord], output_path: Path) -> dict[str, FastaRecord]:
        id_to_record: dict[str, FastaRecord] = {}
        lines: list[str] = []

        for idx, record in enumerate(records, start=1):
            seq_id = f"SEQ_{idx:07d}"
            id_to_record[seq_id] = record
            lines.append(f">{seq_id}")
            lines.append(record.sequence)

        output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return id_to_record

    @staticmethod
    def _parse_tsv(stdout_text: str) -> list[FimoRecord]:
        rows: list[FimoRecord] = []
        filtered_lines = [line for line in stdout_text.splitlines() if line and not line.startswith("#")]
        if not filtered_lines:
            return rows

        reader = csv.DictReader(filtered_lines, delimiter="\t")
        for row in reader:
            try:
                sequence_name = (row.get("sequence_name") or "").strip()
                start = int(row.get("start") or 0)
                end = int(row.get("stop") or 0)
                strand_symbol = (row.get("strand") or "+").strip() or "+"
                score = float(row.get("score") or 0.0)
                p_value = float(row.get("p-value") or row.get("pvalue") or 1.0)
                matched_sequence = (row.get("matched_sequence") or "").strip().upper()
            except Exception:
                continue

            if not sequence_name or start <= 0 or end <= 0:
                continue

            rows.append(
                FimoRecord(
                    sequence_name=sequence_name,
                    start=start,
                    end=end,
                    strand_symbol=strand_symbol,
                    score=score,
                    p_value=p_value,
                    matched_sequence=matched_sequence,
                )
            )

        return rows

    @classmethod
    def scan(
        cls,
        motif_pattern: str,
        records: list[FastaRecord],
        bg_freqs: Dict[str, float],
        p_value_threshold: float,
        scan_reverse_complement: bool,
    ) -> list[dict]:
        pfm = PositionWeightMatrix.parse_motif_pattern_to_pfm(motif_pattern)
        local_fimo_bin = cls._resolve_local_binary()
        docker_container = None if local_fimo_bin else cls._resolve_docker_container()

        if not local_fimo_bin and docker_container and not cls._check_fimo_in_container(docker_container):
            docker_container = None

        if not local_fimo_bin and not docker_container:
            raise RuntimeError(
                "FIMO nao encontrado no PATH local nem em container Docker. "
                "Defina FIMO_DOCKER_CONTAINER=memesuite ou disponibilize 'fimo' no sistema."
            )

        with tempfile.TemporaryDirectory(prefix="motif_fimo_") as tmp_dir:
            tmp_path = Path(tmp_dir)
            motif_file = tmp_path / "motif.meme"
            fasta_file = tmp_path / "sequences.fasta"

            cls._build_meme_motif_file(motif_pattern, pfm, bg_freqs, motif_file)
            id_to_record = cls._build_fasta_file(records, fasta_file)

            if local_fimo_bin:
                command = [
                    local_fimo_bin,
                    "--text",
                    "--thresh",
                    str(p_value_threshold),
                ]
                if not scan_reverse_complement:
                    command.append("--norc")
                command.extend([str(motif_file), str(fasta_file)])
                completed = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )
            else:
                container_dir = f"/tmp/motif_fimo_{os.getpid()}_{uuid.uuid4().hex[:10]}"
                motif_container_path = f"{container_dir}/motif.meme"
                fasta_container_path = f"{container_dir}/sequences.fasta"

                create_dir = subprocess.run(
                    ["docker", "exec", docker_container, "mkdir", "-p", container_dir],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )
                if create_dir.returncode != 0:
                    raise RuntimeError(f"Falha ao preparar diretório no container: {(create_dir.stderr or '').strip()}")

                cp_motif = subprocess.run(
                    ["docker", "cp", str(motif_file), f"{docker_container}:{motif_container_path}"],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )
                if cp_motif.returncode != 0:
                    raise RuntimeError(f"Falha ao copiar motif para container: {(cp_motif.stderr or '').strip()}")

                cp_fasta = subprocess.run(
                    ["docker", "cp", str(fasta_file), f"{docker_container}:{fasta_container_path}"],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )
                if cp_fasta.returncode != 0:
                    raise RuntimeError(f"Falha ao copiar FASTA para container: {(cp_fasta.stderr or '').strip()}")

                run_parts = ["fimo", "--text", "--thresh", str(p_value_threshold)]
                if not scan_reverse_complement:
                    run_parts.append("--norc")
                run_parts.extend([motif_container_path, fasta_container_path])
                run_command = " ".join(shlex.quote(part) for part in run_parts)

                completed = subprocess.run(
                    ["docker", "exec", docker_container, "sh", "-lc", run_command],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )

                subprocess.run(
                    ["docker", "exec", docker_container, "rm", "-rf", container_dir],
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    check=False,
                )

            if completed.returncode != 0:
                stderr = (completed.stderr or completed.stdout or "").strip()
                raise RuntimeError(f"Falha ao executar FIMO: {stderr}")

            raw_rows = cls._parse_tsv(completed.stdout)
            matches: list[dict] = []

            for row in raw_rows:
                record = id_to_record.get(row.sequence_name)
                if not record:
                    continue

                record_length = len(record.sequence)
                strand_symbol = row.strand_symbol if row.strand_symbol in {"+", "-"} else "+"
                strand_name = "forward" if strand_symbol == "+" else "reverse"
                source_strand_symbol = strand_symbol

                if record.is_reverse:
                    source_start = record_length - row.end + 1
                    source_end = record_length - row.start + 1
                    source_strand_symbol = "-" if strand_symbol == "+" else "+"
                else:
                    source_start = row.start
                    source_end = row.end

                source_strand = "forward" if source_strand_symbol == "+" else "reverse"
                neg_log10_pval = round(-np.log10(row.p_value), 2) if row.p_value > 0 else None

                matches.append(
                    {
                        "gene": record.gene_id,
                        "record_orientation": "reverse_complement" if record.is_reverse else "forward",
                        "strand": strand_name,
                        "strand_symbol": strand_symbol,
                        "source_strand": source_strand,
                        "source_strand_symbol": source_strand_symbol,
                        "matched_sequence": row.matched_sequence
                        or record.sequence[row.start - 1:row.end],
                        "position": row.start,
                        "end_position": row.end,
                        "source_position": source_start,
                        "source_end_position": source_end,
                        "score": row.score,
                        "p_value": f"{row.p_value:.2e}",
                        "p_value_numeric": row.p_value,
                        "-log10_pval": neg_log10_pval,
                    }
                )

            matches.sort(key=lambda item: (item["p_value_numeric"], -item["score"]))
            return matches
