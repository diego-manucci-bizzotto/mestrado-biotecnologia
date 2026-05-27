from __future__ import annotations

IUPAC: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def normalize_iupac(pattern: str) -> str:
    return pattern.upper().replace(" ", "").replace("'", "").replace("-", "")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def matches_iupac(window: str, pattern: str) -> bool:
    if len(window) != len(pattern):
        return False
    return all(base in IUPAC.get(code, set()) for base, code in zip(window.upper(), pattern))


def iupac_to_pwm(pattern: str, pseudocount: float = 0.01) -> list[dict[str, float]]:
    matrix: list[dict[str, float]] = []
    for code in normalize_iupac(pattern):
        allowed = IUPAC.get(code)
        if not allowed:
            raise ValueError(f"Unsupported IUPAC code: {code}")
        raw = {base: pseudocount for base in "ACGT"}
        for base in allowed:
            raw[base] += 1.0 / len(allowed)
        total = sum(raw.values())
        matrix.append({base: raw[base] / total for base in "ACGT"})
    return matrix


def probability_for_iupac(pattern: str, background: dict[str, float]) -> float:
    probability = 1.0
    for code in normalize_iupac(pattern):
        allowed = IUPAC.get(code)
        if not allowed:
            raise ValueError(f"Unsupported IUPAC code: {code}")
        probability *= sum(background[base] for base in allowed)
    return probability
