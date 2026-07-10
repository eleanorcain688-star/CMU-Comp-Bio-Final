"""
bio_utils.py
Shared helpers for the alignment scripts: FASTA parsing and the
BLOSUM62 substitution matrix (the standard scoring matrix used by
real tools like BLAST for protein alignment).
"""

import re

# ---------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------

def parse_fasta(path):
    """
    Parse a FASTA file using this project's header convention:
        >CLASS|SPECIES|GENE|ACCESSION description
    Returns a list of dicts: {id, cls, species, gene, accession, desc, seq}
    Skips PLACEHOLDER entries by default is NOT done here — filtering
    happens in run_analysis.py so you can choose whether to include them.
    """
    records = []
    header = None
    seq_chunks = []

    def flush():
        if header is not None:
            seq = "".join(seq_chunks).upper()
            seq = re.sub(r"[^A-Z]", "", seq)  # strip whitespace/stray chars
            parts = header.split("|")
            cls = parts[0] if len(parts) > 0 else "UNKNOWN"
            species = parts[1] if len(parts) > 1 else "unknown"
            gene = parts[2] if len(parts) > 2 else "unknown"
            rest = parts[3] if len(parts) > 3 else ""
            rest_fields = rest.split(None, 1)
            accession = rest_fields[0] if rest_fields else "unknown"
            desc = rest_fields[1] if len(rest_fields) > 1 else ""
            rec_id = f"{cls}|{species}|{gene}|{accession}"
            records.append({
                "id": rec_id,
                "cls": cls,
                "species": species,
                "gene": gene,
                "accession": accession,
                "desc": desc,
                "seq": seq,
                "is_placeholder": "PLACEHOLDER" in header,
            })

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    flush()
    return records


# ---------------------------------------------------------------
# BLOSUM62 substitution matrix
# (standard 20-amino-acid scoring matrix used by BLAST and most
#  alignment tools; higher score = more conservative substitution)
# ---------------------------------------------------------------

_BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"

_BLOSUM62_ROWS = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]

BLOSUM62 = {}
for i, a in enumerate(_BLOSUM62_ORDER):
    for j, b in enumerate(_BLOSUM62_ORDER):
        BLOSUM62[(a, b)] = _BLOSUM62_ROWS[i][j]


def blosum_score(a, b):
    """Look up the BLOSUM62 score for a pair of residues.
    Unknown/ambiguous residues (X, *, etc.) score 0."""
    return BLOSUM62.get((a, b), 0)
