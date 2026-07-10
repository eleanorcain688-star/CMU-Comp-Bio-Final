"""
run_analysis.py
Runs Needleman-Wunsch, Smith-Waterman, and mini-BLAST on every pair of
sequences in data/sequences.fasta, and writes a single results.json
that the website reads to render heatmaps and tables.

Usage:
    python3 run_analysis.py [--include-placeholders]

By default, PLACEHOLDER records are skipped so you don't accidentally
publish results computed on filler sequences. Pass --include-placeholders
while testing the pipeline itself.
"""

import argparse
import json
import os
import sys
import time

from bio_utils import parse_fasta
from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman
from mini_blast import mini_blast

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "sequences.fasta")
OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "website", "results.json")


def compute_matrix(records, align_fn, score_key="percent_identity"):
    """Returns (matrix, alignments) where alignments[i][j] holds the
    aligned strings for i<j (upper triangle only, to keep JSON small;
    the matrix itself is still fully populated and symmetric)."""
    n = len(records)
    matrix = [[None] * n for _ in range(n)]
    alignments = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i][j] = 100.0 if score_key == "percent_identity" else 0
                continue
            if j < i:
                matrix[i][j] = matrix[j][i]  # symmetric, reuse
                continue
            result = align_fn(records[i]["seq"], records[j]["seq"])
            matrix[i][j] = round(result[score_key], 2)

            a1 = result.get("aligned_seq1") or result.get("seq1_region", "")
            a2 = result.get("aligned_seq2") or result.get("seq2_region", "")
            alignments[f"{i}-{j}"] = {"a1": a1, "a2": a2}
    return matrix, alignments


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-placeholders", action="store_true")
    args = parser.parse_args()

    records = parse_fasta(DATA_PATH)
    if not args.include_placeholders:
        records = [r for r in records if not r["is_placeholder"]]

    if len(records) < 2:
        print("Need at least 2 non-placeholder sequences to run. "
              "Re-run with --include-placeholders to test the pipeline, "
              "or add real sequences to data/sequences.fasta.")
        sys.exit(1)

    print(f"Loaded {len(records)} sequences:")
    for r in records:
        print(f"  {r['cls']:10s} {r['species']:20s} {r['gene']:10s} ({len(r['seq'])} aa)")

    labels = [f"{r['cls']}:{r['species']}" for r in records]
    classes = [r["cls"] for r in records]

    t0 = time.time()
    print("\nRunning Needleman-Wunsch (global)...")
    nw_matrix, nw_align = compute_matrix(records, needleman_wunsch)

    print("Running Smith-Waterman (local)...")
    sw_matrix, sw_align = compute_matrix(records, smith_waterman)

    print("Running mini-BLAST (seed-and-extend)...")
    bl_matrix, bl_align = compute_matrix(records, mini_blast)
    elapsed = time.time() - t0
    print(f"Done in {elapsed:.2f}s")

    output = {
        "labels": labels,
        "classes": classes,
        "species": [r["species"] for r in records],
        "genes": [r["gene"] for r in records],
        "accessions": [r["accession"] for r in records],
        "lengths": [len(r["seq"]) for r in records],
        "matrices": {
            "needleman_wunsch": nw_matrix,
            "smith_waterman": sw_matrix,
            "mini_blast": bl_matrix,
        },
        "alignments": {
            "needleman_wunsch": nw_align,
            "smith_waterman": sw_align,
            "mini_blast": bl_align,
        },
        "metric": "percent_identity",
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nWrote results to {OUT_PATH}")


if __name__ == "__main__":
    main()
