"""
mini_blast.py
A simplified, from-scratch version of BLAST's seed-and-extend strategy
(Altschul et al., 1990) — NOT a reimplementation of real BLAST, but
enough to demonstrate the core idea and compare its behavior to the
full dynamic-programming methods above.

Algorithm:
  1. SEED: index all short "words" (default length 3) in sequence 2.
  2. MATCH: for every word in sequence 1, look up exact-word hits in
     sequence 2 (real BLAST also allows near-matches via a scoring
     threshold — noted as a future direction below).
  3. EXTEND: from each seed hit, extend left and right one residue at
     a time using BLOSUM62 scores, stopping when the running score
     drops more than X below the best score seen so far (an "X-drop"
     rule, same concept real BLAST uses to avoid extending forever
     through noise).
  4. Report the highest-scoring extended hit (the best HSP,
     "high-scoring segment pair") as this tool's similarity score.

This is intentionally ungapped within each extension (no insertions/
deletions), which is the main way it differs from Smith-Waterman —
and is exactly why it's fast: no full O(n*m) matrix is ever built.
"""

from bio_utils import blosum_score

WORD_SIZE = 3
X_DROP = 10  # stop extending once score falls this far below the best-so-far


def _build_word_index(seq, word_size):
    index = {}
    for i in range(len(seq) - word_size + 1):
        word = seq[i:i + word_size]
        index.setdefault(word, []).append(i)
    return index


def _extend(seq1, seq2, i, j, word_size, x_drop):
    """Ungapped extension from a seed at seq1[i:i+w] == seq2[j:j+w]."""
    # score of the seed itself
    seed_score = sum(blosum_score(seq1[i + k], seq2[j + k]) for k in range(word_size))

    # extend right
    right_score = 0
    best_right = 0
    best_right_ext = 0
    a, b = i + word_size, j + word_size
    while a < len(seq1) and b < len(seq2):
        right_score += blosum_score(seq1[a], seq2[b])
        if right_score > best_right:
            best_right = right_score
            best_right_ext = a - (i + word_size) + 1
        if best_right - right_score > x_drop:
            break
        a += 1
        b += 1

    # extend left
    left_score = 0
    best_left = 0
    best_left_ext = 0
    a, b = i - 1, j - 1
    while a >= 0 and b >= 0:
        left_score += blosum_score(seq1[a], seq2[b])
        if left_score > best_left:
            best_left = left_score
            best_left_ext = i - a
        if best_left - left_score > x_drop:
            break
        a -= 1
        b -= 1

    total_score = seed_score + best_left + best_right
    start1 = i - best_left_ext
    end1 = i + word_size + best_right_ext
    start2 = j - best_left_ext
    end2 = j + word_size + best_right_ext

    return total_score, start1, end1, start2, end2


def mini_blast(seq1, seq2, word_size=WORD_SIZE, x_drop=X_DROP):
    index2 = _build_word_index(seq2, word_size)

    best_hsp = None
    for i in range(len(seq1) - word_size + 1):
        word = seq1[i:i + word_size]
        for j in index2.get(word, []):
            total_score, s1, e1, s2, e2 = _extend(seq1, seq2, i, j, word_size, x_drop)
            if best_hsp is None or total_score > best_hsp["score"]:
                best_hsp = {
                    "score": total_score,
                    "seq1_region": seq1[s1:e1],
                    "seq2_region": seq2[s2:e2],
                    "seq1_start": s1,
                    "seq1_end": e1,
                    "seq2_start": s2,
                    "seq2_end": e2,
                }

    if best_hsp is None:
        return {
            "score": 0,
            "seq1_region": "",
            "seq2_region": "",
            "percent_identity": 0.0,
            "aligned_length": 0,
            "note": "No shared word of length {} found — no seed to extend from.".format(word_size),
        }

    r1, r2 = best_hsp["seq1_region"], best_hsp["seq2_region"]
    matches = sum(1 for x, y in zip(r1, r2) if x == y)
    aligned_length = max(len(r1), len(r2))
    percent_identity = 100.0 * matches / aligned_length if aligned_length else 0.0

    best_hsp["percent_identity"] = percent_identity
    best_hsp["aligned_length"] = aligned_length
    return best_hsp


if __name__ == "__main__":
    result = mini_blast("HEAGAWGHEE", "PAWHEAE")
    print("Best HSP score:", result["score"])
    print(result["seq1_region"])
    print(result["seq2_region"])
    print(f"% identity: {result['percent_identity']:.1f}")
