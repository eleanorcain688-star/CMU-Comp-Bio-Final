"""
smith_waterman.py
Local pairwise alignment (Smith & Waterman, 1981).

Finds the single best-matching SUBREGION between two sequences,
rather than forcing a full end-to-end alignment. This is the right
tool when you expect only a domain or motif to be conserved (e.g.
comparing a structurally-similar but functionally-different protein
pair), because unrelated flanking regions don't get penalized.
"""

from bio_utils import blosum_score

GAP_PENALTY = -8


def smith_waterman(seq1, seq2, gap_penalty=GAP_PENALTY):
    n, m = len(seq1), len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    best_score = 0
    best_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + blosum_score(seq1[i - 1], seq2[j - 1])
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            best = max(0, diag, up, left)  # the "0" floor is what makes this local
            score[i][j] = best

            if best == 0:
                trace[i][j] = None
            elif best == diag:
                trace[i][j] = 'D'
            elif best == up:
                trace[i][j] = 'U'
            else:
                trace[i][j] = 'L'

            if best > best_score:
                best_score = best
                best_pos = (i, j)

    # Traceback from the highest-scoring cell until we hit a 0
    aligned1, aligned2 = [], []
    i, j = best_pos
    while i > 0 and j > 0 and score[i][j] != 0:
        move = trace[i][j]
        if move == 'D':
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif move == 'U':
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        elif move == 'L':
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1
        else:
            break

    aligned1.reverse()
    aligned2.reverse()
    a1, a2 = "".join(aligned1), "".join(aligned2)

    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
    aligned_positions = len(a1)
    percent_identity = 100.0 * matches / aligned_positions if aligned_positions else 0.0

    return {
        "score": best_score,
        "aligned_seq1": a1,
        "aligned_seq2": a2,
        "percent_identity": percent_identity,
        "aligned_length": aligned_positions,
    }


if __name__ == "__main__":
    result = smith_waterman("HEAGAWGHEE", "PAWHEAE")
    print("Score:", result["score"])
    print(result["aligned_seq1"])
    print(result["aligned_seq2"])
    print(f"% identity: {result['percent_identity']:.1f}")
