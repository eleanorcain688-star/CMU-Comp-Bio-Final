"""
needleman_wunsch.py
Global pairwise alignment (Needleman & Wunsch, 1970).

Forces the ENTIRE length of both sequences into the alignment, using
a linear gap penalty. Good baseline for "how similar are these two
proteins overall, end to end?" — but it's easily dominated by
insertions/deletions when comparing proteins of different length or
with only a locally-conserved domain.
"""

from bio_utils import blosum_score

GAP_PENALTY = -8  # linear gap cost, applied per gap position


def needleman_wunsch(seq1, seq2, gap_penalty=GAP_PENALTY):
    n, m = len(seq1), len(seq2)

    # DP score matrix: (n+1) x (m+1)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap_penalty
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap_penalty

    # traceback matrix: 'D' diag, 'U' up (gap in seq2), 'L' left (gap in seq1)
    trace = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        trace[i][0] = 'U'
    for j in range(1, m + 1):
        trace[0][j] = 'L'

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + blosum_score(seq1[i - 1], seq2[j - 1])
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                trace[i][j] = 'D'
            elif best == up:
                trace[i][j] = 'U'
            else:
                trace[i][j] = 'L'

    # Traceback to build the alignment strings
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
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
        else:  # 'L'
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()
    a1, a2 = "".join(aligned1), "".join(aligned2)

    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
    aligned_positions = len(a1)
    percent_identity = 100.0 * matches / aligned_positions if aligned_positions else 0.0

    return {
        "score": score[n][m],
        "aligned_seq1": a1,
        "aligned_seq2": a2,
        "percent_identity": percent_identity,
        "aligned_length": aligned_positions,
    }


if __name__ == "__main__":
    # quick smoke test
    result = needleman_wunsch("HEAGAWGHEE", "PAWHEAE")
    print("Score:", result["score"])
    print(result["aligned_seq1"])
    print(result["aligned_seq2"])
    print(f"% identity: {result['percent_identity']:.1f}")
