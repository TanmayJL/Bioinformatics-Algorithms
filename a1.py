from typing import List, Set, Tuple, Optional
def needleman_wunsch(dp, seq1, seq2, match=2, miss=-1, d=-2, *, pointer_table: Optional[List[List[Optional[List[Tuple[int, int]]]]]]=None):
    """
    Uses dynamic programming to compute the optimal alignment score of two sequences.
    Returns optimal score (int)
    """
    m = len(seq1)
    n = len(seq2)
    dp[0][0] = 0
    if not pointer_table:
        pointer_table = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    # initialize sequence aligned to gaps
    for i in range(1, m + 1):
        pointer_table[i][0] = [(i - 1, 0)]
        dp[i][0] = dp[i - 1][0] + d
        pointer_table[i][0] = [(i - 1, 0)]
    for j in range(1, n + 1):
        pointer_table[0][j] = [(0, j - 1)]
        dp[0][j] = dp[0][j - 1] + d

    # compute full dp table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                pointer_table[i][j] = [(i - 1, j - 1)]
                dp[i][j] = max(dp[i - 1][j - 1] + match, dp[i - 1][j] + d, dp[i][j - 1] + d)
            else:
                pointer_table[i][j] = [(i - 1, j - 1)]
                dp[i][j] = max(dp[i - 1][j - 1] + miss, dp[i - 1][j] + d, dp[i][j - 1] + d)
    return dp[m][n], pointer_table

def print_dp_table(dp):
    for row in dp:
        print(" ".join(map(str, row)))

def print_pointer_table(pointer_table):
    for row in pointer_table:
        print(" ".join(str(cell) if cell is not None else "None" for cell in row))


def get_optimal_alignments_with_pointer_table(dp, seq1, seq2, pointer_table):
    """
    Backtrack through the dp table using the pointer table to find the optimal alignment of the
    two sequences. Start from bottom right corner and follow pointers to top left corners. If
    adjacent cells have the same score, then there are multiple optimal alignments. Recursively
    explore paths with same score.
    
    Returns a list of tuples of the aligned sequences
    """
    alignments = []
    def backtrack(i, j, aligned_seq1, aligned_seq2):
        if i == 0 and j == 0:
            alignments.append((aligned_seq1[::-1], aligned_seq2[::-1]))
            return
        if pointer_table[i][j] is None:
            return
        prev_i, prev_j = pointer_table[i][j]

        # debugging print statements
        print(f"Backtracking from ({i}, {j}) to ({prev_i}, {prev_j})")
        print(f"Current aligned sequences: {aligned_seq1[::-1]}, {aligned_seq2[::-1]}")
        # derived from a (mis)match
        if prev_i == i - 1 and prev_j == j - 1:
            backtrack(prev_i, prev_j, aligned_seq1 + seq1[i - 1], aligned_seq2 + seq2[j - 1])
        # derived from a deletion (gap in seq2)
        if prev_i == i - 1 and prev_j == j:
            backtrack(prev_i, prev_j, aligned_seq1 + seq1[i - 1], aligned_seq2 + "-")
        # derived from an insertion (gap in seq1)
        if prev_i == i and prev_j == j - 1:
            backtrack(prev_i, prev_j, aligned_seq1 + "-", aligned_seq2 + seq2[j - 1])

    backtrack(len(seq1), len(seq2), "", "")
    
    return alignments

def multiple_alignments(dp, seq1, seq2, match=2, miss=-1, d=-2):
    """
    Brack through dp table to find all optimal alignments.
    Returns a list of tuples of the aligned sequences 
    """
    m = len(seq1)
    n = len(seq2)
    alignments = []
    def backtrack(i, j, aligned_seq1, aligned_seq2):
        if i == 0 and j == 0:
            alignments.append((aligned_seq1[::-1], aligned_seq2[::-1]))
            return
        inc = match if seq1[i - 1] == seq2[j - 1] else miss
        # derived from a match/mismatch
        if i > 0 and j > 0 and dp[i][j] == (dp[i - 1][j - 1] + inc):
            backtrack(i - 1, j - 1, aligned_seq1 + seq1[i - 1], aligned_seq2 + seq2[j - 1])
        # derived from a deletion (gap in seq2)
        if i > 0 and dp[i][j] == (dp[i - 1][j] - d):
            backtrack(i - 1, j, aligned_seq1 + seq1[i - 1], aligned_seq2 + "-")
        # derived from an insertion (gap in seq1)
        if j > 0 and dp[i][j] == (dp[i][j - 1] - d):
            backtrack(i, j - 1, aligned_seq1 + "-", aligned_seq2 + seq2[j - 1])
    backtrack(m, n, "", "")
    return alignments

def main():
    """
    Main function to run the Needleman-Wunsch algorithm.
    """
    # Read sequences from file
    with open("sequences.txt", "r") as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()

    # Initialize DP table
    m = len(seq1)
    n = len(seq2)
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    pointer_table: list[list[Optional[tuple[int, int]]]] = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    # Compute optimal score
    optimal_score = needleman_wunsch(dp, seq1, seq2, pointer_table=pointer_table)[0]
    print(f"Optimal alignment score: {optimal_score}")
    print("Dynamic Programming Table:")
    print_dp_table(dp)

    # Get optimal alignment
    aligned_seq1, aligned_seq2 = get_optimal_alignment(dp, seq1, seq2)
    print("Optimal Alignment:")
    print(aligned_seq1)
    print(aligned_seq2)
    print("Pointer Table:")
    print_pointer_table(pointer_table)

    # Get all optimal alignments
    all_alignments = get_optimal_alignments_with_pointer_table(dp, seq1, seq2, pointer_table)
    all_alignments = multiple_alignments(dp, seq1, seq2)
    print("All Optimal Alignments:")
    for alignment in all_alignments:
        print(alignment[0])
        print(alignment[1])
        print()
    print(f"Total number of optimal alignments: {len(all_alignments)}")


if __name__ == "__main__":
    main()