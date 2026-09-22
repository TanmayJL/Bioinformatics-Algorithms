def needleman_wunsch(dp, seq1, seq2, match=2, miss=-1, d=-2):
    """
    Reads 2 sequences from sequences.txt.
    Returns optimal score (int)
    """
    m = len(seq1)
    n = len(seq2)
    dp[0][0] = 0

    # initialize sequence aligned to gaps
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + d
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + d

    # compute full dp table
    for i in range(2, m + 1):
        for j in range(2, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = max(dp[i - 1][j - 1] + match, dp[i - 1][j] + d, dp[i][j - 1] + d)
            else:
                dp[i][j] = max(dp[i - 1][j - 1] + miss, dp[i - 1][j] + d, dp[i][j - 1] + d)
    return dp[m][n]

def print_dp_table(dp):
    for row in dp:
        print(" ".join(map(str, row)))

def get_optimal_alignment(dp, seq1, seq2, match=2, miss=-1, d=-2):
    """
    Backtrack through the dp table to find the optimal alignment of the
    two sequences. Returns a tuple of the aligned sequences
    """
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len(seq1)
    j = len(seq2)

    
    while i > 0 or j > 0:
        inc = match if seq1[i - 1] == seq2[j - 1] else miss
        # derived from a match/mismatch
        if i > 0 and j > 0 and dp[i][j] == (dp[i - 1][j - 1] + inc):
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += seq2[j - 1]
            i -= 1
            j -= 1
        # derived from a deletion (gap in seq2)
        elif i > 0 and dp[i][j] == (dp[i - 1][j] - d):
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += "-"
            i -= 1
        else: # derived from an insertion (gap in seq1)
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j - 1]
            j -= 1
    return aligned_seq1[::-1], aligned_seq2[::-1]

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
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Compute optimal score
    optimal_score = needleman_wunsch(dp, seq1, seq2)
    print(f"Optimal alignment score: {optimal_score}")
    print("Dynamic Programming Table:")
    print_dp_table(dp)

    # Get optimal alignment
    aligned_seq1, aligned_seq2 = get_optimal_alignment(dp, seq1, seq2)
    print("Optimal Alignment:")
    print(aligned_seq1)
    print(aligned_seq2)


if __name__ == "__main__":
    main()