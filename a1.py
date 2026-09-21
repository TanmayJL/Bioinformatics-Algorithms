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

    # Print the optimal score
    print(f"Optimal alignment score: {optimal_score}")

if __name__ == "__main__":
    main()