from typing import Callable, Tuple

DEBUG = False


def score_fun(a: str,
              b: str,
              match_score: int = 5,
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score


def pick_max(a, b, c):
    if a >= b and a >= c:
        return a, '\\'
    if c >= b and c >= a:
        return c, '<'
    if b >= a and b >= c:
        return b, '^'


class NeedlemanWunsch(object):

    def __init__(self, seq_a, seq_b, gap_open=-10, gap_extend=-1):
        self.seq_a = seq_a
        self.seq_b = seq_b

        self.gap_open = gap_open
        self.gap_extend = gap_extend

        self.min_value = -float("inf")

        self.rows = len(seq_a) + 1
        self.cols = len(seq_b) + 1

        # Algorithm is based on three matrices X (gap in X), Y (gap in Y), M (match between X and Y) #
        self.I = [[0 for _ in range(self.cols)] for _ in range(self.rows)]
        self.M = [[0 for _ in range(self.cols)] for _ in range(self.rows)]
        self.D = [[0 for _ in range(self.cols)] for _ in range(self.rows)]
        self.res = [[0 for _ in range(self.cols)] for _ in range(self.rows)]

    def compute_matrix(self):
        self.M[0][0], self.I[0][0], self.D[0][0] = 0, float('-inf'), float('-inf')

        for i in range(1, self.rows):
            self.M[i][0], self.D[i][0] = float('-inf'), float('-inf')
            self.I[i][0] = self.gap_open + (i - 1) * self.gap_extend

        for j in range(1, self.cols):
            self.M[0][j], self.I[0][j] = float('-inf'), float('-inf')
            self.D[0][j] = self.gap_open + (j - 1) * self.gap_extend

        for i in range(1, self.rows):
            for j in range(1, self.cols):
                self.M[i][j] = max(self.M[i - 1][j - 1] + score_fun(self.seq_a[i - 1], self.seq_b[j - 1]),
                                   self.I[i - 1][j - 1] + score_fun(self.seq_a[i - 1], self.seq_b[j - 1]),
                                   self.D[i - 1][j - 1] + score_fun(self.seq_a[i - 1], self.seq_b[j - 1]))

                self.I[i][j] = max(self.I[i][j - 1] + self.gap_extend,
                                   self.M[i][j - 1] + self.gap_open)

                self.D[i][j] = max(self.D[i - 1][j] + self.gap_extend,
                                   self.M[i - 1][j] + self.gap_open)

        for i in range(self.rows):
            for j in range(self.cols):
                self.res[i][j] = pick_max(self.M[i][j], self.I[i][j], self.D[i][j])

    def get_align(self):
        aln1 = ''
        aln2 = ''
        i, j = self.rows - 1, self.cols - 1
        is_extended = False

        while i > 0 or j > 0:
            a, b = '-', '-'
            # (A, B)
            if not is_extended and i > 0 and j > 0 and self.res[i][j][0] == self.res[i - 1][j - 1][0] + score_fun(
                    self.seq_a[i - 1], self.seq_b[j - 1]):
                a = self.seq_a[i - 1]
                b = self.seq_b[j - 1]
                i -= 1
                j -= 1

            # (A, -)
            elif i > 0 and self.res[i][j][0] == self.res[i - 1][j][0] + self.gap_open:
                is_extended = False
                a = self.seq_a[i - 1]
                i -= 1

            # (-, A)
            elif j > 0 and self.res[i][j][0] == self.res[i][j - 1][0] + self.gap_open:
                is_extended = False
                b = self.seq_b[j - 1]
                j -= 1

            elif i > 0 and self.res[i][j][0] == self.res[i - 1][j][0] + self.gap_extend:
                is_extended = True
                a = self.seq_a[i - 1]
                i -= 1

            elif j > 0 and self.res[i][j][0] == self.res[i][j - 1][0] + self.gap_extend:
                is_extended = True
                b = self.seq_b[j - 1]
                j -= 1

            aln1 += a
            aln2 += b
        return aln1, aln2


def needleman_wunsch_affine(seq1: str,
                            seq2: str,
                            score_fun: Callable = score_fun,
                            gap_open: int = -10,
                            gap_extend: int = -1) -> Tuple[str, str, int]:
    '''
    Inputs:
    seq1 - first sequence
    seq2 - second sequence
    score_fun - function that takes two characters and returns score
    gap_open - gap open penalty
    gap_extend - gap extend penalty
    Outputs:
    aln1 - first aligned sequence
    aln2 - second aligned sequence
    score - score of the alignment
    '''
    needlemanWunschAffine = NeedlemanWunsch(seq1, seq2, gap_open, gap_extend)

    needlemanWunschAffine.compute_matrix()

    aln1, aln2 = needlemanWunschAffine.get_align()

    return aln1[::-1], aln2[::-1], needlemanWunschAffine.res[-1][-1][0]


def print_array(matrix: list):
    for row in matrix:
        for element in row:
            print(f"{element:6}", end="")
        print()


def main():
    aln1, aln2, score = needleman_wunsch_affine("ACGT", "TAGT", gap_open=-10, gap_extend=-1)
    print(f'str 1: {aln1}')
    print(f'str 2: {aln2}')
    print(f'score: {score}')


if __name__ == "__main__":
    main()
