from bio.format import *
import argparse
import sys


class NeedlemanWunsch(object):

    def __init__(self, seq_a, seq_b, substution_matrix=None, position_dict=None, gap_penalty=-5, func_match=None):
        self.seq_a = seq_a
        self.seq_b = seq_b

        self.substution_matrix = substution_matrix
        self.position_dict = position_dict

        self.gap_penalty = gap_penalty

        self.func_match = func_match

        self.rows = len(seq_a) + 1
        self.cols = len(seq_b) + 1

        self.matrix = [[0 for _ in range(self.cols)] for _ in range(self.rows)]

    @property
    def score(self):
        return self.matrix[self.rows - 1][self.cols - 1]

    def s(self, i, j):
        if self.substution_matrix is not None:
            index_i = self.position_dict[self.seq_a[i - 1]]
            index_j = self.position_dict[self.seq_b[j - 1]]
            return self.substution_matrix[index_i][index_j]
        elif self.func_match is not None:
            return self.func_match(self.seq_a[i - 1], self.seq_b[j - 1])
        else:
            return 1 if self.seq_a[i - 1] == self.seq_b[j - 1] else -1

    def compute_matrix(self):

        self.matrix[0][0] = 0
        for i in range(self.rows):
            self.matrix[i][0] = i * self.gap_penalty
        for j in range(self.cols):
            self.matrix[0][j] = j * self.gap_penalty
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                self.matrix[i][j] = max(self.matrix[i - 1][j - 1] + self.s(i, j),  # Diag
                                        self.matrix[i - 1][j] + self.gap_penalty,  # Up
                                        self.matrix[i][j - 1] + self.gap_penalty)  # Left


    def get_align(self):
        seq_a_align, seq_b_align = "", ""

        i = self.rows - 1
        j = self.cols - 1

        while i > 0 and j > 0:
            if self.matrix[i][j] == (self.matrix[i - 1][j - 1] + self.s(i, j)):
                seq_a_align = self.seq_a[i - 1] + seq_a_align
                seq_b_align = self.seq_b[j - 1] + seq_b_align
                i -= 1
                j -= 1
            elif self.matrix[i][j] == self.matrix[i - 1][j] + self.gap_penalty:
                seq_a_align = self.seq_a[i - 1] + seq_a_align
                seq_b_align = "-" + seq_b_align
                i -= 1
            elif self.matrix[i][j] == self.matrix[i][j - 1] + self.gap_penalty:
                seq_a_align = "-" + seq_a_align
                seq_b_align = self.seq_b[j - 1] + seq_b_align
                j -= 1
        while i > 0:
            seq_a_align = self.seq_a[i - 1] + seq_a_align
            seq_b_align = "-" + seq_b_align
            i -= 1
        while j > 0:
            seq_a_align = "-" + seq_a_align
            seq_b_align = self.seq_b[j - 1] + seq_b_align
            j -= 1

        return seq_a_align, seq_b_align


def new_parser():
    parser = argparse.ArgumentParser(description="Description of using script")
    parser.add_argument('-m', dest="match", type=int, help="Value of match")
    parser.add_argument('-mm', dest="mismatch", type=int, help="Value of mismatch")
    parser.add_argument('-f', dest="matrix", type=str, help="File, which contains table of match/mismatch")
    parser.add_argument('-g', dest="gap", type=int, help="Value of gap penalty")
    parser.add_argument("file", type=str, help="Name of fasta file, which contains sequences")
    return parser


if __name__ == "__main__":
    arg_parse = new_parser()
    args = arg_parse.parse_args()

    if args.file is not None:

        seq = read_file(args.file)
        a = seq[0].seq
        b = seq[1].seq

        if args.matrix is not None:
            matrix, pos_dict = load_matrix(args.matrix)

            gap = -10 if args.gap is None else args.gap

            needlemanWinsch = NeedlemanWunsch(a, b, substution_matrix=matrix, position_dict=pos_dict, gap_penalty=gap)
            needlemanWinsch.compute_matrix()

            output(needlemanWinsch.get_align(), needlemanWinsch.score)

        elif args.match is not None and args.mismatch is not None:
            func_match = lambda x, y: args.match if x == y else args.mismatch

            gap = -10 if args.gap is None else args.gap

            needlemanWinsch = NeedlemanWunsch(a, b, func_match=func_match, gap_penalty=gap)
            needlemanWinsch.compute_matrix()

            output(needlemanWinsch.get_align(), needlemanWinsch.score)
        else:
            sys.exit(1)

    else:
        sys.exit(1)
