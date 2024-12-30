import re
import copy
from tqdm import tqdm
from enum import IntEnum
import numpy as np
import blosum as bl


def parse_fasta(filename, n=3):
    with open(filename, 'r') as f:
        if n == -1:
            seqs = re.split(r'>.*?\n', f.read())[1:]
            print(seqs)
        else:
            seqs = re.split(r'>', f.read())[1:n + 1]
            name, dd = seqs[0].split('\n', maxsplit=1)
            print(name)
            print(dd)


    return [s.replace('\n', '') for s in seqs]


class Trace(IntEnum):
    STOP, LEFT, UP, DIAGONAL = 0, 1, 2, 3


def smith_waterman(seq1, seq2):
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)
    tracing_matrix = np.zeros(shape=(row, col), dtype=int)

    max_score = -1
    max_index = (-1, -1)

    for i in range(1, row):
        for j in range(1, col):
            match_value = bl.BLOSUM(62).get(str(seq1[i - 1] + seq2[j - 1]))
            diagonal_score = matrix[i - 1, j - 1] + match_value

            vertical_score = matrix[i - 1, j] + gap

            horizontal_score = matrix[i, j - 1] + gap

            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)

            if matrix[i, j] == 0:
                tracing_matrix[i, j] = Trace.STOP

            elif matrix[i, j] == horizontal_score:
                tracing_matrix[i, j] = Trace.LEFT

            elif matrix[i, j] == vertical_score:
                tracing_matrix[i, j] = Trace.UP

            elif matrix[i, j] == diagonal_score:
                tracing_matrix[i, j] = Trace.DIAGONAL

            if matrix[i, j] >= max_score:
                max_index = (i, j)
                max_score = matrix[i, j]

    aligned_seq1 = ""
    aligned_seq2 = ""
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""
    (max_i, max_j) = max_index

    while tracing_matrix[max_i, max_j] != Trace.STOP:
        if tracing_matrix[max_i, max_j] == Trace.DIAGONAL:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1

        elif tracing_matrix[max_i, max_j] == Trace.UP:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            max_i = max_i - 1

        elif tracing_matrix[max_i, max_j] == Trace.LEFT:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            max_j = max_j - 1

        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return aligned_seq1, aligned_seq2, max_score









class Index:
    def __init__(self, sequence, k=2):
        self.data, self.k, self.s = {}, k, sequence
        for i in range(len(sequence) - k + 1):
            if sequence[i:i + k] not in self.data:
                self.data[sequence[i:i + k]] = []
            self.data[sequence[i:i + k]].append(i)


class Diff:
    def __init__(self, i_index, db_index):
        self.data, self.raw_lines = {}, {}
        for window in list(i_index.data.keys()):
            for i_pos in i_index.data.get(window, []):
                for db_pos in db_index.data.get(window, []):
                    diff = i_pos - db_pos
                    if diff not in self.data:
                        self.data[diff] = 0
                        self.raw_lines[diff] = []
                    self.raw_lines[diff].append(Line((i_pos, db_pos),
                                                     (i_pos + i_index.k, db_pos + db_index.k),
                                                     i_index.s, db_index.s))
                    self.data[diff] += 1

    def merge_lines(self):
        self.lines = []
        for k, lines in self.raw_lines.items():
            lines.sort(key=lambda x: x.begin)
            cur_lines = [lines[0]]
            for i in range(1, len(lines)):
                if lines[i].begin < cur_lines[-1].end:
                    cur_lines[-1].end = lines[i].end
                else:
                    cur_lines.append(lines[i])
            self.lines += cur_lines
        self.lines = list(filter(lambda line: line.end[0] - line.begin[0] > min_diag_len, self.lines))

    def get_rects(self):
        self.rects = [Rect(line) for line in self.lines]
        is_continue = True
        while is_continue:
            is_continue = False
            for rect in self.rects:
                if rect.is_used:
                    continue
                for line in self.lines:
                    new_rect = rect.add_line(line)
                    if new_rect:
                        is_continue = True
                        self.rects.append(new_rect)


class Line:
    def __init__(self, begin, end, i_s, db_s):
        self.begin, self.end = begin, end
        self.i_s, self.db_s = i_s, db_s
        self.k = abs(begin[0] - begin[1])
        self.score = sum([bl.BLOSUM(62).get(str(i_s[begin[0] + i] + db_s[begin[1] + i]))
                          for i in range(min(end[0] - begin[0], end[1] - begin[1]))])

    def __eq__(self, other):
        return self.begin == other.begin and self.end == other.end

    def gap_score(self, line):
        return gap * abs(line.k - self.k)


class Rect:
    def __init__(self, line):
        self.lines, self.is_used = [line], False

    @property
    def end(self):
        return self.lines[-1].end

    def is_line_valid(self, line):
        return (line.begin[0] < self.end[0] or line.begin[1] < self.end[1]) and \
            ((line.end[0] >= self.end[0] and line.end[1] > self.end[1]) or
             (line.end[0] > self.end[0] and line.end[1] >= self.end[1]))

    def add_line(self, line):
        if line in self.lines:
            return None
        if self.is_line_valid(line):
            delta = abs(line.k - self.lines[-1].k)
            if self.lines[-1].k == line.k:
                return None
            if self.lines[-1].k < line.k:
                new_begin = (self.lines[-1].end[0], self.lines[-1].end[1] + delta)
            else:
                new_begin = (self.lines[-1].end[0] + delta, self.lines[-1].end[1])
            new_line = Line(new_begin, line.end, line.i_s, line.db_s)
            if new_line.score + new_line.gap_score(self.lines[-1]) > 0:
                new_rect = copy.deepcopy(self)
                new_rect.lines.append(new_line)
                self.is_used = True
                new_rect.is_used = False
                return new_rect
        return None


def get_compare_data(query, seqs):
    i_index = Index(query)
    ans = [-float('inf')]
    for seq in tqdm(seqs):
        db_index = Index(seq)
        diff = Diff(i_index, db_index)
        diff.merge_lines()
        diff.get_rects()

        for rect in diff.rects:
            seq1 = query[rect.lines[-1].begin[0]:rect.lines[-1].end[0]]
            seq2 = seq[rect.lines[-1].begin[1]:rect.lines[-1].end[1]]
            ans_tmp = (seq,) + smith_waterman(seq1, seq2)
            if ans_tmp[-1] > ans[-1]:
                ans = copy.copy(ans_tmp)
    return ans


gap = -10
min_diag_len = 5
seqs = parse_fasta('uniprot_sprot.fasta', n=5000)
query = 'EKGLIVGHFSGIKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD'
ans = get_compare_data(query, seqs)
print('seq1:', ans[1])
print('seq2:', ans[2])
print('score:', ans[3])
print('full answer:', ans[0])
