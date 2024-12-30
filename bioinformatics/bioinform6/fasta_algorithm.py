import re
from tqdm import tqdm
import blosum as bl


class Seq:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence


def parse(file_name):
    res = []
    with open(file_name) as f:
        sequences = re.split(r'\n>', f.read())
        for elem in sequences:
            name, seq = elem.split('\n', maxsplit=1)
            res.append(Seq(name, seq.replace('\n', '')))
    return res


def score_fun(a: str,
              b: str,
              match_score: int = 5,
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score


def smith_waterman(seq1: str, seq2: str, gap_score = -10):
    m, n = len(seq1) + 1, len(seq2) + 1
    matrix = [[0] * n for _ in range(m)]

    for i in range(m):
        matrix[i][0] = i * gap_score
    for j in range(n):
        matrix[0][j] = j * gap_score

    for i in range(1, m):
        for j in range(1, n):
            matrix[i][j] = max(matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]),
                               matrix[i - 1][j] + gap_score,
                               matrix[i][j - 1] + gap_score)
    score = gap_score * max(n, m)
    i, j = -1, -1
    for k in range(1, n):
        for l in range(1, m):
            if matrix[k][l] > score:
                score = matrix[k][l]
                i, j = k, l
    aln1 = ""
    aln2 = ""
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]):
            a, b = seq1[i - 1], seq2[j - 1]
            i -= 1
            j -= 1

        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_score:
            a, b = seq1[i - 1], '-'
            i -= 1

        elif j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_score:
            a, b = '-', seq2[j - 1]
            j -= 1

        aln1 += a
        aln2 += b
    return aln1[::-1], aln2[::-1], score


def group_up(seq, size=2):
    groups = {}
    for ind in range(len(seq) - size + 1):
        if seq[ind:ind + size] in groups:
            groups[seq[ind: ind + size]].append(ind)
        else:
            groups[seq[ind:ind + size]] = [ind]
    return groups


if __name__ == "__main__":
    fasta_path = 'uniprot_sprot.fasta'
    search_path = 'test.txt'
    gap_score = -10
    filter = 50

    records = parse(fasta_path)
    search = parse(search_path)[1]
    print(search.name, search.sequence)
    search_groups = group_up(search.sequence)

    for i in tqdm(range(len(records))):
        record = records[i]

        record_group = group_up(record.sequence)
        diagonal = {}
        for item_a in search_groups:
            if item_a in record_group:
                for start_a in search_groups[item_a]:
                    for start_b in record_group[item_a]:
                        d = start_b - start_a
                        if d not in diagonal:
                            diagonal[d] = {
                                'match': 1,
                                'interval': [[start_a, start_a]]
                            }
                        else:
                            diagonal[d]['match'] += 1
                            dop_var = 0
                            for elem in diagonal[d]['interval']:
                                if elem[0] == start_a + 1:
                                    elem[0] -= 1
                                    dop_var += 1
                                elif elem[1] == start_a - 1:
                                    elem[1] += 1
                                    dop_var += 1
                            if dop_var == 0:
                                diagonal[d]['interval'].append([start_a, start_a])

        diagonals = {key: diagonal[key] for key, value in diagonal.items() if value['match'] >= 10}
        dop = {}
        for key in diagonals:
            for subdiag in diagonals[key]['interval']:
                dop[(subdiag[0], key + subdiag[0])] = {
                    'length': subdiag[1] - subdiag[0] + 1,
                    'score': sum(bl.BLOSUM(62).get(str(search.sequence[subdiag[0] + p] +
                                                       record.sequence[key + subdiag[0] + p]))
                                 for p in range(subdiag[1] - subdiag[0] + 1))
                }
        diagonals = dop
        diagonals = {key: diagonals[key] for key in diagonals if diagonals[key]['score'] > filter}
        
        res = {}
        for diagonal_1 in list(diagonals.keys()):
            count = 0
            for diagonal_2 in list(diagonals.keys()):
                if diagonal_1 != diagonal_2:
                    dist_1 = diagonal_2[0] - diagonal_1[0] + diagonals[diagonal_1]['length'] - 1
                    dist_2 = diagonal_2[1] - diagonal_1[1] + diagonals[diagonal_1]['length'] - 1

                    if dist_1 >= 0 and dist_2 >= 0:
                        score = (dist_1 + dist_2) * gap_score + \
                                diagonals[diagonal_1]['score'] + \
                                diagonals[diagonal_2]['score']
                        if score > diagonals[diagonal_1]['score']:
                            diagonals[(diagonal_1[0], diagonal_1[1])] = {
                                'length': diagonals[diagonal_1]['length'] + dist_1 + diagonals[diagonal_2]['length'],
                                'score': score
                            }
                            count += 1
            if count == 0:
                res[diagonal_1] = diagonals[diagonal_1]
            diagonals = {key: diagonals[key] for key in diagonals if key != diagonal_1}
        diagonals = res

        result = [smith_waterman(search.sequence[seq[0]: (seq[0] + diagonals[seq]['length'])],
                     record.sequence[seq[1]: (seq[1] + diagonals[seq]['length'])], gap_score)
                  for seq in diagonals]
        score = max(s[2] for s in result) if [s[2] for s in result] else -1

        if score > 0:
            print(score, record.name)