def score_fun(a: str,
              b: str,
              match_score: int = 0,
              mismatch_score: int = 1) -> int:
    return match_score if a == b else mismatch_score


def needleman_wunsch(seq1: str, seq2: str, gap_score: int = 2):
    m, n = len(seq1) + 1, len(seq2) + 1
    matrix = [[0] * n for _ in range(m)]

    for i in range(m):
        matrix[i][0] = i * gap_score
    for j in range(n):
        matrix[0][j] = j * gap_score

    for i in range(1, m):
        for j in range(1, n):
            matrix[i][j] = min(matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]),
                               matrix[i - 1][j] + gap_score,
                               matrix[i][j - 1] + gap_score)
    return matrix[-1][-1]


def dict_min(dict):
    min = float('inf')
    min_i = None
    for i in dict:
        splitted = i.split('|')
        if splitted[0] != splitted[1] and dict[i] < min:
            min = dict[i]
            min_i = i

    return min_i


def get_dict(seqs, seq_names):
    dist = {}

    for i in range(len(seqs)):
        for j in range(len(seqs)):
            res = needleman_wunsch(seqs[i], seqs[j])
            dist[seq_names[i] + "|" + seq_names[j]] = res
    return dist


def get_from_dict(dict, str1, str2):
    if dict.get(str1 + '|' + str2) is not None:
        return dict[str1 + '|' + str2]
    return dict[str2 + '|' + str1]


def parse_fasta(path):
    clusters = []
    strings = []
    f = open(path, 'r')
    data = f.read()
    if path.count('large'):
        splitted = data.split('\n\n')
        for str in splitted:
            clusters.append(str.split(' ')[0][1:])
            strings.append(''.join(str.split('\n')[1:]))
        return strings, clusters
    else:
        splitted = data.split('\n')
        for i in range(0, len(splitted), 2):
            clusters.append(splitted[i][1:])
            strings.append(splitted[i + 1])
        return strings, clusters


def run_seq(str):
    spl = str.split('|')
    left = spl[0]
    right = spl[1]
    print("{} ---- {}".format(left, right))


def concat_dict(dict, min_i, answers):

    splitted = min_i.split('|')
    run_seq(min_i)

    if len(splitted[0]) > 2 and len(splitted[1]) > 2:
        answers.append('(' + '(' + str(splitted[0]) + ')' + ',' + '(' + str(splitted[1]) + ')' + ')' + ' => ' + str(dict[min_i] / 2))
    elif len(splitted[0]) > 2:
        answers.append('(' + '(' + str(splitted[0]) + ')' + ',' + str(splitted[1]) + ')' + ' => ' + str(dict[min_i] / 2))
    elif len(splitted[1]) > 2:
        answers.append('(' + str(splitted[0]) + ',' + '(' + str(splitted[1]) + ')' + ')' + ' => ' + str(dict[min_i] / 2))
    else:
        answers.append('(' + str(splitted[0]) + ',' + str(splitted[1]) + ')' + ' => ' + str(dict[min_i] / 2))

    name_new_cluster = str(splitted[0]) + ',' + str(splitted[1])
    cluster_cardinality[name_new_cluster] = cluster_cardinality[splitted[0]] + cluster_cardinality[splitted[1]]


    new_dict = {splitted[0] + ',' + splitted[1] + '|' + splitted[0] + ',' + splitted[1]: 0}
    for i in dict:
        splitted_i = i.split('|')
        if not splitted_i.count(splitted[0]) and not splitted_i.count(splitted[1]):
            new_dict[i] = dict[i]

    #print(new_dict)
    for cluster in clusters:
        if cluster != splitted[0] and cluster != splitted[1]:
            #print(splitted[0] + ',' + splitted[1] + '|' + cluster)
            new_dict[splitted[0] + ',' + splitted[1] + '|' + cluster] = \
                (get_from_dict(dict, splitted[0], cluster) *
                 cluster_cardinality[splitted[0]] +
                 get_from_dict(dict, splitted[1], cluster) *
                 cluster_cardinality[splitted[1]]) / \
                (cluster_cardinality[splitted[0]] + cluster_cardinality[splitted[1]])
    clusters.remove(splitted[0])
    clusters.remove(splitted[1])
    clusters.append(name_new_cluster)
    #print(new_dict)
    return new_dict, answers


answers = []
strings, clusters = parse_fasta('small.fasta')
dict = get_dict(strings, clusters)
cluster_cardinality = {}
for cluster in clusters:
    cluster_cardinality[cluster] = 1
#print(dict)

while len(dict.keys()) > 1:
    dict, answers = concat_dict(dict, dict_min(dict), answers)

for a in answers:
    print(a)

