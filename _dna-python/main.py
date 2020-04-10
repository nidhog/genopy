# Introduction to Dynamic Programming for Bioinformatics
# Comparing DNA strands
import numpy as np

# Global Alignment
# 'global': entire sequence
# 'alignment': maximise score of match/mismatch/gaps between strands
# Pairwise: create gaps for better alignment
# Note there are other ways

DNA_nuc = {'A': 0b00,
           'C': 0b11,
           'G': 0b00,
           'T': 0b01}

# Needlman-Wunsch algorithm for similarity scoring, alignment
# not good:
"""
Sequence 1 ==> G T C C A T A C A
Sequence 2 ==> T C A T A T C A G
"""
# +1 for each match
# -1 for mismatch
# -1 for gap (insert/delete)
# Above score: 2-7-0=-5
"""
Sequence 1 ==> G T C C A T A - C A -
Sequence 2 ==> - T C - A T A T C A G
"""
# Score: 7-4=3
# insertion/deletion makes sense in biology because of mutations
"""
Matrix
 G T C C ...
T -1 . . ...
C  .       .
A    .     .
.      .   .
.        . .
.  . . . ...
"""
X = 'GTCCATACA'
Y = 'TCATATCAG'

PENAL = -1
RWARD = +1


def get_score(n1, n2, penalty=PENAL, reward=RWARD):
    """NW algorithm"""
    if n1 == n2:
        return reward
    return penalty


def get_align(x=X, y=Y):
    mat_size = (len(x) + 1, len(y) + 1)
    score_matrix_xy = np.ndarray(mat_size)
    for j in range(len(y) + 1):
        score_matrix_xy[0, j] = PENAL*j
    for i in range(len(x) + 1):
        score_matrix_xy[i, 0] = PENAL*i
    # start from 1 or X, Y won't have indices
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            its_a_match = score_matrix_xy[i-1, j-1] + get_score(x[i - 1], y[i - 1])
            del_s = score_matrix_xy[i-1, j] + PENAL
            add_s = score_matrix_xy[i, j-1] + PENAL

            score_matrix_xy[i,j] = max([its_a_match, del_s, add_s])

    n_x, n_y = len(x), len(y)
    alignment = {0:'', 1:''}
    while n_x or n_y:
        score = score_matrix_xy[n_x, n_y]
        left_score = score_matrix_xy[n_x-1, n_y]

        if n_x and n_y and x[n_x - 1] == y[n_y - 1]:
            alignment[0]=(x[n_x-1])+alignment[0]
            alignment[1]=(y[n_x-1])+alignment[1]
            n_x -= 1
            n_y -= 1
        elif n_x>0 and score == left_score + PENAL:
            alignment[0]=(x[n_x - 1])+alignment[0]
            alignment[1]='-'+alignment[1]
            n_x -= 1
        else:
            alignment[0]= '-'+alignment[0]
            alignment[1]= y[n_y - 1]+alignment[1]
            n_y -= 1
    return alignment


if __name__=="__main__":
    al=get_align(X, Y)
    [print('{}: {}'.format(x, al[x])) for x in al]





