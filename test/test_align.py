from opt_align import align
from opt_align import io

def test_align_score():
    score_mat = io.read_score("score_matrices/BLOSUM50")
    s1 = "MEANLY"
    s2 = "PLEASANTLY"
    assert int(align.local_alignment(s1,s2, score_mat, gap=10, ext=2)[0]) == 15

def test_alignment():
    score_mat = io.read_score("score_matrices/BLOSUM62")
    s1 = "PRTEINS"
    s2 = "PRTWPSEIN"
    score, back, max_i, max_j, middle = align.local_alignment(s1, s2, score_mat, gap=11, ext=1)
    s1, s2 = align.output_lcs(back,max_i, max_j, s1, s2, "","",middle)
    assert [s1,s2] == ["PRT---EIN","PRTWPSEIN"]
