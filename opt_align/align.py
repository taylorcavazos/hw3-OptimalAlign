# File with methods for aligning two strings 

import numpy as np, pandas as pd, re, os
import matplotlib.pyplot as plt
from .roc import *
import progressbar
from time import sleep
import sys
sys.setrecursionlimit(15000)

def local_alignment(s1, s2, score_mat, gap, ext):
    n, m = len(s1), len(s2)
    max_score, max_i, max_j = float("-inf"), None, None
    # initialize 3-level manhattan graph to track scores of all possible alignments
    # Lower and upper level will track gap extenstions for each string
    # Middle for matches and mismatches
    lower, middle, upper = np.zeros((n+1, m+1)), np.zeros((n+1, m+1)), np.zeros((n+1, m+1))
    lower[:,0], upper[:,0] = float("-inf"), float("-inf")
    lower[0,:], upper[0,:] = float("-inf"), float("-inf")

    # Also keep three arrays for tracking alignments
    backtrack = np.zeros((n+1, m+1))

    for i in range(1, n+1):
        for j in range(1, m+1):
            score = score_mat.loc[s1[i-1], s2[j-1]]
            lower[i,j] = max(lower[i-1,j]-ext, middle[i-1, j]-(gap+ext))

            upper[i,j] = max(upper[i,j-1]-ext, middle[i,j-1]-(gap+ext))

            middle[i,j] = max(0,lower[i,j], upper[i,j],
                              middle[i-1,j-1]+score)
            
            if middle[i][j] == middle[i-1][j-1]+score or middle[i][j] == 0:
                backtrack[i][j]= 1
            elif middle[i,j] == upper[i,j]: backtrack[i,j] = 2
            else: backtrack[i,j] = 0
            
            if middle[i,j] > max_score: 
                max_score = middle[i,j]
                max_i = i
                max_j = j

    return max_score, backtrack, max_i, max_j, middle

def output_lcs(backtrack, i, j, s1, s2, str1, str2, score):
    if score[i][j] == 0:
        return str1, str2
    
    if backtrack[i][j] == 2:
        str1, str2 = output_lcs(backtrack, i, j-1, s1, s2, '-'+str1, s2[j-1]+str2, score)
        
    if backtrack[i][j] == 0:
        str1, str2 = output_lcs(backtrack, i-1, j, s1, s2, s1[i-1]+str1, '-'+str2, score)

    if backtrack[i][j] == 1:
        str1, str2 = output_lcs(backtrack, i-1, j-1, s1, s2, s1[i-1]+str1, s2[j-1]+str2, score)

    return str1, str2

def write_align_output(scores, alignments, outfile):
    out = open(outfile, "w")
    for i in range(0, len(scores)):
        out.write(str(int(scores[i]))+"\t"+"\t".join(alignments[i])+"\n")
    out.close()

def run_alignment(seqs, score_mat, gap=4, ext=3, normalize=False, 
                 out_align=False, outfile=None):
    scores = np.zeros(len(seqs))
    alignments = ['']*len(seqs)
    for i in range(0, len(seqs)):
        s, backtrack, max_i, max_j, middle = local_alignment(seqs[i][0], seqs[i][1], score_mat, gap, ext)
        if normalize == False:
            scores[i] = s
        else:
            scores[i] = (s/min(len(seqs[i][0]), len(seqs[i][1])))
        if out_align == True:
            align_s = output_lcs(backtrack, max_i, max_j,seqs[i][0], seqs[i][1],"","", middle)
            alignments[i] = align_s
    if out_align == True:
        write_align_output(scores, alignments, outfile)
    return scores

def output_align_scores_mats(pos_seqs, neg_seqs, matrix_dict, gap=4, ext=3):
    pos_scores, neg_scores = {}, {}
    for mat, scores in matrix_dict.items():
        pos_scores[mat] = run_alignment(pos_seqs, scores, gap, ext)
        neg_scores[mat] = run_alignment(neg_seqs, scores, gap, ext)
    return pos_scores, neg_scores

def find_best_align_penalties(pos,neg,gap_range, ext_range, score_mat):
    bar = progressbar.ProgressBar(maxval=gap_range*ext_range, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    
    f = open("part1_1.txt", "w")
    best_fpr, best_gap, best_ext = float("inf"), None, None
    
    bar.start()
    i = 0
    for gap in range(1, gap_range+1):
        for ext in range(1, ext_range+1):
            scores_pos = np.array(run_alignment(pos, score_mat, gap, ext))
            scores_neg = np.array(run_alignment(neg, score_mat, gap, ext))
        
            thresh = calc_thresh(0.7, scores_pos)
            fpr = len(scores_neg[scores_neg >= thresh])/len(scores_neg)
            if fpr < best_fpr:
                best_fpr, best_gap, best_ext = fpr, gap, ext
            bar.update(i+1)
            i+=1
            sleep(0.1)
    bar.finish()
    f.write("fpr\tgap\text\n"+str(best_fpr)+"\t"+str(best_gap)+"\t"+str(best_ext))
    f.close()

def find_best_scoring_matrix(pos_scores_dict, neg_scores_dict, tpr=0.7):
    best_fpr, best_mat = float("inf"), None
    mat_fprs = []
    for mat, scores in pos_scores_dict.items():
        thresh = calc_thresh(tpr, scores)
        scores_neg = neg_scores_dict.get(mat)
        fpr = len(scores_neg[scores_neg >= thresh])/len(scores_neg)
        mat_fprs.append(str(fpr))
        if fpr < best_fpr:
            best_fpr, best_mat = fpr, mat
    out = open("part1_2a.txt", "w")
    out.write("\t".join(pos_scores_dict.keys())+"\n")
    out.write("\t".join(mat_fprs))
    out.close()
    return best_mat
