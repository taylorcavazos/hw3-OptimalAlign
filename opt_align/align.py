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
    # Lower and upper level will track gaps for each string
    # Middle for matches and mismatches
    lower, middle, upper = np.zeros((n+1, m+1)), np.zeros((n+1, m+1)), np.zeros((n+1, m+1))
    # initialize first row and column to -infinity because gap must be opened before it can be extended
    lower[:,0], upper[:,0] = float("-inf"), float("-inf")
    lower[0,:], upper[0,:] = float("-inf"), float("-inf")

    # Also keep arrays for tracking alignments
    backtrack = np.zeros((n+1, m+1))
    # Loop through score matrix 
    for i in range(1, n+1):
        for j in range(1, m+1):
            # get score for string1 and string2 character
            score = score_mat.loc[s1[i-1], s2[j-1]]
            # either extend gap or open gap for string1
            lower[i,j] = max(lower[i-1,j]-ext, middle[i-1, j]-(gap+ext))
            # either extend gap or open gap for string2
            upper[i,j] = max(upper[i,j-1]-ext, middle[i,j-1]-(gap+ext))
            # choose to have a gap, match, or end alignment (0)
            middle[i,j] = max(0,lower[i,j], upper[i,j],
                              middle[i-1,j-1]+score)
            # fill backtrack matrix
            if middle[i][j] == middle[i-1][j-1]+score or middle[i][j] == 0:
                backtrack[i][j]= 1 # match
            elif middle[i,j] == upper[i,j]: backtrack[i,j] = 2 # gap in string2
            else: backtrack[i,j] = 0 # gap in string1
            # maintain position of max score to start backtrack
            if middle[i,j] > max_score: 
                max_score = middle[i,j]
                max_i = i
                max_j = j

    return max_score, backtrack, max_i, max_j, middle

def output_lcs(backtrack, i, j, s1, s2, str1, str2, score):
    """
    Function to output alignment recursively
    INPUT: matrix with backtrack info, position i of max alignment,
    position j of max alignment, alignment s1, alignment s2, input s1, input s2,
    score matrix for when to stop local alignment
    OUTPUT: alignment of s1 and s2
    """
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
    """
    Write alignment and score to file
    INPUT: array of scores for pairs of sequences, array of alignments, and output file for writing
    OUTPUT: Tab separated file with score, alignment1, alignment2
    """
    out = open(outfile, "w")
    for i in range(0, len(scores)):
        out.write(str(int(scores[i]))+"\t"+"\t".join(alignments[i])+"\n")
    out.close()

def run_alignment(seqs, score_mat, gap=5, ext=3, normalize=False, 
                 out_align=False, outfile=None):
    """
    Run local alignment for a list of sequence pairs
    INPUT: list of sequences, scoring matrix, gap and extension penalty, whether to normalize
    the scores, whether to output the alignment and an output file if alignment set to true
    """
    scores = np.zeros(len(seqs)) # create empty array of scores for each sequence pair
    alignments = ['']*len(seqs) # create empty array of sequence pair alignments 
    # Loop through sequence pairs
    for i in range(0, len(seqs)): 
        # align sequences
        s, backtrack, max_i, max_j, middle = local_alignment(seqs[i][0], seqs[i][1], score_mat, gap, ext)
        # add score to matrix, normalize if set to True
        if normalize == False:
            scores[i] = s
        else:
            scores[i] = (s/min(len(seqs[i][0]), len(seqs[i][1])))
        # Output alignment if set to True
        if out_align == True:
            align_s = output_lcs(backtrack, max_i, max_j,seqs[i][0], seqs[i][1],"","", middle)
            alignments[i] = align_s
    # write alignments to file
    if out_align == True:
        write_align_output(scores, alignments, outfile)
    return scores

def output_align_scores_mats(pos_seqs, neg_seqs, matrix_dict, gap=5, ext=3):
    """
    Function that outputs alignment scores for all matrices
    INPUT: positive and negative sequence pairs, list of matrices, gap and ext penalty
    OUTPUT: dictionary where keys are matrices and value is array of scores for all sequence pairs
    """
    pos_scores, neg_scores = {}, {}
    for mat, scores in matrix_dict.items():
        pos_scores[mat] = run_alignment(pos_seqs, scores, gap, ext)
        neg_scores[mat] = run_alignment(neg_seqs, scores, gap, ext)
    return pos_scores, neg_scores

def find_best_align_penalties(pos,neg,gap_range, ext_range, score_mat):
    """
    Function to determine the best gap and extension penalities by finding the 
    smallest fpr with tpr of 0.7
    INPUT: positive and negative sequences, max gap and ext, dict of scores
    OUTPUT: file with optimal gap and extension 
    """
    bar = progressbar.ProgressBar(maxval=gap_range*ext_range, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    
    f = open("part1_1.txt", "w")
    best_fpr, best_gap, best_ext = float("inf"), None, None
    
    bar.start()
    i = 0
    for gap in range(1, gap_range+1):
        for ext in range(1, ext_range+1):
            # get alignment scores for a gap and extension penality
            scores_pos = np.array(run_alignment(pos, score_mat, gap, ext))
            scores_neg = np.array(run_alignment(neg, score_mat, gap, ext))
            # determine threshold for tpr of 0.7 in positive sequences
            thresh = calc_thresh(0.7, scores_pos)
            # get fpr using threshold
            fpr = len(scores_neg[scores_neg >= thresh])/len(scores_neg)
            # if this is the best fpr encountered, save it 
            if fpr < best_fpr:
                best_fpr, best_gap, best_ext = fpr, gap, ext
            bar.update(i+1)
            i+=1
            sleep(0.1)
    bar.finish()
    # write best gap and extension penalty to file
    f.write("fpr\tgap\text\n"+str(best_fpr)+"\t"+str(best_gap)+"\t"+str(best_ext))
    f.close()

def find_best_scoring_matrix(pos_scores_dict, neg_scores_dict, tpr=0.7):
    """
    Determine best scoring matrix by finding minimum fpr under tpr of 0.7
    INPUT: dictionary of alignment scores for each matrix
    OUTPUT: best performing matrix
    """
    best_fpr, best_mat = float("inf"), None
    mat_fprs = []
    # Loop through matrix and their alignment scores
    for mat, scores in pos_scores_dict.items():
        # determine threshold using fpr and positive sequences
        thresh = calc_thresh(tpr, scores)
        scores_neg = neg_scores_dict.get(mat)
        # calculate fpr
        fpr = len(scores_neg[scores_neg >= thresh])/len(scores_neg)
        mat_fprs.append(str(fpr))
        # track best fpr
        if fpr < best_fpr:
            best_fpr, best_mat = fpr, mat
    # write matrix performance in fpr to file
    out = open("part1_2a.txt", "w")
    out.write("\t".join(pos_scores_dict.keys())+"\n")
    out.write("\t".join(mat_fprs))
    out.close()
    return best_mat
