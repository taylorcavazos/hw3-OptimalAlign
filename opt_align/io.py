# File with functions for processing of input files
import numpy as np, pandas as pd, re, os
import matplotlib.pyplot as plt

def read_score(filename):
    """
    Read in scoring matrix and return dataframe of pairwise scores
    between possible elements in expected sequence
    """
    mat_file = open(filename).read().splitlines()
    mat_file = [line for line in mat_file if line[0] != "#"]
    aas = [m.group(0) for m in re.finditer("([A-Z,*])",mat_file[0])]
    df = pd.DataFrame(columns=aas, index=range(len(aas)))
    for i in range(len(mat_file[1:])):
        line = mat_file[1:][i]
        df.loc[i,:] = [int(m.group(0)) for m in re.finditer("-?[0-9]+", line)]
    df.index = aas
    return df

def read_mult_scores(dir_scores):
    """
    Given a directory of scoring matrices, create score dataframes
    output dictionary of score matrices
    """
    score_dict = {}
    for mat in os.listdir(dir_scores):
        score_dict[mat] = read_score(dir_scores+mat)
    return score_dict

def read_fasta(filename):
    """
    Open and read fasta file and output sequence string
    """
    file_fa = open(filename).read().splitlines()
    return "".join(file_fa[1:])

def read_seqs(filename):
    """
    Given a list of sequence pair filepaths, extract sequences
    and return list
    """
    lines = open(filename).read().splitlines()
    seqs = []
    for i in range(0, len(lines)):
        # get filepaths
        seq1_fa, seq2_fa = lines[i].split(" ")[0], lines[i].split(" ")[1]
        # read and extract sequences
        seq1, seq2 = read_fasta(seq1_fa), read_fasta(seq2_fa)
        # replace no-code elements with * in string to match lowest
        # scoring symbol in score matrix
        seq1, seq2 = seq1.replace("x","*"), seq2.replace("x","*")
        seqs.append((seq1, seq2))
    return seqs

