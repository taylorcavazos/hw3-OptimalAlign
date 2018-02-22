# File for optimizing scoring matrix
import numpy as np 
import pandas as pd
from opt_align import roc
from opt_align import align
import progressbar
from time import sleep


def mutate_matrix(mat):
	"""
	Choose random indices in score matrix to swap
	"""
	# create temporary matrix to mutate
	mut_mat = mat.copy()
	A,B = np.random.choice(range(len(mat.columns)), 2)
	index = list(mat.columns)
	index[A], index[B] = index[B], index[A]
	mut_mat.columns = index
	mut_mat.index = index
	return mut_mat

def calc_tpr(fpr, pos, neg):
	"""
	Calculate the TPR given a FPR
	"""
	thresh = int(roc.calc_thresh(fpr, neg))
	return len(pos[pos >= thresh])/len(pos)

def calc_obj(pos, neg):
	"""
	Calculate the objective function
	"""
	return calc_tpr(0, pos, neg) + calc_tpr(0.1, pos, neg) \
			 + calc_tpr(0.2, pos, neg) + calc_tpr(0.3, pos, neg)

def optimize(mat, pos_seqs, neg_seqs, pos_s0, neg_s0, max_idx=100):
	"""
	Optimize the scoring matrix using a genetic algorithm optimization
	"""
	bar = progressbar.ProgressBar(maxval=max_idx, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()

	obj_func = calc_obj(pos_s0, neg_s0)
	idx = 0
	while (obj_func < 4 and idx < max_idx):
		temp_mat = mutate_matrix(mat)
		pos_s = align.run_alignment(pos_seqs,temp_mat)
		neg_s = align.run_alignment(neg_seqs, temp_mat)
		obj_new = calc_obj(pos_s, neg_s)
		if obj_new > obj_func:
			obj_func = obj_new
			mat = temp_mat
		idx = idx+1
		bar.update(idx)
		sleep(0.1)
	bar.finish()
	return mat, obj_func


