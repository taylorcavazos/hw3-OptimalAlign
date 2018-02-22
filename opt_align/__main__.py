# Run code to complete assignment questions
from .align import *
from .io import *
from .roc import *
import sys
from .optimize import *

#################################### PART 1 ##########################################
if sys.argv[1] == "--align":
	score_dir = sys.argv[2]
	pos_file = sys.argv[3]
	neg_file = sys.argv[4]
	print("Loading scoring matrices")
	score_mats = read_mult_scores(score_dir)
	print("Reading sequences")
	pos_seqs = read_seqs(pos_file)
	neg_seqs = read_seqs(neg_file)

	# Find best gap and extension pentalities using BLOSUM50
	print("Finding best gap and extension penalties using BLOSUM50")
	find_best_align_penalties(pos_seqs, neg_seqs, 20, 5, score_mats.get("BLOSUM50"))

	# Output dictionary of scores for positive and negative sequences given 
	print("Calculate alignment scores for all matrices")
	pos_scores, neg_scores = output_align_scores_mats(pos_seqs, neg_seqs, score_mats)

	# Determine matrix that performs best with TPR=0.7
	print("Find best scoring matrix")
	best_mat = find_best_scoring_matrix(pos_scores, neg_scores)

	# Plot ROC for methods
	print("Plot ROC performance of all scoring matrices")
	performance = output_ROC_all(pos_scores, neg_scores)

	# Plot ROC for best method with normalization
	print("Plot ROC for best scoring matrix normalized vs not")
	output_ROC_best(pos_seqs, neg_seqs, performance, score_mats, best_mat)

	# Create alignments for each positive and negative pair of alignments
	run_alignment(pos_seqs, score_mats.get(best_mat), out_align=True, outfile="alignments_pos.txt")
	run_alignment(neg_seqs, score_mats.get(best_mat), out_align=True, outfile="alignments_neg.txt")

#################################### PART 2 ##############################################
elif sys.argv[1] == "--optimize":
	score_file = sys.argv[2]
	pos_seq_file = sys.argv[3]
	neg_seq_file = sys.argv[4]
	pos_align = sys.argv[5]
	neg_align = sys.argv[6]

	print("Loading scoring matrix")
	score_mat = read_score(score_file)

	print("Reading sequences")
	pos_seqs = read_seqs(pos_seq_file)
	neg_seqs = read_seqs(neg_seq_file)

	print("Reading initial pre-computed alginments")
	scores_pos_initial = read_alignment(pos_align)
	scores_neg_initial = read_alignment(neg_align)

	print("Optimizing the score matrix")
	init_fitness = calc_obj(scores_pos_initial, scores_neg_initial)
	opt_mat, opt_fitness = optimize(score_mat, pos_seqs, neg_seqs, scores_pos_initial, scores_neg_initial)

	print("Calculate alignment scores for optimized and original matrix")
	mat_name = score_file.split("/")[1]
	name1 = mat_name+" (fitness = "+str(round(init_fitness,2))+")"
	name2 = mat_name+"_OPT (fitness = "+str(round(opt_fitness,2))+")"
	score_matrices = {name1: score_mat, name2: opt_mat}
	pos_scores_all, neg_scores_all = output_align_scores_mats(pos_seqs, neg_seqs, score_matrices)

	print("Plot ROC performance for optimized and origional matrix")
	output_ROC_all(pos_scores_all, neg_scores_all, outfile="ROC_OPT_"+mat_name+".pdf")
	opt_mat.to_csv(score_file+"_OPT", sep = "\t", index = False)





