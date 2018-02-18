# Run code to complete assignment questions
from .align import *
from .io import *
from .roc import *
import sys

score_dir = sys.argv[1]
pos_file = sys.argv[2]
neg_file = sys.argv[3]

#################################### PART 1 ##########################################

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
print("Plot ROC coming performance of all scoring matrices")
performance = output_ROC_all(pos_scores, neg_scores)

# Plot ROC for best method with normalization
print("Plot ROC for best scoring matrix normalized vs not")
output_ROC_best(pos_seqs, neg_seqs, performance, score_mats, best_mat)

# Create alignments for each positive and negative pair of alignments
run_alignment(pos_seqs, score_mats.get(best_mat), out_align=True, outfile="alignments_pos.txt")
run_alignment(neg_seqs, score_mats.get(best_mat), out_align=True, outfile="alignments_neg.txt")
