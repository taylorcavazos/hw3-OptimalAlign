from opt_align import align
from opt_align import io

def test_matrix_load():
    score_mats = io.read_mult_scores("score_matrices/")
    assert sorted(score_mats.keys()) == sorted(["BLOSUM50","BLOSUM62", "MATIO", "PAM100", "PAM250"])

def test_sequence_load_first():
    first = io.read_seqs("PosPairs.txt")[0]
    assert first[0] == "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
    assert first[1] == "ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEINK"

def test_sequence_load_last():
    sequences = io.read_seqs("PosPairs.txt")
    last = sequences[len(sequences)-1]
    assert last[0] == "MDSVCPQGKYIHPQNNSICCTKCHKGTYLYNDCPGPGQDTDCRECESGSFTASENHLRHC"
    assert last[1] == "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"

