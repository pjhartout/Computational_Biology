test.align_local = function() {
  score_gap = -2
  score_match = +3
  score_mismatch = -1
  
  #local, one solution possible
  alignment = R.utils::withTimeout(align("TCACACGT", "AGCACACT", score_gap, score_match, score_mismatch, T), timeout = 30, onTimeout = "error")
  checkEquals(16, alignment$score, "Wrong score for local alignment")
  checkEquals("CACACGT", alignment$alignment[1], "Wrong sequence 1 for local alignment")
  checkEquals("CACAC-T", alignment$alignment[2], "Wrong sequence 2 for local alignment")
  
  #local, three solutions possible
  alignment = R.utils::withTimeout(align("TCTGAGTA", "ACGTGCTA", score_gap, score_match, score_mismatch, T), timeout = 30, onTimeout = "error")
  checkEquals(10, alignment$score, "Wrong score for local alignment with multiple solutions")
  checkTrue(alignment$alignment[1] %in% c("C-TGAGTA","C-TGAGTA","CTGAG-TA"), "Wrong sequence 1 for local alignment with multiple solutions")
  checkTrue(alignment$alignment[2] %in% c("CGTGC-TA","CGTG-CTA","C-GTGCTA"), "Wrong sequence 2 for local alignment with multiple solutions")
  checkTrue(which(c("CGTGC-TA","CGTG-CTA","C-GTGCTA")==alignment$alignment[2]) %in% which(c("C-TGAGTA","C-TGAGTA","CTGAG-TA")==alignment$alignment[1]), 
            "Mismatched sequences 1 and 2 when multiple solutions are possible (local alignment)")
}

test.align_global = function() {
    score_gap = -2
    score_match = +3
    score_mismatch = -1

    #global, one solution possible
    alignment = R.utils::withTimeout(align("TCTGAGTA", "ACGAGCTA", score_gap, score_match, score_mismatch, F), timeout = 30, onTimeout = "error")
    checkEquals(13, alignment$score, "Wrong score for global alignment")
    checkEquals("TCTGAG-TA", alignment$alignment[1], "Wrong sequence 1 for global alignment")
    checkEquals("AC-GAGCTA", alignment$alignment[2], "Wrong sequence 2 for global alignment")

    #global, two solutions possible
    alignment = R.utils::withTimeout(align("TGAGAGTA", "ACGAGAGA", score_gap, score_match, score_mismatch, F), timeout = 30, onTimeout = "error")
    checkEquals(13, alignment$score, "Wrong score for global alignment with multiple solutions")
    checkTrue(alignment$alignment[1] %in% c("-TGAGAGTA","T-GAGAGTA"), "Wrong sequence 1 for global alignment with multiple solutions")
    checkTrue(alignment$alignment[2] == "ACGAGAG-A", "Wrong sequence 2 for global alignment with multiple solutions")
}