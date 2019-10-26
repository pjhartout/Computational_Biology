test.get_best_alignment_local = function() {
  load("test_matrices_prov.RData")
  
  alignment = R.utils::withTimeout(get_best_alignment("TGCTCTGT", "TACTATCAT", score_matrix_l, path_matrix_l, T), timeout = 30, onTimeout = "error")
  checkEquals(10, alignment$score, "Wrong score for the local alignment")
  checkTrue(alignment$alignment[1] %in% c("TGCTCT", "TGCTCTG-T", "TGCTCT-GT"), "Wrong sequence 1 in the local alignment")
  checkTrue(alignment$alignment[2] %in% c("TACTAT", "TACTATCAT"), "Wrong sequence 2 in the local alignment")
  
  alignment = R.utils::withTimeout(get_best_alignment("TAGAC","ATGACT", score_matrix_l_2, path_matrix_l_2, T), timeout = 30, onTimeout = "error")
  checkEquals(10, alignment$score, "Wrong score for the local alignment")
  checkTrue(alignment$alignment[1] %in% c("A-GAC", "TAGAC"), "Wrong sequence 2 in the local alignment")
  checkTrue(alignment$alignment[2] %in% c("ATGAC", "T-GAC"), "Wrong sequence 1 in the local alignment")
}
