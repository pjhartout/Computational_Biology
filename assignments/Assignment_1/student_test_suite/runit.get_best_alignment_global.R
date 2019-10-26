test.get_best_alignment_global = function() {
  load("test_matrices_prov.RData")
  #global
  alignment = R.utils::withTimeout(get_best_alignment("TGCTCTGT", "TACTATCAT", score_matrix_g, path_matrix_g, F), timeout = 30, onTimeout = "error")
  checkEquals(10, alignment$score, "Wrong score for the global alignment")
  checkTrue(alignment$alignment[1] %in% c("TGCTCTG-T", "TGCTCT-GT"), "Wrong sequence 1 in the global alignment")
  checkEquals("TACTATCAT", alignment$alignment[2], "Wrong sequence 2 in the global alignment")
  
  alignment = R.utils::withTimeout(get_best_alignment("TAGAC", "ATGACT", score_matrix_g_2, path_matrix_g_2, F), timeout = 30, onTimeout = "error")
  checkEquals(6, alignment$score, "Wrong score for the global alignment")
  checkTrue(alignment$alignment[1] %in% c("-TAGAC-", "-TA-GAC-"), "Wrong sequence 1 in the global alignment")
  checkTrue(alignment$alignment[2] %in% c("AT-GACT", "-ATGACT"), "Wrong sequence 2 in the global alignment")
}
