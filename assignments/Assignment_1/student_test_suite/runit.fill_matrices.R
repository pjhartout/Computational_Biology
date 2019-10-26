test.fill_matrices = function() {
    load("test_matrices_prov.RData")

    # local
    matrices = fill_matrices("TGCTCTGT","TACTATCAT", score_gap, score_match, score_mismatch, TRUE, score_matrix_l_init, path_matrix_l_init)
    checkEqualsNumeric(matrices$score_matrix, score_matrix_l, "Filled score matrix does not match for local alignment")
    for(i in 2:length(matrices$path_matrix)) {
      if(score_matrix_l[i] != 0) checkTrue(matrices$path_matrix[i] %in% lpm[i][[1]], "Filled path matrix does not match for local alignment")
    }
    # global
    matrices = fill_matrices("TGCTCTGT","TACTATCAT", score_gap, score_match, score_mismatch, FALSE, score_matrix_g_init, path_matrix_g_init)
    checkEqualsNumeric(matrices$score_matrix, score_matrix_g, "Filled score matrix does not match for global alignment")
    for(i in 2:length(matrices$path_matrix)) {
      checkTrue(matrices$path_matrix[i] %in% gpm[i][[1]], "Filled path matrix does not match for global alignment")
    }
}
