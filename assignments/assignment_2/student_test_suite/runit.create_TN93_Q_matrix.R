test.create_TN93_Q_matrix = function() {
  load("test_matrices_prov.RData")
  
  pi = c(0.22, 0.26, 0.33, 0.19)
  alpha1 = 450/8.406; alpha2 = 220/8.406; beta = 10/8.406
  
  matrix = create_TN93_Q_matrix(pi, alpha1, alpha2, beta)
  checkEquals(c(4,4), dim(matrix), "Wrong dimensions for TN93 matrix")
  checkEqualsNumeric(TN93_Qfull, matrix, "Wrong TN93 rate matrix")
}
