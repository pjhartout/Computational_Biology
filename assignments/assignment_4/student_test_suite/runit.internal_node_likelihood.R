test.calculate_likelihood_from_subtree_likelihoods = function() {
  load("test_matrices_prov.RData")
  
  lik1 = matrix(c(0,0,1,0),nrow = 1)
  lik2 = matrix(c(1,0,0,0),nrow = 1)
  lik = calculate_likelihood_from_subtree_likelihoods(1, TN93_Qfull, lik1, 0.5, lik2, 0.5)
  checkEqualsNumeric(matrix(c(0.05223156,0.05223124,0.04938197,0.0492957), nrow = 1), lik, tolerance = 1e-4, 
                     "Wrong internal node likelihood")
  
  lik = calculate_likelihood_from_subtree_likelihoods(1, TN93_Qfull, lik1, 0.2, lik2, 0.4)
  checkEqualsNumeric(matrix(c(0.0258762,0.02587411,0.04960063,0.04459762), nrow = 1), lik, tolerance = 1e-4, 
                     "Wrong internal node likelihood")
  
  lik1 = matrix(c(0.015,0.045,0.022,0.04), nrow = 1)
  lik2 = matrix(c(0.017,0.025,0.042,0.031), nrow = 1)
  lik = calculate_likelihood_from_subtree_likelihoods(1, TN93_Qfull, lik1, 0.9, lik2, 1.4)
  checkEqualsNumeric(matrix(c(0.0008607208,0.0008607208,0.000926254,0.0009262558), nrow = 1), lik, tolerance = 1e-4, 
                     "Wrong internal node likelihood")
}