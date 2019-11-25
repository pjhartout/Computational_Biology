test.Felsensteins_pruning_loglikelihood = function() {
  tree = "((a:5,b:5):3,(c:2,d:8):1);"
  sequences = list(a = "TTAAG", b = "TCATG", c = "CCTAG", d = "TATGC")
  pi = c(0.22, 0.26, 0.33, 0.19)
  alpha1 = 53; alpha2 = 26; beta = 1.2
  
  logL = Felstensteins_pruning_loglikelihood(pi, alpha1, alpha2, beta, tree, sequences)
  checkEqualsNumeric(-28.173393875, logL, tolerance = 1e-7, "Wrong log likelihood for the full tree")
}