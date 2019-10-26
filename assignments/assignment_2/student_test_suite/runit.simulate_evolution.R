test.simulate_evolution = function() {

  newick_tree = "(orangutan:13,(gorilla:10.25,(human:5.5,chimp:5.5):4.75):2.75);"
  N = 100
  alpha1 = 450/8.406; alpha2 = 220/8.406; beta = 10/8.406
  pi = c(0.22, 0.26, 0.33, 0.19)
  result = simulate_evolution(newick_tree, pi, alpha1, alpha2, beta, N)
  
  checkEquals(4, length(result), "Wrong number of sequences in final result")
  for(i in 1:4) checkEquals(N, nchar(result[[i]]), "Wrong number of characters in final sequence")
  
  old_get_evolved_sequence <<- get_evolved_sequence
  gev_has_been_called <- FALSE
  get_evolved_sequence <<- function(sequence, branch_length, Q) {
      gev_has_been_called <<- TRUE
      stop()
  }
  try(simulate_evolution(newick_tree, pi, alpha1, alpha2, beta, N), silent=TRUE)
  checkEquals(TRUE, gev_has_been_called, "Function simulate_evolution does not use get_evolved_sequence")
  get_evolved_sequence <<- old_get_evolved_sequence
}