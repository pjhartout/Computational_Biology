test.get_subtree_likelihood = function() {
  load("test_matrices_prov.RData")
  tree = read.tree(text ="((a:2,b:2):3,(c:2,d:8):1);")
  sequences = list(a = "TTAAG", b = "TCATG", c = "CCTAG", d = "TATGC")
  
  lik = get_subtree_likelihood(1, tree, sequences, TN93_Qfull)
  checkEqualsNumeric(lik_tip, lik, "Wrong likelihood at tip")
  
  lik = get_subtree_likelihood(6, tree, sequences, TN93_Qfull)
  checkEqualsNumeric(lik_node, lik, "Wrong likelihood at internal node")
}