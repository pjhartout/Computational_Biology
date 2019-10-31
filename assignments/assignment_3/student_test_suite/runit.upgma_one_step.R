test.upgma_one_step = function() {
  load("matrices_prov.RData")
  sequences = list(a = "GACCAGCGGTGACTAGAACGCAAGCAAATGTAGTCGAGCT",
                   b = "GTTCAGGGTTTATTAGCGCTCTATCCACTCCAGCCTAGCT",
                   c = "GATCAGGCCGTGTTAGCACCCACTCCGATCAAGCCTAGCT",
                   d = "GATCAGGACGGGTTAGCGCTCACTCCATTCCAGCCTAGCT")
  
  df = data.frame(node_sizes = rep(1,4), node_heights = rep(0,4) ,row.names = c("a","b","c","d"))
  
  result = upgma_one_step(df,JCdm,matrix(nrow = 0, ncol = 2),c())
  checkEqualsNumeric(c(0.09963869,0.09963869), result$edge_lengths, "Wrong edge lengths, step 1", tolerance = 1e-5)
  checkEqualsNumeric(2, result$node_description$node_sizes[5], "Wrong size for internal node, step 1")
  checkEqualsNumeric(0.09963869, result$node_description$node_heights[5], "Wrong height for internal node, step 1", tolerance = 1e-5)
  
  result = upgma_one_step(result$node_description,result$distance_matrix,result$edges,result$edge_lengths)
  checkEqualsNumeric(c(0.06301767,0.16265636), sort(result$edge_lengths[3:4]), "Wrong edge lengths, step 2", tolerance = 1e-5)
  checkEqualsNumeric(3, result$node_description$node_sizes[6], "Wrong size for internal node, step 2")
  checkEqualsNumeric(0.16265636, result$node_description$node_heights[6], "Wrong height for internal node, step 2", tolerance = 1e-5)
}
