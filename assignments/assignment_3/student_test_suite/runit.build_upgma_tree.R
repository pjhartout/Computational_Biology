test.build_upgma_tree = function() {
  load("matrices_prov.RData")
  sequences = list(a = "GACCAGCGGTGACTAGAACGCAAGCAAATGTAGTCGAGCT",
                   b = "GTTCAGGGTTTATTAGCGCTCTATCCACTCCAGCCTAGCT",
                   c = "GATCAGGCCGTGTTAGCACCCACTCCGATCAAGCCTAGCT",
                   d = "GATCAGGACGGGTTAGCGCTCACTCCATTCCAGCCTAGCT")
  
  tree = R.utils::withTimeout(build_upgma_tree(sequences,"hamming"), timeout = 60, onTimeout = "error")
  checkEquals(ape::reorder.phylo(Htree), ape::reorder.phylo(tree), "Wrong UPGMA tree for Hamming distance")
  
  tree = R.utils::withTimeout(build_upgma_tree(sequences,"JC69"), timeout = 60, onTimeout = "error")
  checkEquals(ape::reorder.phylo(JCtree), ape::reorder.phylo(tree), "Wrong UPGMA tree for JC69 distance")
  
  tree = R.utils::withTimeout(build_upgma_tree(sequences,"K80"), timeout = 60, onTimeout = "error")
  checkEquals(ape::reorder.phylo(Ktree), ape::reorder.phylo(tree), "Wrong UPGMA tree for K80 distance")
}
