test.get_starting_sequence = function() {
  pi = c(0.22, 0.26, 0.33, 0.19)
  seq = get_starting_sequence(pi, 100000)
  
  checkEquals(100000, length(seq), "Wrong length for starting sequence")
  
  checkEquals(pi[1], sum(seq == 1)/100000, msg = "Wrong proportion of nucleotide 1", tolerance = 2e-2)
  checkEquals(pi[2], sum(seq == 2)/100000, msg = "Wrong proportion of nucleotide 2", tolerance = 2e-2)
  checkEquals(pi[3], sum(seq == 3)/100000, msg = "Wrong proportion of nucleotide 3", tolerance = 2e-2)
  checkEquals(pi[4], sum(seq == 4)/100000, msg = "Wrong proportion of nucleotide 4", tolerance = 2e-2)
}
