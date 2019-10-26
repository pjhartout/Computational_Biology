test.get_starting_nucleotide = function() {
  pi = c(0.22, 0.26, 0.33, 0.19)
  samples = numeric(100000)
  for(i in 1:100000) samples[i] = get_starting_nucleotide(pi)
  
  checkEquals(pi[1], sum(samples == 1)/100000, msg = "Wrong proportion of nucleotide 1", tolerance = 2e-2)
  checkEquals(pi[2], sum(samples == 2)/100000, msg = "Wrong proportion of nucleotide 2", tolerance = 2e-2)
  checkEquals(pi[3], sum(samples == 3)/100000, msg = "Wrong proportion of nucleotide 3", tolerance = 2e-2)
  checkEquals(pi[4], sum(samples == 4)/100000, msg = "Wrong proportion of nucleotide 4", tolerance = 2e-2)
}
