test.get_evolved_sequence = function() {
  load("test_matrices_prov.RData")

  cat("\n(This test can take up to a minute or so. Please be patient!)\n")

  testseq = get_evolved_sequence(rep(1, 100000),branch_length = 6, Q = TN93_Qfull)
  checkEquals(100000, length(testseq), "Wrong length of evolved sequence")
  
  checkEquals(0.22, sum(testseq == 1)/100000, msg = "Wrong proportion of nucleotide 1", tolerance = 4e-2)
  checkEquals(0.26, sum(testseq == 2)/100000, msg = "Wrong proportion of nucleotide 2", tolerance = 4e-2)
  checkEquals(0.33, sum(testseq == 3)/100000, msg = "Wrong proportion of nucleotide 3", tolerance = 4e-2)
  checkEquals(0.19, sum(testseq == 4)/100000, msg = "Wrong proportion of nucleotide 4", tolerance = 4e-2)
}