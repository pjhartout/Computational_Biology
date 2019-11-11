test.get_likelihood_from_sequence = function() {
  lik = get_likelihood_from_sequence(22, rep(2, 22))
  checkEquals(c(22,4), dim(lik), "Wrong dimensions for likelihood matrix")
  
  lik = get_likelihood_from_sequence(8, rep(4, 8))
  checkEqualsNumeric(rep(1,8), lik[,4], "Wrong likelihood for sequence GGGGGGGG")
  
  lik = get_likelihood_from_sequence(5, c(3,4,2,1,3))
  for(i in 1:20) {
    if(i %in% c(4,8,11,15,17)) {
      checkEqualsNumeric(1, lik[i], "Wrong likelihood for sequence GCATC")
    } else {
      checkEqualsNumeric(0, lik[i], "Wrong likelihood for sequence GCATC")
    }
  }
}