test.compute_initial_distance_matrix = function() {
  load("matrices_prov.RData")
  sequences = list(a = "GACCAGCGGTGACTAGAACGCAAGCAAATGTAGTCGAGCT",
                   b = "GTTCAGGGTTTATTAGCGCTCTATCCACTCCAGCCTAGCT",
                   c = "GATCAGGCCGTGTTAGCACCCACTCCGATCAAGCCTAGCT",
                   d = "GATCAGGACGGGTTAGCGCTCACTCCATTCCAGCCTAGCT")
  
  matrix = compute_initial_distance_matrix(sequences, "hamming")
  for(i in 1:(length(sequences)-1)) {
    for(j in (i+1):length(sequences)) {
        namei = rownames(Hdm)[i]
        namej = colnames(Hdm)[j]
        checkEqualsNumeric(Hdm[namei,namej], matrix[namei,namej], "Wrong Hamming distance matrix")
    }
  }
  
  matrix = compute_initial_distance_matrix(sequences, "JC69")
  for(i in 1:(length(sequences)-1)) {
    for(j in (i+1):length(sequences)) {
        namei = rownames(JCdm)[i]
        namej = colnames(JCdm)[j]
        checkEqualsNumeric(JCdm[namei,namej], matrix[namei,namej], "Wrong JC69 distance matrix")
    }
  }
  
  matrix = compute_initial_distance_matrix(sequences, "K80")
  for(i in 1:(length(sequences)-1)) {
    for(j in (i+1):length(sequences)) {
        namei = rownames(Kdm)[i]
        namej = colnames(Kdm)[j]
        checkEqualsNumeric(Kdm[namei,namej], matrix[namei,namej], "Wrong K80 distance matrix")
    }
  }
}