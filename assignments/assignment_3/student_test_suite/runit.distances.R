test.get_hamming_distance = function() {
  d = get_hamming_distance("GCCAAT","GCCAAT")
  checkEqualsNumeric(0,d,"Wrong Hamming distance for equal sequences")
  
  d = get_hamming_distance("GCCAAT","GCCTTT")
  checkEqualsNumeric(2,d,"Wrong Hamming distance for different sequences")
}

test.get_JC69_distance = function() {
  d = get_JC69_distance("GCCAAT","GCCAAT")
  checkEqualsNumeric(0,d,"Wrong JC69 distance for equal sequences")
  
  d = get_JC69_distance("GCCAATT","GCCTTTA")
  checkEqualsNumeric(0.6354,d,"Wrong JC69 distance for different sequences", tolerance = 1e-3)
}

test.get_K80_distance = function() {
  d = get_K80_distance("GCCAAT","GCCAAT")
  checkEqualsNumeric(0,d,"Wrong K80 distance for equal sequences")
  
  d = get_K80_distance("GCCAATGC","GCCTATCC")
  checkEqualsNumeric(0.317,d,"Wrong K80 distance for different sequences with transversions", tolerance = 1e-3)
  
  d = get_K80_distance("GCCAAATTC","ACCGGATCC")
  checkEqualsNumeric(1.099,d,"Wrong K80 distance for different sequences with transitions", tolerance = 1e-3)
  
  d = get_K80_distance("ATTCAATTGCC","AAAGGATCCCA")
  checkEqualsNumeric(1.452,d,"Wrong K80 distance for different sequences with transitions and transversions", tolerance = 1e-3)
}