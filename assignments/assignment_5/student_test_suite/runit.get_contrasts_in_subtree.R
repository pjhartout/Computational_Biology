test.get_contrasts_in_subtree = function() {
  load("test_data_prov.RData")
  
  ftree = read.tree(text = "((a:5,b:5):3,(c:2,d:8):6);")
  
  lst2 = get_contrasts_in_subtree(7,ftree, lst)
  checkEqualsNumeric(7.6, lst2$corrected_branch_lengths[7], "Wrong corrected branch length")
  checkEqualsNumeric(c(8.8,9.2), lst2$traits[7,], "Wrong traits values")
  checkEqualsNumeric(c(-1.264911,-0.3162278), lst2$normalized_contrasts[7,], "Wrong contrast values", tolerance = 1e-6)
  
  lst2 = get_contrasts_in_subtree(5,ftree, lst)
  checkEqualsNumeric(5.5, lst2$corrected_branch_lengths[6], "Wrong corrected branch length")
  checkEqualsNumeric(c(5.725191,5.89313), lst2$traits[5,], "Wrong traits values", tolerance = 1e-6)
  checkEqualsNumeric(c(-1.4643343,-1.5748500), lst2$normalized_contrasts[5,], "Wrong contrast values", tolerance = 1e-6)
}
