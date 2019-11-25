test.calculate_trait_and_contrast = function() {
  values = calculate_trait_and_contrast(1,2,c(4,12),c(4,4,3))
  checkEqualsNumeric(8, values$trait_value, "Wrong trait value with equal branch lengths")
  checkEqualsNumeric(-2.828427, values$contrast, "Wrong contrast value with equal branch lengths", tolerance = 1e-6)
  
  values = calculate_trait_and_contrast(1,2,c(8,5),c(2,5,6))
  checkEqualsNumeric(7.142857, values$trait_value, "Wrong trait value with unequal branch lengths", tolerance = 1e-6)
  checkEqualsNumeric(1.133893, values$contrast, "Wrong contrast value with unequal branch lengths", tolerance = 1e-6)
}