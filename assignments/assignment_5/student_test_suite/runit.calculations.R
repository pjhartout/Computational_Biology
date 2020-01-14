test.calculate_trait_at_internal_node = function() {
  t = calculate_trait_at_internal_node(3,2,1.5,1.5)
  checkEqualsNumeric(2.5, t, "Wrong trait value with equal branch lengths")
  
  t = calculate_trait_at_internal_node(5,2,1.5,3.5)
  checkEqualsNumeric(4.1, t, "Wrong trait value with unequal branch lengths")
}

test.calculate_contrast = function() {
  c = calculate_contrast(7,5,3,3)
  checkEqualsNumeric(0.8164966, c, "Wrong contrast value with equal branch lengths", tolerance = 1e-6)
  
  c = calculate_contrast(8,5,4,2)
  checkEqualsNumeric(1.224745, c, "Wrong contrast value with unequal branch lengths", tolerance = 1e-6)
}

test.calculate_corrected_branch_length = function() {
  b = calculate_corrected_branch_length(1,2,3,c(5,3,3))
  checkEqualsNumeric(6.5, b, "Wrong branch length value with unequal branch lengths")
  
  b = calculate_corrected_branch_length(1,2,3,c(5,9,3))
  checkEqualsNumeric(7.25, b, "Wrong branch length value with unequal branch lengths")
}

test.linear_regression = function() {
  values = linear_regression(seq(3,2,-0.01), c(rep(8,25), rep(4,25), rep(2,51)))
  checkEqualsNumeric(c(7.396622,-14.51136,4.287764e-34,0.7776760), values, "Wrong outcome of linear regression", tolerance = 1e-6)
}