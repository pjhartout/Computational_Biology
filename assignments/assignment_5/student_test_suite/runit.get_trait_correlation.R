test.get_trait_correlation = function() {
  tree = "((a:5,b:5):3,(c:2,d:8):6);"
  traits = list(trait1 = c(a = 2, b = 3, c = 8, d = 12),
                trait2 = c(a = 3, b = 4, c = 9, d = 10))
  result = get_trait_correlation(tree, traits)
  
  target = list(linear_regression_raw_data = c(slope = 0.73359073, intercept = 1.91505792, p.value = 0.02955068, R2 = 0.94177189),
                linear_regression_contrasts = c(slope = 0.7579596, intercept = 0.1034854, p.value = 0.4539161, R2 = 0.5721358))
  checkEquals(target$linear_regression_raw_data, result$linear_regression_raw_data, tolerance = 1e-6, "Wrong results of linear regression on raw data")
  checkEquals(target$linear_regression_contrasts, result$linear_regression_contrasts, tolerance = 1e-6, "Wrong results of linear regression on contrasts")
}
