test.update_distance_matrix = function() {
  load("matrices_prov.RData")
  df = data.frame(node_sizes = c(2,4,1,1) ,row.names = c("a","b","c","d"))
  
  matrix = update_distance_matrix(df, Hdm, c("c","d"), "c.d")
  for(i in 1:2) {
    for(j in (i+1):3) {
      namei = rownames(Hdm_one)[i]
      namej = colnames(Hdm_one)[j]
      checkEqualsNumeric(Hdm_one[namei,namej], matrix[namei,namej], "Wrong update when merging two tips")
    }
  }
  
  matrix = update_distance_matrix(df, Hdm, c("b","d"), "b.d")
  for(i in 1:2) {
    for(j in (i+1):3) {
      namei = rownames(Hdm_two)[i]
      namej = colnames(Hdm_two)[j]
      checkEqualsNumeric(Hdm_two[namei,namej], matrix[namei,namej], "Wrong update when merging a tip and internal node")
    }
  }
  
  matrix = update_distance_matrix(df, Hdm, c("a","b"), "a.b")
  for(i in 1:2) {
    for(j in (i+1):3) {
      namei = rownames(Hdm_three)[i]
      namej = colnames(Hdm_three)[j]
      checkEqualsNumeric(Hdm_three[namei,namej], matrix[namei,namej], "Wrong update when merging two internal nodes")
    }
  }
}