test.get_merge_node_distance = function() {
  load("matrices_prov.RData")
  df = data.frame(node_sizes = c(1,1,3,2) ,row.names = c("a","b","c","d"))
  
  md = get_merge_node_distance(df,Hdm,c("a","b"),"d")
  checkEqualsNumeric(13.5, md, "Wrong merge distance when merging tips")
  
  md = get_merge_node_distance(df,Hdm,c("a","c"),"d")
  checkEqualsNumeric(9.75, md, "Wrong merge distance when merging a tip and an internal node")
  
  md = get_merge_node_distance(df,Hdm,c("d","c"),"a")
  checkEqualsNumeric(18, md, "Wrong merge distance when merging internal nodes")
}