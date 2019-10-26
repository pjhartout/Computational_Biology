test.init_score_matrix = function() {
  score_gap = -2
  test_mat = init_score_matrix(nrow=4,ncol=12,F,score_gap)
  checkEquals(c(4,12), dim(test_mat), "Wrong dimensions in the score matrix")
  
  checkEquals(seq(0, score_gap * 11, score_gap), test_mat[1,], "Wrong initial values in the upper row for local=F")
  checkEquals(seq(0, score_gap * 3, score_gap), test_mat[,1], "Wrong initial values in the left column for local=F")
  
  test_mat = init_score_matrix(nrow=4,ncol=12,T,score_gap)
  checkEquals(rep(0, 11), test_mat[1,][-1], "Wrong initial values in the upper row for local=T")
  checkEquals(rep(0, 3), test_mat[,1][-1], "Wrong initial values in the left column for local=T")
}

test.init_path_matrix = function() {
  test_mat = init_path_matrix(nrow=4,ncol=12, F)
  checkEquals(c(4,12), dim(test_mat), "Wrong dimensions in the path matrix for local=F")
  
  checkEquals(rep("left", 11), test_mat[1,][-1], "Wrong initial values in the upper row of the path matrix for local=F")
  checkEquals(rep("up", 3), test_mat[,1][-1], "Wrong initial values in the left column of the path matrix for local=F")
  checkEquals(sum(test_mat[-1,-1] == "left") + sum(test_mat[-1,-1] == "up") + sum(test_mat[-1,-1] == "diag"), 0, "Wrong initial values (not empty string) outside of upper row and left column for local=F")
  
  test_mat = init_path_matrix(nrow=4,ncol=12, T)
  checkEquals(c(4,12), dim(test_mat), "Wrong dimensions in the path matrix for local=T")
  checkEquals(sum(test_mat == "left") + sum(test_mat == "up") + sum(test_mat == "diag"), 0, "Path matrix is not empty for local=T")
}