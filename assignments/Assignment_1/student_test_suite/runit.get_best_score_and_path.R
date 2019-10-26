test.get_best_score_and_path_global = function() {
  score_gap = -2
  score_match = 3
  score_mismatch = -1
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(0,1,1,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(3,score$score, "Wrong score in case of match")
  checkTrue(grepl(paste0("^",score$path), "diagonal"), "Wrong path in case of diagonal path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(5,1,1,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(4,score$score, "Wrong score in case of mismatch")
  checkTrue(grepl(paste0("^",score$path), "diagonal"), "Wrong path in case of diagonal path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,-2,3,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(1,score$score, "Wrong score in case of vertical gap")
  checkTrue(grepl(paste0("^",score$path), "up"), "Wrong path in case of up path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,3,-2,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(1,score$score, "Wrong score in case of horizontal gap")
  checkTrue(grepl(paste0("^",score$path), "left"), "Wrong path in case of left path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,0,0,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(-1,score$score, "Wrong score in case of mismatch, negative score")
  checkTrue(grepl(paste0("^",score$path), "diagonal"), "Wrong path in case of diagonal path, negative score")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(-1,-5,1,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(-1,score$score, "Wrong score in case of vertical gap, negative score")
  checkTrue(grepl(paste0("^",score$path), "up"), "Wrong path in case of up path, negative score")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(-1,1,-5,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(-1,score$score, "Wrong score in case of horizontal gap, negative score")
  checkTrue(grepl(paste0("^",score$path), "left"), "Wrong path in case of left path, negative score")
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(-1,4,0,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(2,score$score, "Wrong score when diagonal and left paths are possible")
  checkTrue(any(grepl(paste0("^",score$path), c("left","diagonal"))),"Wrong path when diagonal and left paths are possible")
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(-5,4,4,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(2,score$score, "Wrong score when left and up paths are possible")
  checkTrue(any(grepl(paste0("^",score$path), c("left","up"))),"Wrong path when left and up paths are possible")
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(-1,0,4,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(2,score$score, "Wrong score when diagonal and up paths are possible")
  checkTrue(any(grepl(paste0("^",score$path), c("up","diagonal"))),"Wrong path when diagonal and up paths are possible")

  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,0,0,0),nrow = 2),score_gap,score_match,score_mismatch,F)
  checkEquals(-1,score$score, "Wrong score in case of mismatch")
  checkTrue(grepl(paste0("^",score$path), "diagonal"), "Wrong path in case of diagonal path")
}
 
test.get_best_score_and_path_local = function() {
  score_gap = -2
  score_match = 3
  score_mismatch = -1
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(0,1,1,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(3,score$score, "Wrong score in case of match")
  checkTrue(grepl(paste0("^",score$path), "diagonal"), "Wrong path in case of diagonal path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(5,1,1,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(4,score$score, "Wrong score in case of mismatch")
  checkTrue(grepl(paste0("^",score$path), "diagonal"), "Wrong path in case of diagonal path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,-2,3,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(1,score$score, "Wrong score in case of vertical gap")
  checkTrue(grepl(paste0("^",score$path), "up"), "Wrong path in case of up path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,3,-2,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(1,score$score, "Wrong score in case of horizontal gap")
  checkTrue(grepl(paste0("^",score$path), "left"), "Wrong path in case of left path")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,0,0,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(0,score$score, "Wrong score in case of mismatch, negative score")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(-1,-5,1,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(0,score$score, "Wrong score in case of vertical gap, negative score")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(-1,1,-5,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(0,score$score, "Wrong score in case of horizontal gap, negative score")
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(-1,4,0,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(2,score$score, "Wrong score when diagonal and left paths are possible")
  checkTrue(any(grepl(paste0("^",score$path), c("left","diagonal"))),"Wrong path when diagonal and left paths are possible")
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(-5,4,4,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(2,score$score, "Wrong score when left and up paths are possible")
  checkTrue(any(grepl(paste0("^",score$path), c("left","up"))),"Wrong path when left and up paths are possible")
  
  score = get_best_score_and_path(2,2,"C","C",matrix(c(-1,0,4,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(2,score$score, "Wrong score when diagonal and up paths are possible")
  checkTrue(any(grepl(paste0("^",score$path), c("up","diagonal"))),"Wrong path when diagonal and up paths are possible")
  
  score = get_best_score_and_path(2,2,"C","G",matrix(c(0,0,0,0),nrow = 2),score_gap,score_match,score_mismatch,T)
  checkEquals(0,score$score, "Wrong score in case of mismatch")
}