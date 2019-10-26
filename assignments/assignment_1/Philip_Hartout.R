#####################################################
#######        COMPUTATIONAL BIOLOGY         ########
#######             HOMEWORK 1               ########
#####################################################
#                                                   #
# Implement the pairwise alignment algorithms       #
# Needleman-Wunsch and Smith-Waterman.              #
#                                                   #
#####################################################
#####################################################

# In all functions the following parameters are the same:
# seqA: the first sequence to align
# seqB: the second sequence to align
# score_gap: score for a gap
# score_match: score for a character match
# score_mismatch: score for a character mismatch
# local: True if alignment is local, False otherwise

library(utils, lib.loc = "/usr/lib/R/library")

init_score_matrix = function(nrow, ncol, local, score_gap) {
    # Initialize the score matrix with zeros.
    # If the alignment is global, the leftmost column and the top row will have incremental gap scores,
    # i.e. if the gap score is -2 and the number of columns is 4, the top row will be [0, -2, -4, -6].
    # nrow: number of rows in the matrix
    # ncol: number of columns in the matrix
    
    # Initialize matrix.
    score_matrix = matrix(0, nrow, ncol)
    
    # Initialize gap score incrementally along columns and axes when the alignment is global.
    if (local == F){
      for (i in 2:nrow) {
        score_matrix[i, 1] =  score_matrix[i-1, 1] + score_gap
      }
      for (i in 2:ncol) {
        score_matrix[1, i] =  score_matrix[1, i-1] + score_gap
      }
    }
    
    # Return the initialized empty score matrix
    # score_matrix: nrow by ncol matrix
    return(score_matrix)
}

init_path_matrix = function(nrow, ncol, local) {
    # Initialize the path matrix with empty values, except the top row and the leftmost column if global alignment.
    # If global alignment, the top row has "left" on all positions except 1st.
    # Similarly, leftmost column has "up" on all positions except 1st.
    # nrow: number of rows in the matrix
    # ncol: number of columns in the matrix

    path_matrix <- matrix("", nrow, ncol)
    
    if (local == F) {
      path_matrix[1, ] = "left"
      path_matrix[, 1] = "up"
      path_matrix[1, 1] = ""
    }
    
    # Return the initialized empty path matrix
    # path_matrix: nrow by ncol matrix
    return(path_matrix)
}

get_best_score_and_path = function(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) {
    # Compute the score and the best path for a particular position in the score matrix
    # nucA: nucleotide in sequence A
    # nucB: nucleotide in sequence B
    # row: row-wise position in the matrix
    # col: column-wise position in the matrix
    # score_matrix: the score_matrix that is being filled out
    
    # Compute score for each possibility
    # Calculate diagonal score
    if (nucA == nucB){
      score_diag = score_matrix[row - 1, col - 1] + score_match
    } else {
      score_diag = score_matrix[row - 1, col - 1] + score_mismatch
    }
  
    # Calculate vertical and horizontal scores
    score_vertical = score_matrix[row - 1, col] + score_gap
    score_horizontal = score_matrix[row, col - 1] + score_gap
    
    # Possible scores
    if (local == T) {
      possible_scores = c(score_horizontal, score_diag, score_vertical, 0)  
    }
    else {
      possible_scores = c(score_horizontal, score_diag, score_vertical)
    }
    
    # Get maximum score 
    score = max(possible_scores)
    # Get the position fot the score to determine the path
    max_score_index = which(possible_scores==max(possible_scores))[1]
    
    if (max_score_index == 1) {
      path = "left"
    }
    else if (max_score_index == 2){
      path = "diag"    
    }
    else if (max_score_index == 3) {
      path = "up"
    }
    else if (max_score_index == 4 & local == T) {
      path = "-"
    }
    
    # Return the best score for the particular position in the score matrix
    # In the case that there are several equally good paths available, return any one of them.
    # score: best score at this position
    # path: path corresponding to the best score, one of ["diag", "up", "left"] in the global case and of ["diag", "up", "left", "-"] in the local case
    return(list("score"=score, "path"=path))
}

fill_matrices = function(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix) {
    # Compute the full score and path matrices
    # score_matrix: initial matrix of the scores
    # path_matrix: initial matrix of paths

    for (i in 2:dim(score_matrix)[1]) {
      for (j in 2:dim(score_matrix)[2]) {
        nucA = substr(seqA, i-1, i-1)
        nucB = substr(seqB, j-1, j-1)
        score_path = get_best_score_and_path(row=i, col=j, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local)
        score_matrix[i,j] = score_path[["score"]]
        path_matrix[i,j] = score_path[["path"]]  
      }
    }
  
    # Return the full score and path matrices
    # score_matrix: filled up matrix of the scores
    # path_matrix: filled up matrix of paths
    return(list("score_matrix"=score_matrix, "path_matrix"=path_matrix))
}

get_best_move = function(nucA, nucB, path, row, col) {
    # Compute the aligned characters at the given position in the score matrix and return the new position,
    # i.e. if the path is diagonal both the characters in seqA and seqB should be added,
    # if the path is up or left, there is a gap in one of the sequences.
    # nucA: nucleotide in sequence A
    # nucB: nucleotide in sequence B
    # path: best path pre-computed for the given position
    # row: row-wise position in the matrix
    # col: column-wise position in the matrix
  
    # print(path)
    if (path == "left") {
      # gap in seqA
      newrow = row
      newcol = col - 1
      char1 = '-'
      char2 = nucB
    }
    else if (path == "up") {
      # gap in seqB
      newrow = row - 1
      newcol = col
      char1 = nucA
      char2 = '-'
    }
    else {
      # match/mismatch
      newrow = row - 1
      newcol = col - 1
      char1 = nucA
      char2 = nucB
    }

    # Return the new row and column and the aligned characters
    # newrow: row if gap in seqA, row - 1 otherwise
    # newcol: col if gap in seqB, col - 1 otherwise
    # char1: '-' if gap in seqA, appropriate character if a match
    # char2: '-' if gap in seqB, appropriate character if a match
    return(list("newrow"=newrow, "newcol"=newcol, "char1"=char1, "char2"=char2))
}

get_best_alignment = function(seqA, seqB, score_matrix, path_matrix, local) {
    # Return the best alignment from the pre-computed score matrix
    # score_matrix: filled up matrix of the scores
    # path_matrix: filled up matrix of paths
    
    alignment = c("", "")
    
    if (local == T) {
      # Optimize
      start_position = which(score_matrix == max(score_matrix), arr.ind = TRUE)[1,]
      score = score_matrix[start_position[1], start_position[2]]
      i = start_position[1] 
      j = start_position[2]
      while (i > 1 & j > 1 & score_matrix[i, j] != 0) {
        nucA = substr(seqA, i-1, i-1)
        nucB = substr(seqB, j-1, j-1)
        path = path_matrix[i,j]
        row = i
        col = j
        best_move_info = get_best_move(nucA, nucB, path, row, col)
        i = best_move_info[["newrow"]]
        j = best_move_info[["newcol"]]
        alignment[1] = paste0(alignment[1], best_move_info["char1"])
        alignment[2] = paste0(alignment[2], best_move_info["char2"])
      }
    } 
    else {
      score = score_matrix[dim(score_matrix)[1],dim(score_matrix)[2]]
      i = dim(score_matrix)[1] 
      j = dim(score_matrix)[2]
      while (i >= 1 & j >= 1) {
        nucA = substr(seqA, i-1, i-1)
        nucB = substr(seqB, j-1, j-1)
        path = path_matrix[i,j]
        row = i
        col = j
        best_move_info = get_best_move(nucA, nucB, path, row, col)
        i = best_move_info[["newrow"]]
        j = best_move_info[["newcol"]]
        alignment[1] = paste0(alignment[1], best_move_info["char1"])
        alignment[2] = paste0(alignment[2], best_move_info["char2"])
      }
    }
    
    #Reverse string to get alignment
    for (i in 1:length(alignment)) {
      tmp_sequence = strsplit(alignment[i], NULL)[[1]]
      tmp_sequence = rev(tmp_sequence)
      alignment[i] = paste(tmp_sequence, collapse = '')
    }
  
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # score: score of the best alignment
    # alignment: the actual alignment in the form of a vector of two strings
    return(list("score"=score, "alignment"=alignment))
}

align = function(seqA, seqB, score_gap, score_match, score_mismatch, local) {
    # Align the two sequences given the scoring scheme
    # For testing purposes, use seqA for the rows and seqB for the columns of the matrices
  
    # Initialize score and path matrices
    score_matrix = init_score_matrix(nrow = nchar(seqA)+1, ncol = nchar(seqB)+1, local = local, score_gap = score_gap)
    path_matrix = init_path_matrix(nrow = nchar(seqA)+1, ncol = nchar(seqB)+1, local = local)
 
    # Fill in the matrices with scores and paths using dynamic programming
    matrices = fill_matrices(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix)
    
    score_matrix = matrices[["score_matrix"]]
    path_matrix = matrices[["path_matrix"]]
    
    # Get the best score and alignment (or one thereof if there are multiple with equal score)
    result = get_best_alignment(seqA, seqB, score_matrix, path_matrix, local)
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # Returns the same value types as get_best_alignment
    return(result)
}

test_align = function() {
    seqA = "TCTGAGTA" 
    seqB = "ACGAGCTA"
    score_gap = -2
    score_match = +3
    score_mismatch = -1
    local = F
    result = align(seqA, seqB, score_gap, score_match, score_mismatch, local)
    print(result$alignment)
    print(result$score)
}

test_align()