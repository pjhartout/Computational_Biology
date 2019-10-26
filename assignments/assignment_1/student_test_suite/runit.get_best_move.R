test.get_best_move = function() {
  path = get_best_move("A","T", "up", 3, 5)
  checkEquals(2, path$newrow, "Wrong new row index with up move")
  checkEquals(5, path$newcol, "Wrong new column index with up move")
  checkEquals("A", path$char1, "Wrong char1 with up move")
  checkEquals("-", path$char2, "Wrong char2 with up move")
  
  path = get_best_move("A","G", "diag",  3, 4)
  checkEquals(2, path$newrow, "Wrong new row index with diagonal move")
  checkEquals(3, path$newcol, "Wrong new column index with diagonal move")
  checkEquals("A", path$char1, "Wrong char1 with diagonal move")
  checkEquals("G", path$char2, "Wrong char2 with diagonal move")
  
  path = get_best_move("G","C","left", 3, 2)
  checkEquals(3, path$newrow, "Wrong new row index with left move")
  checkEquals(1, path$newcol, "Wrong new column index with left move")
  checkEquals("-", path$char1, "Wrong char1 with left move")
  checkEquals("C", path$char2, "Wrong char2 with left move")
}
