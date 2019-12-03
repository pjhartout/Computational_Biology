tests_dir <- "/home/pjhartout/Polybox/ETHZ/Semester_1/Computational_Biology/assignments/assignment_5/student_test_suite/"

library("RUnit")

testsuite <- defineTestSuite("HW", tests_dir)
currentdir <- getwd()
setwd(tests_dir)

out <- runTestSuite(testsuite)
printTextProtocol(out)

setwd(currentdir)