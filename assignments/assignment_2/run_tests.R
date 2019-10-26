tests_dir <- "~/Dropbox/ETHZ/Semester_1/Computational_Biology/Assignments/Assignment_2/student_test_suite/"
if(tests_dir == "~/Dropbox/ETHZ/Semester_1/Computational_Biology/Assignments/Assignment_2/run_tests.R") stop("tests_dir needs to be set to a proper path")

library("RUnit")

testsuite <- defineTestSuite("Philip_Hartout.R", tests_dir)
currentdir <- getwd()
setwd(tests_dir)

out <- runTestSuite(testsuite)
printTextProtocol(out)

setwd(currentdir)
