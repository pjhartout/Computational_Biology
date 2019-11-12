tests_dir <- "C:/Users/phili/Dropbox/ETHZ/Semester_1/Computational_Biology/assignments/assignment_4/student_test_suite"
library("RUnit")

testsuite <- defineTestSuite("HW", tests_dir)
currentdir <- getwd()
setwd(tests_dir)

out <- runTestSuite(testsuite)
printTextProtocol(out)

setwd(currentdir)
