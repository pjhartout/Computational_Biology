tests_dir <- "/path/to/tests"
if(tests_dir == "path/to/tests") stop("tests_dir needs to be set to a proper path")

library("RUnit")

testsuite <- defineTestSuite("HW", tests_dir)
currentdir <- getwd()
setwd(tests_dir)

out <- runTestSuite(testsuite)
printTextProtocol(out)

setwd(currentdir)
