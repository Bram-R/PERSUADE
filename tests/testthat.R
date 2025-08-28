# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# This file tells R to run your tests when someone calls devtools::test()
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(PERSUADE)

test_check("PERSUADE")
