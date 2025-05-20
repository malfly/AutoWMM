# This file is part of the standard setup for testthat.

Sys.setenv("NOT_CRAN" = "true")  # Enable full tests locally
library(testthat)
library(AutoWMM)

test_check("AutoWMM")
