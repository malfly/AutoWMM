library(testthat)
library(AutoWMM)

# Only create expensive objects if not on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  tree <- create_tree()
  Zhats <- create_Zs()
} else {
  tree <- NULL
  Zhats <- NULL
}

# Test creation of tree data type
test_that("makeTree function returns correct s3 class", {
  skip_on_cran()
  expect_s3_class(tree, "Node")
})

# Test whether the output of confInts contains the right number of values
test_that("confInts() returns correct number of values", {
  skip_on_cran()
  expect_equal(length(confInts(Zhats$estimates)), 2)
})

# Test mmEstimate returns NA for non-leaf
test_that("mmEstimate() returns correct value for non-leaf", {
  skip_on_cran()
  mmEstimate(tree$A)
  expect_equal(tree$A$targetEst, NA)
})

# Test mmEstimate returns value for leaf at end of informative path
test_that("mmEstimate() returns correct value for leaf at end of informative path", {
  skip_on_cran()
  mmEstimate(tree$B)
  expect_type(tree$B$targetEst, "double")
})

# Test sampleBeta returns a number < 1
test_that("sampleBeta() returns double less than 1", {
  skip_on_cran()
  expect_type(sampleBeta(10, 55, pop = FALSE, tree$A), "double")
  expect_lt(sampleBeta(10, 55, pop = FALSE, tree$A), 1)
})

# Test ssEstimate returns NA for non-leaf
test_that("ssEstimate() returns correct value for non-leaf", {
  skip_on_cran()
  ssEstimate(tree$A)
  expect_equal(tree$A$targetEst, NA)
  expect_equal(tree$A$variance, NA)
})

# Test ssEstimate returns valid double values for informative leaf
test_that("ssEstimate() returns correct value for a valid leaf endpoint", {
  skip_on_cran()
  ssEstimate(tree$B)
  expect_type(tree$B$targetEst, "double")
  expect_type(tree$B$variance, "double")
})

# Test countTree returns list
test_that("countTree function returns correct type", {
  skip_on_cran()
  Zhats <- wmmTree(tree, sample_length = 2)
  expect_type(countTree(tree), "list")
})

# Test drawTree returns list
test_that("drawTree function returns correct type", {
  skip_on_cran()
  expect_type(drawTree(tree), "list")
})

# Test estTree returns list
test_that("estTree function returns correct type", {
  skip_on_cran()
  Zhats <- wmmTree(tree, sample_length = 2)
  expect_type(estTree(tree), "list")
})

# ---------------------
# Optional Heavy Tests — Keep Commented for CRAN
# ---------------------

# Test ko.weights returns correct number of values
# test_that("ko.weights() returns correct number of values", {
#   skip_on_cran()
#   tree <- makeTree(treeData1)
#   Zhats <- wmmTree(tree, sample_length = 2)
#   getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf, traversal = 'post-order'))
#   expect_equal(length(ko.weights(tree)), length(getleaves))
# })

# Test logEstimates returns correct number of values
# test_that("logEstimates() returns correct number of values", {
#   skip_on_cran()
#   tree <- makeTree(treeData1)
#   Zhats <- wmmTree(tree, sample_length = 2)
#   getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf, traversal = 'post-order'))
#   expect_equal(length(logEstimates(tree)), length(getleaves))
# })

# Test Nhats returns data frame with expected number of rows
# test_that("Nhats() returns correct number of values", {
#   skip_on_cran()
#   tree <- makeTree(treeData1)
#   Zhats <- wmmTree(tree, sample_length = 2)
#   expect_equal(nrow(Nhats(tree)), 2)
# })

# Test root.confInt returns length-2 vector
# test_that("root.confInt() returns correct number of values", {
#   skip_on_cran()
#   tree <- makeTree(treeData1)
#   Zhats <- wmmTree(tree, sample_length = 2)
#   expect_equal(length(root.confInt(tree)), 2)
# })

# Test ss.confInts returns length-2 vector
# test_that("ss.confInts() returns correct number of values", {
#   skip_on_cran()
#   tree <- makeTree(treeData1)
#   Zhats <- wmmTree(tree, sample_length = 2)
#   expect_equal(length(ss.confInts(tree$B)), 2)
# })

# Test wmmTree output has 4 components
# test_that("wmmTree function generates list with four entries", {
#   skip_on_cran()
#   tree <- makeTree(treeData1)
#   Zhats <- wmmTree(tree, sample_length = 2)
#   expect_equal(length(Zhats), 4)
# })
