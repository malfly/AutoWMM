library(testthat)
library(AutoWMM)

# create tree and Zhats
tree <- create_tree()
Zhats <- create_Zs()

# Test creation of tree data type
test_that("makeTree function returns correct s3 class", {
  expect_s3_class(tree, "Node")
})

# Test whether the output of confInts contains the right number of values
test_that("confInts() returns correct number of values", {
  expect_equal(length(confInts(Zhats$estimates)), 2)
})

# Test whether the output of ko.weights contains the right number of values
#test_that("ko.weights() returns correct number of values", {
# tree <- makeTree(treeData1)
#  Zhats <- wmmTree(tree, sample_length = 2)
#  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf,
#                              traversal = 'post-order'))
#  expect_equal(length(ko.weights(tree)), length(getleaves))
#})

# Test whether the output of logEstimates contains the right number of values
#test_that("logEstimates() returns correct number of values", {
#  tree <- makeTree(treeData1)
#  Zhats <- wmmTree(tree, sample_length = 2)
#  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf,
#                              traversal = 'post-order'))
#  expect_equal(length(logEstimates(tree)), length(getleaves))
#})

# Test whether the output of mmEstimate returns NA for non-leaf
test_that("mmEstimate() returns correct value for non-leaf", {
  mmEstimate(tree$A)
  expect_equal(tree$A$targetEst, NA)
})

# Test whether the output of mmEstimate returns NA for leaf at end of informative path
test_that("mmEstimate() returns correct value for leaf at end of informative path", {
  mmEstimate(tree$B)
  expect_type(tree$B$targetEst, "double")
})

# Test whether the output of Nhats contains the right number of values
#test_that("Nhats() returns correct number of values", {
#  tree <- makeTree(treeData1)
#  Zhats <- wmmTree(tree, sample_length = 2)
#  expect_equal(nrow(Nhats(tree)), 2)
#})

# Test whether the output of root.confInt contains the right number of values
#test_that("root.confInt() returns correct number of values", {
#  tree <- makeTree(treeData1)
#  Zhats <- wmmTree(tree, sample_length = 2)
#  expect_equal(length(root.confInt(tree)), 2)
#})

# Test whether the output of sampleBeta returns double less than 1
test_that("sampleBeta() returns double less than 1", {
  expect_type(sampleBeta(10, 55, pop = FALSE, tree$A), "double")
  expect_lt(sampleBeta(10, 55, pop = FALSE, tree$A), 1)
})

# Test whether the output of ss.confInts contains the right number of values
#test_that("ss.confInts() returns correct number of values", {
#  tree <- makeTree(treeData1)
#  Zhats <- wmmTree(tree, sample_length = 2)
#  expect_equal(length(ss.confInts(tree$B)), 2)
#})

# Test whether the output of ssEstimate returns NA for non-leaf
test_that("ssEstimate() returns correct value for non-leaf", {
  ssEstimate(tree$A)
  expect_equal(tree$A$targetEst, NA)
  expect_equal(tree$A$variance, NA)
})

# Test whether the output of ssEstimate returns NA for leaf endpoint of informative path
test_that("ssEstimate() returns correct value for a valid leaf endpoint", {
  ssEstimate(tree$B)
  expect_type(tree$B$targetEst, "double")
  expect_type(tree$B$variance, "double")
})

# Test creation of correct type with countTree
test_that("countTree function returns correct type", {
  skip_on_cran()
  Zhats <- wmmTree(tree, sample_length = 2)
  expect_type(countTree(tree), "list")
})

# Test creation of correct type with drawTree
test_that("drawTree function returns correct type", {
  expect_type(drawTree(tree), "list")
})

# Test creation of correct type with estTree
test_that("estTree function returns correct type", {
  skip_on_cran()
  Zhats <- wmmTree(tree, sample_length = 2)
  expect_type(estTree(tree), "list")
})

# Test wmmTree generates list with 4 entries
#test_that("wmmTree function generates list with four entries", {
  #tree <- makeTree(treeData1)
#  Zhats <- wmmTree(tree, sample_length = 2)
#  expect_equal(length(Zhats), 4)
#})
