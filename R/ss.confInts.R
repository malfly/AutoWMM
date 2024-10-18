#' @title ss.confInts
#' @description Estimation helper function:
#'  Method that takes samples and generates confidence intervals for
#'  nodes in single source sibling tree (single.source = TRUE)
#' @param o A node of a makeTree object
#' @param digits The number of significant digits to report
#' @return A confidence interval for nodes of the tree that only use a single
#'   source of sibling data.
#' @examples \donttest{
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  Zhats <- wmmTree(tree, sample_length = 3)
#'  ss.confInts(tree$B)
#' }
#' @export
#' @import data.tree
#' @importFrom stats var

ss.confInts <- function(o, digits = 3){
  # get estimate, variance of root from node
  m <- o$targetEst
  o$variance <- var(o$targetEst_samples)
  v <- o$variance

  # generate endpoints
  lc <- signif(m-1.96*sqrt(v), digits = digits)
  uc <- signif(m+1.96*sqrt(v), digits = digits)
  int <- c(lc, uc)
  names(int) <- c("lower", "upper")

  return(int)
}
