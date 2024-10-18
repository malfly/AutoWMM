#' @title sampleBeta
#' @description Helper function:
#'  Method for sampling from a Beta distribution given the survey estimates
#' @param x An integer; typical use case is survey numerator
#' @param n An integer; typical use case is survey sample size
#' @param pop A logical value which takes TRUE if sample size is population size
#' @param node A node from a makeTree object; carried forever from
#' @return Returns sample from Beta distribution with parameters dependent on x, n
#' @examples \donttest{
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  Zhats <- wmmTree(tree, sample_length = 3)
#'  sampleBeta(10, 55, pop = FALSE, tree$A)
#' }
#' @export
#' @importFrom stats "rbeta"

sampleBeta <- function(x,n,pop,node){
  result <- NA
  # if x,n are population values and not sample estimates, set probability exactly
  if(pop){
    result <- x/n
  }
  # else, sample a probability from a beta distribution with parameters below
  else{
    result <- rbeta(1,x + 1,n - x + 1)
  }
  return(result)
}
