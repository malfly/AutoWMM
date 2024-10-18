#' @title confInts
#' @description Estimation helper function:
#'  Method that takes samples and generates confidence intervals for
#'  nodes other than the root. Assume raw data (not log), with normal distributed
#'  log data for confidence interval construction
#' @param v A vector
#' @return Returns a confidence interval. For non-root nodes
#' @examples \donttest{
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  message("note - longer run time example")
#'  Zhats <- wmmTree(tree, sample_length = 10)
#'  confInts(Zhats$estimates)
#' }
#' @export
#' @importFrom stats quantile

confInts <- function(v){
  v <- log(v)
  uc <- quantile(v,0.975,na.rm = TRUE)
  lc <- quantile(v,0.025,na.rm = TRUE)
  int <- c(exp(lc), exp(uc))
  names(int) <- c("lower", "upper")

  return(int)
}
