#' @title logEstimates
#' @description Estimation helper function:
#'  Creates a vector of mean estimate values given by each informative path
#' @param tree A makeTree object
#' @return Returns a vector of mean log estimate values of the root population size
#'  from each informative path
#' @examples \donttest{
#'   message("note - longer run time example")
#'   data(treeData1)
#'   tree <- makeTree(treeData1)
#'   Zhats <- wmmTree(tree, sample_length = 3)
#'   logEstimates(tree)
#' }
#' @export
#' @import data.tree
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom MASS "ginv"
#' @importFrom dplyr "select"

logEstimates <- function(tree){
  x <- tree$Get('targetEst_samples', filterFun = function(node) node$isLeaf,
                traversal = 'post-order')

  ## incase the above is not in matrix form...
  if(is.list(x)){
    mat.x <- NULL
    for (i in 1:length(x)) {
      if (length(x[[i]]) > 0) {
        mat.x <- cbind(mat.x, x[[i]])
        colnames(mat.x)[dim(mat.x)[2]] <- names(x)[[i]]
      }
    }
    x <- mat.x
  }

  # set as data frame
  x <- as.data.frame(log(x))

  ## select only leaves with marginal counts
  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf,
                              traversal = 'post-order'))
  # if number of columns of x is greater than number of leaves, choose leaves only
  if(dim(x)[2]>length(getleaves)){
    x <- x %>%
      select(all_of(getleaves))
  }

  weights <- as.matrix(ko.weights(tree)) #weights are internally calculated with log(data)
  logNhats <- as.matrix(x) %*% t(weights)

  return(logNhats)
}

