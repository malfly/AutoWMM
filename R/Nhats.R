#' @title Nhats
#' @description Estimation helper function:
#'  Returns Nhat for each sample from the WMM (rather than the aggregate
#'  value given by the average, this calculates weights and applies the
#'  weighted sum to each of the samples)
#' @param tree A makeTree object
#' @return Returns a vector of root population size estimates from each sample
#'  run of the wmmTree function
#' @examples \donttest{
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  Zhats <- wmmTree(tree, sample_length = 3)
#'  Nhats(tree)
#' }
#' @export
#' @importFrom dplyr "select"
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @import data.tree

Nhats <- function(tree){
  weights <- as.matrix(ko.weights(tree))

  # extract sample values for each path with terminal node counts and
  # informative paths
  x <- tree$Get('targetEst_samples', filterFun = function(node) node$isLeaf,
                traversal = 'post-order')

  ## incase the above is not in matrix form
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

  # select only leaves with marginal counts
  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf,
                              traversal = 'post-order'))

  # if number of columns of x is greater than number of leaves, choose leaves only
  if(dim(x)[2]>length(getleaves)){
    x <- x %>%
      select(all_of(getleaves))
  }

  # calculate logN for each sample
  logNhats <- as.data.frame(as.matrix(log(x)) %*% t(weights))

  # return unlogged values, Nhat
  return(round(exp(logNhats)))
}
