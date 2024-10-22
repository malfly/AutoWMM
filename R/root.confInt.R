#' @title root.confInt
#' @description Estimation helper function:
#'  Method for generating confidence interval for root node. Assumes unlogged
#'  input data and normally distributed logged data.  Completes conversion
#'  internally.
#' @param tree A makeTree object
#' @param int.type A string specifying interval type, passed from the wmmTree
#'  function.
#' @return Returns a confidence interval for the root population size estimate
#'  in un-logged form.
#' @examples \donttest{
#'  message("note - longer run time example")
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  Zhats <- wmmTree(tree, sample_length = 3)
#'  root.confInt(tree)
#' }
#' @export
#' @importFrom dplyr "select"
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom stats quantile
#' @importFrom stats cov
#' @import data.tree

root.confInt <- function(tree, int.type='quantiles'){
  # extract sample values for each path (column) (only leaves are required for
  # because that is where the estimates are stored)
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

  # get mean estimate
  m <- tree$Get('Estimate', filterFun = isRoot)

  # calculate interval based on int.type argument from wmmTree
  if(int.type =='var'){   #provides variance based interval
    message('using variance-weighted confidence interval - compare to quantiles')

    # calculate variance
    sig <- cov(log(x)) # covariance matrix
    #prec <- Inverse(sig) # precision matrix
    prec <- ginv(sig) # precision matrix
    e <- seq(1,1,length.out = dim(sig)[1])
    den.w <- t(e)%*%prec%*%e
    var <- 1/den.w[1]

    # generate endpoints
    lc <- log(m)-2*sqrt(var)
    uc <- log(m)+2*sqrt(var)
    lc <- exp(lc)
    uc <- exp(uc)
  }else if(int.type == 'cox'){        # provides cox interval
    message('using cox interval for uncertainty - compare to quantiles')

    # get weights
    weights <- as.matrix(ko.weights(tree))

    # get logN for each sample
    logrootEsts <- as.data.frame(as.matrix(x) %*% t(weights))

    # calculate sample variance of estimates
    nsamps <- length(logrootEsts)
    s2 <-sum((logrootEsts - log(m))^2)
    s2 <- s2/(nsamps-1)

    # generate interval endpoints
    m <- exp(log(m) + s2/2)
    lc <- m - 1.96*sqrt((s2/nsamps) + (s2^2/(2*(nsamps-1))))
    uc <- m + 1.96*sqrt((s2/nsamps) + (s2^2/(2*(nsamps-1))))
  }else{        # provides central 95% using quantiles
    # get weights
    weights <- as.matrix(ko.weights(tree))

    # get logN for each sample
    logrootEsts <- as.data.frame(as.matrix(x) %*% t(weights))

    uc <- exp(quantile(logrootEsts,0.975,na.rm = TRUE))
    lc <- exp(quantile(logrootEsts,0.025,na.rm = TRUE))
  }

  # return printed interval
  int <- c(round(lc),round(uc))
  names(int) <- c("lower", "upper")
  return(int)
}
