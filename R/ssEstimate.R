#' @title ssEstimate
#' @description Helper function:
#'  Performs the closed form calculation of variances and means based on a
#'  "single-source sibling" tree.  Engages when single.source = TRUE
#'  in wmmTree function.  See documentation for further details.
#' @param o A node from a makeTree object
#' @return Returns node with modified attributes
#' @examples \donttest{
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  Zhats <- wmmTree(tree, sample_length = 3)
#'  ssEstimate(tree$B)
#' }
#' @export
#' @import data.tree

ssEstimate <- function(o){
  if(o$isLeaf){
    # first, multiply the marginal count by first branch segment's contribution to product
    meanCalc <- o$Count*((o$Total + o$Estimate - 1)/(o$Estimate - 1))
    varCalc <- (o$Count)^2*(((o$Total + o$Estimate - 1)*(o$Total + o$Estimate - 2))/
                              ((o$Estimate - 1)*(o$Estimate - 2)))
    current_node <- o$parent # set the current node to now be that parent

    # while current node exists (ie the path back to the root is not yet
    # done), keep moving backwards and multiplying by the quotient for successively
    # higher nodes, until you reach the root, and complete the mean and var calculation
    # this while loop performs the multiplicative back-calculation of one path
    # which goes from node o back to the root
    while(!is.null(current_node$parent)){
      meanCalc <- meanCalc*((current_node$Total + current_node$Estimate - 1)/(current_node$Estimate - 1))
      varCalc <- varCalc*(((current_node$Total + current_node$Estimate - 1)*
                             (current_node$Total + current_node$Estimate - 2))/
                            ((current_node$Estimate - 1)*(current_node$Estimate - 2)))
      current_node <- current_node$parent
    }
    o$targetEst <- meanCalc  # o$targetEst is the path-specific estimate
    o$variance <- varCalc
  }

  # The following assures an estimate of the root if not generated if
  # the node with marginal count is not a leaf.
  else{
    o$targetEst <- NA
    o$variance <- NA
  }
}


