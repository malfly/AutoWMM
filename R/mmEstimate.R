#' @title mmEstimate
#' @description Helper function:
#'  Performs the multiplier method from a single terminal node (o) and
#'  returns the root estimate given that path, and the probabilities
#'  of each branch on that path.
#' @param o A node from a makeTree object
#' @return Returns node with modified attributes
#' @examples \donttest{
#'  data(treeData1)
#'  tree <- makeTree(treeData1)
#'  mmEstimate(tree$A)
#'  mmEstimate(tree$B)
#'  tree$A$targetEst
#'  tree$B$targetEst
#' }
#' @export
#' @import data.tree

mmEstimate <- function(o){
  if(o$isLeaf){
    # first, generate an estimate of the root of the leaf
    estimate <- o$Count*(1/o$probability)
    current_node <- o$parent # set the current node to now be that parent

    # while current node exists (ie the path back to the root is not yet
    # done), keep moving backwards and generating estimates for successively
    # higher nodes, until you reach the root, and generate a root estimate
    # this while loop performs the multiplicative back-calculation of one path
    # which goes from node o back to the root
    while(!is.null(current_node$parent)){
      estimate <- estimate *(1/current_node$probability)
      current_node <- current_node$parent
    }
    o$targetEst <- estimate  # o$targetEst is the path-specific estimate
  }

  # The following assures an estimate of the root if not generated if
  # the node with marginal count is not a leaf.
  else{
    o$targetEst <- NA
  }
}
