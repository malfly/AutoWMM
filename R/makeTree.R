#' @title makeTree
#' @description Assuming a specific structure, create a tree with the following columns:
#'  from (node label), to (node label), Estimate (+ integer), Total (+ integer),
#'  and Count (for terminal nodes with marginal counts).
#'  'from' and 'to' describe the edge for that row of data, where 'Estimate'
#'  and 'Total' are assumed to come from surveys of size 'Total' (a sample of
#'  the population at node 'from'), and observe 'Estimate' number of those
#'  individuals at 'Total' which move to the node described by 'to'.
#'  'Estimate' and 'Total' columns are used for branching probabilities only.
#'  'Count' column is NA for rows where 'to' nodes are not leaves; and also
#'  for all leaves without a marginal count.
#'  A Population (logical) column is not needed, but can be added if 'Estimate'
#'  and 'Total' come from population numbers, rather than samples.
#'  A 'Description' column (string) is also possible to include if particulars
#'  are desired on the tree diagram.
#'  'TerminalCount' (binary) will be created for functional purposes, where
#'  marginal counts are included on leaves.
#' @param data A dataframe object
#' @return Returns a makeTree object
#' @examples
#' data(treeData1)
#' tree <- makeTree(treeData1)
#' @export
#' @import data.tree

makeTree <- function(data){
  # check structure of data here
  if(!is.data.frame(data)){
    message('data must be dataframe type.')
  }
  # if no column x, print 'no column x'
  if(is.null(data$from)){
    message('data must have \'from\' column')
  }
  if(is.null(data$to)){
    message('data must have \'to\' column')
  }
  if(is.null(data$Estimate)){
    message('data must have \'Estimate\' column')
  }
  if(is.null(data$Total)){
    message('data must have \'Total\' column')
  }
  if(is.null(data$Count)){
    message('data must have \'Count\' column')
  }
  if(all(is.na(data$Total))){
    message('\'Total\' values cannot all be NA')
  }
  if(all(is.na(data$Estimate))){
    message('\'Estimate\' values cannot all be NA')
  }
  # if no 'Population' column exists, assume all values are samples
  if(is.null(data$Population)){
    data$Population <- rep(FALSE, times = dim(data)[1])
    message('WARNING: No \'Population\' column exists. We assume all
          values are sample estimates')
  }

  # Now create tree structure
  tree <- FromDataFrameNetwork(data)
  tree$Do(function(node){
    node$probability <- round(node$Estimate/node$Total,2)
    # Create a binary value called 'TerminalCount' which is true for leaves with
    # marginal count data only
    if(node$isLeaf & !is.null(node$Count)){
      node$TerminalCount <- TRUE
    }else{
      node$TerminalCount <- FALSE
    }
  })
  return(tree)
}
