#' @title countTree
#' @description Visualize post-wmmTree tree with root estimate and marginal counts
#'  Also displays average of probability samples on each branch
#' @param tree A makeTree object
#' @return Returns a tree plot
#' @examples \donttest{
#'   message("note - longer run time example")
#'   data(treeData1)
#'   tree <- makeTree(treeData1)
#'   Zhats <- wmmTree(tree, sample_length = 3)
#'   countTree(tree)
#'  }
#' @export
#' @import data.tree
#' @import DiagrammeR

countTree <- function(tree){
  # check if probability samples are empty everywhere - this indicates
  # weightedTree has not yet been used
  if(is.na(tree$Get('Estimate', filterFun = isRoot))){
    message('weightedTree() function has not yet been applied. conduct root estimation first.')
  }else{
    SetGraphStyle(tree,scale=2)
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2,
                 label = function(node) round(mean(node$probability_samples), digits = 2))
    SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                 label = function(node) if(isRoot(node)){round(node$Estimate)}else{node$Count})
    plot(tree)
  }
}
