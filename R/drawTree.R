#' @title drawTree
#' @description Visualize tree with descriptions and probabilities; can be used
#'  pre-WMM analysis
#' @param tree A makeTree object
#' @param probs A logical with default TRUE to specify whether to display
#'  probabilities on branches
#' @param desc A logical with default TRUE to specify whether to display
#'  node descriptions
#' @return Returns a descriptive tree plot
#' @examples
#' data(treeData1)
#' tree <- makeTree(treeData1)
#' drawTree(tree)
#' @export
#' @import data.tree
#' @import DiagrammeR

drawTree <- function(tree, probs=TRUE, desc = TRUE){
  SetGraphStyle(tree,scale=2)
  if(probs){
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2,
                 label = function(node) round(mean(node$probability), 2))
  }else{
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2)
  }
  tree$Do(function(node){
    if(!is.null(node$Description) && desc){ # incorporate description on visual diagram
      SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                   label = function(node) node$Description)
    }else{
      SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                   label = function(node) node$to)
    }
  })
  plot(tree)
}

