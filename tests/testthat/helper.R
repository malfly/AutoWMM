# helper.R

create_tree <- function() {
  return(makeTree(treeData1))
}

create_Zs <- function() {
  tree_obj <- makeTree(treeData1)
  return(wmmTree(tree_obj, sample_length = 2))
}
