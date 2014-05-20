#' @name getChildren
#' @title Return descendant species from a node
#' @description Return all tip labels that descend from a specifed node.
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @param display if true, a tree will be plotted with specifed node highlighted
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getChildren <- function(tree, node, display = FALSE) {
  if (!is.numeric (node)) {
    stop("Node is not numeric!")
  }
  if (node > tree$Nnode + length (tree$tip.label)) {
    stop("Node is greater than the number of nodes in tree!")
  }
  if (node <= length (tree$tip.label)) {
    term.nodes <- node
  } else {
    term.nodes <- vector ()
    temp.nodes <- node
    while (length (temp.nodes) > 0) {
      connecting.nodes <- tree$edge[tree$edge[,1] %in% temp.nodes, 2]
      term.nodes <- c(term.nodes, connecting.nodes[connecting.nodes <=
                                                     length (tree$tip.label)])
      temp.nodes <- connecting.nodes[connecting.nodes > length(tree$tip.label)]
    }
  }
  children <- tree$tip.label[term.nodes]
  if (display) {
    tip.cols <- ifelse(tree$tip.label %in% children, "black", "grey")
    plot.phylo(tree, tip.color = tip.cols, show.tip.label = TRUE)
    nodelabels("node", node)
  }
  return (children)  
}

#' @name getAge
#' @title Return age of node
#' @description Return age of tip or internal node (so long as the tree provided is time calibrated with branch lengths).
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getAge <- function (tree, node) {
  term.node <- length (tree$tip.label) + 1
  # if running for multiple nodes, define tree.age first to
  #  save time
  if (!exists ('tree.age')) {
    tree.age <- max (diag (vcv.phylo (tree)))
  }
  # if it's the root node, return tree.age
  if (term.node == node) {
    return (tree.age)
  }
  # if it's a tip return 0
  if (node < length (tree$tip.label) & is.ultrametric (tree)) {
    return (0)
  }
  # else find all its edge.lengths and subtract from tree.age
  edges <- c ()
  while (node != term.node) {
    edges <- c (edges, which (tree$edge[ ,2] == node))
    node <- tree$edge[tree$edge[ ,2] == node, 1]
  }
  return (tree.age - sum (tree$edge.length[edges]))
}