#' Return descendant species from a node
#'
#' Return all tip labels that descend from a specifed node.
#' 
#' No details
#'
#' @template test_template
#' @param node number indicating node
#' @param display if true, a tree will be plotted with specifed node highlighted
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

nodeDescendants <- function(phylo, node, display = FALSE) {
  if (!is.numeric(node)) {
    stop("node is not numeric!")
  }
  if (node > phylo$Nnode + length(phylo$tip.label)) {
    stop("node is greater than the number of nodes in phylo!")
  }
  if (node <= length(phylo$tip.label)) {
    term.nodes <- node
  } else {
    term.nodes <- vector()
    temp.nodes <- node
    while (length(temp.nodes) > 0) {
      connecting.nodes <- phylo$edge[phylo$edge[,1] %in% temp.nodes, 2]
      term.nodes <- c(term.nodes, connecting.nodes[connecting.nodes <= length(phylo$tip.label)])
      temp.nodes <- connecting.nodes[connecting.nodes > length(phylo$tip.label)]
    }
  }
  descendants <- phylo$tip.label[term.nodes]
  if (display) {
    tip.cols <- ifelse(phylo$tip.label %in% descendants, "black", "grey")
    plot.phylo(phylo, tip.color = tip.cols, show.tip.label = TRUE)
    nodelabels("node", node)
  }
  return (descendants)  
}