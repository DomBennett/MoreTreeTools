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

#' @name getSister
#' @title Return sister node
#' @description Return sister node of node given
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getSister <- function (tree, node) {
  inner.node <- tree$edge[match(node, tree$edge[,2]),1]
  # find all nodes and edges descending from inner node
  d.edges <- which(tree$edge[,1] %in% inner.node)
  d.nodes <- tree$edge[d.edges, 2]
  # remove starting node
  nearest.node <- d.nodes[!d.nodes %in% node][1]
  nearest.node
}

#' @name getEdges
#' @title Return edges
#' @description Return all descendant edges from a node,
#'  or all edges connected to tips based on three different methods
#' @details If node specified, all edges that descend from specified node are returned
#' 
#' Otherwise, if tips are specified edges are returned based on one of three methods:
#' Type:
#' 
#' 1 -- phylogeny consisting solely of the tips, default
#' 
#' 2 -- edges from taxon tips to terminal node
#' 
#' 3 -- edges unique to tips
#' @template base_template
#' @param node number indicating node
#' @param tips tip labels or numbers
#' @param type 1, 2 or 3 indicating which method to use
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getEdges <- function (tree, node = NA, tips = NA, type = 1) {
  if (!is.na (node)) {
    ## Find all edges from given node to tips
    edges <- c ()
    while (TRUE) {
      bool <- tree$edge[ ,1] %in% node
      if (sum (bool) > 1) {
        node <- tree$edge[bool, 2]
        edges <- c (edges, which (bool))
      } else {
        break
      }
    }
    return (edges)
  }
  # otherwise, use tips...
  if (!type %in% c (1,2,3)) {
    stop ("Type must be an integer: 1, 2 or 3.")
  }
  if (length (tips) == length (tree$tip.label)){
    return (tree$edge)
  }
  if (type == 1 & length (tips) == 1){
    stop ("length (tips) == 1 :
         Cannot return a single edge for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  if (is.character (tips)) {
    # if tips are labels
    edges <- match (match (tips, tree$tip.label), tree$edge[,2])
  } else {
    # ... else they're numbers
    edges <- match (tips, tree$edge[,2])
  }
  end.nodes <- tree$edge[edges, 1]
  term.node <- length (tree$tip.label) + 1
  if (all (end.nodes %in% term.node)) {
    return (edges)
  } else {
    if (type == 3){
      while (any (duplicated (end.nodes))){
        start.node <- end.nodes[duplicated(end.nodes)][1]
        if (sum (tree$edge[,1] %in% start.node) == sum (end.nodes %in% start.node)){
          edge <- match (start.node, tree$edge[,2])
          end.node <- tree$edge[edge,1]
          edges <- c(edges, edge)
          end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        } else {
          end.nodes <- end.nodes[end.nodes != start.node]
        }
      }
    } else {
      while (TRUE){
        end.nodes <- sort (end.nodes, TRUE)
        start.node <- end.nodes[1]
        edge <- match (start.node, tree$edge[,2])
        end.node <- tree$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        if (type == 2){
          if (sum (term.node == end.nodes) == length (end.nodes)){
            break
          }
        } else {
          if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
            break
          }
        }
      }
    }
    return (edges)
  }
}

#' @name getNodes
#' @title Return ascendant nodes
#' @description Return all nodes that connect the specified node to the root
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getNodes <- function (tree, node) {
  ## Find all nodes from given node to root
  base.node <- length (tree$tip.label) + 1
  nodes <- c ()
  while (node != base.node) {
    node <- tree$edge[tree$edge[ ,2] == node,1]
    nodes <- c (nodes, node)
  }
  nodes
}

#' @name getClades
#' @title Return all children of each node
#' @description All children for each node in the tree are returned.
#' @details Returns a list of lists: 'clade.children' contains all the children of each clade,
#'  'clade.node' contains the clade number of each clade in clade.children.
#' @template base_template
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getClades <- function (tree) {
  if (is.null (tree$all.node)) {
    nodes <- 1:(tree$Nnode + length (tree$tip.label))
  } else {
    nodes <- phylo$all.node
  }
  clades <- mlply (.data = data.frame (node = nodes),
                   .fun = getChildren, tree)
  if (!is.null (tree$all.node.label)) {
    names (clades) <- tree$all.node.label
  }
  sizes <- ldply (.data = clades, .fun = length)[ ,2]
  clades <- clades[order (sizes, decreasing = TRUE)]
  nodes <- nodes[order (sizes, decreasing = TRUE)]
  list (clade.children = clades, clade.node = nodes)
}