#' @name addTip
#' @title Add a new tip to a phylogenetic tree
#' @description A new tip is added to an existing phylogenetic tree by
#' specifying an edge/branch in the tree, name of the new tip and ages of
#' the new tip and node.
#' @details Return a tree with an extra tip. The new tip is created by
#' finding the \code{edge} in the tree, adding a new node at \code{node.age}
#' and then by adding a new branch to \code{tip.age}.
#' The function is modified code from \code{multi2di} from \code{ape}.
#' N.B. \code{ape} uses node.label to label all internal nodes, \code{MoreTreeTools}
#' uses it to labels all nodes including tip nodes.
#' @template base_template
#' @param edge tree branch where tip will be added
#' @param tip.name name for the new tip to be added
#' @param node.age age of the new node to be created
#' @param tip.age age of the new tip to be added (default 0)
#' @param node.label label for the new
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

addTip <- function (tree, edge, tip.name, node.age,
                    tip.age = 0, node.label) {
  insert <- function (target.vector, source.vector, index) {
    #http://stackoverflow.com/questions/1493969/how-to-insert-elements-into-a-vector
    index <- index - 1
    positions <- c (seq_along (target.vector), index + 0.5)
    res <- c (target.vector, source.vector)
    res[order (positions)]
  }
  # nodes that make edge
  node.1 <- tree$edge[edge, 1]
  node.2 <- tree$edge[edge, 2]
  new.tip.edge.length <- node.age - tip.age
  if (new.tip.edge.length < 0) {
    stop ("Node age must be greater than tip age.")
  }
  tree.node.ages <- getAge (tree)[ ,2]
  new.node.edge.length <- tree.node.ages[node.1] - node.age
  if (new.node.edge.length < 0) {
    stop ("Node age is greater than edge.")
  }
  # build a new edge matrix
  edges.to.replace <- c(edge, getEdges (tree, node.2))
  new.edge.lengths <- c (new.node.edge.length, new.tip.edge.length,
                         tree$edge.length[edge] - new.node.edge.length,
                         tree$edge.length[edges.to.replace][-1])
  target <- node.1 + 1
  tree$edge[tree$edge > length (tree$tip.label)] <- 
    tree$edge[tree$edge > length (tree$tip.label)] + 1
  tree$Nnode <- tree$Nnode + 1
  new.tip.labels <- c (tree$tip.label, tip.name)
  new.tip <- length (new.tip.labels)
  new.node <- new.tip + tree$Nnode
  new.edges <- matrix (nrow = length (edges.to.replace) + 2, ncol = 2)
  new.edges[1, ] <- c (target, new.node)
  new.edges[2, ] <- c (new.node, new.tip)
  for (i in 1:length (edges.to.replace)) {
    new.edge <- tree$edge[edges.to.replace[i],]
    new.edge[new.edge == target] <- new.node
    new.edges[i+2, ] <- new.edge
  }
  # replace old with new
  tree$edge <- rbind (tree$edge[-edges.to.replace, ], new.edges)
  tree$edge.length <- c (tree$edge.length[-edges.to.replace],
                          new.edge.lengths)
  tree$tip.label <- new.tip.labels
  # tidy up new tree (multi2di code)
  if (!is.null (attr (tree, "order"))) 
    attr(tree, "order") <- NULL
  if (!is.null (tree$node.label)) {
    tree$node.label <- insert (tree$node.label, c (tip.name, node.label),
                               c (new.tip, new.node))
  }
  if (!is.null (tree$node.ages)) {
    tree.node.ages <- insert (tree.node.ages, c (tip.age, node.age),
                              c (new.tip, new.node))
  }
  tree <- reorder (tree)
  newNb <- integer (tree$Nnode)
  n <- length (tree$tip.label)
  newNb[1] <- n + 1L
  sndcol <- tree$edge[, 2] > n
  o <- 1 + rank (tree$edge[sndcol, 2])
  int.node.i <- (length (tree$tip.label) + 1):
    (length (tree$tip.label) + tree$Nnode)
  if (!is.null (tree$node.label)) {
    int.node.label <- tree$node.label[int.node.i]
    int.node.label <- int.node.label[c(1, o)]
    tree$node.label[int.node.i] <- int.node.label
  }
  if (!is.null (tree$node.ages)) {
    int.node.age <- tree.node.ages[int.node.i]
    int.node.age <- int.node.age[c(1, o)]
    tree$node.ages[int.node.i] <- int.node.age
  }
  tree$edge[sndcol, 2] <- newNb[tree$edge[sndcol, 2] - n] <- n + 
    2:tree$Nnode
  tree$edge[, 1] <- newNb[tree$edge[, 1] - n]
  tree
}