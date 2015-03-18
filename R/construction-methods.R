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
  if (!is.null (tree$node.ages)) {
    tree.node.ages <- tree$node.ages
    node.1.age <- tree.node.ages[node.1]
  } else {
    node.1.age <- getAge (tree, node=node.1)
  }
  new.node.edge.length <- node.1.age - node.age
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
    tree$node.label <- insert (tree$node.label, node.label,
                                   new.node - length (tree$tip.label))
  }
  if (!is.null (tree$all.node.label)) {
    tree$all.node.label <- insert (tree$all.node.label,
                                   c (tip.name, node.label),
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
  if (!is.null (tree$all.node.label)) {
    int.node.label <- tree$all.node.label[int.node.i]
    int.node.label <- int.node.label[c(1, o)]
    tree$all.node.label[int.node.i] <- int.node.label
  }
  if (!is.null (tree$node.label)) {
    o <- 1 + rank (tree$edge[sndcol, 2])
    tree$node.label <- tree$node.label[c (1, o)]
  }
  if (!is.null (tree$node.ages)) {
    int.node.age <- tree.node.ages[int.node.i]
    int.node.age <- int.node.age[c(1, o)]
    tree$node.ages[int.node.i] <- int.node.age
  }
  tree$edge[sndcol, 2] <- newNb[tree$edge[sndcol, 2] - n] <-
    n + 2:tree$Nnode
  tree$edge[, 1] <- newNb[tree$edge[, 1] - n]
  tree
}

#' @name removeTip
#' @title Remove a tip from a phylogenetic tree
#' @description Remove a tip from a phylogenetic tree.
#' @details Return a tree with specified tip dropped. This function
#' differs from \code{ape}'s \code{drop.tip} as it works with
#' \code{MoreTreeTools}' node labelling convention.
#' @template base_template
#' @param tip.name of tip to be dropped
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

removeTip <- function (tree, tip.name) {
  ## TODO: better to work out a way to work with ape's
  ##  nodelabels rather than do all this!!
  edges.to.drop <- c ()
  # find connecting nodes and edges
  tip.node <- which (tree$tip.label == tip.name)
  tip.edge <- which (tree$edge[ ,2] == tip.node)
  edges.to.drop <- c (edges.to.drop, tip.edge)
  internal.node <- tree$edge[tip.edge, 1]
  # if internal node is root, there is no internal edge
  if (internal.node != length (tree$tip.label) + 1) {
    internal.edge <- which (tree$edge[ ,2] == internal.node)
    edges.to.drop <- c (edges.to.drop, internal.edge)
    # if internal node has more than 1 child edge, must
    #  add this length to the already existing edge
    if (sum (tree$edge[, 1] == internal.node) > 1) {
      corres.edge <- which (tree$edge[, 1] == internal.node)
      corres.edge <- corres.edge[corres.edge != tip.edge]
      corres.node <- tree$edge[corres.edge, 2]
      internal.edge.length <- tree$edge.length[internal.edge]
      tree$edge.length[corres.edge] <-
        tree$edge.length[corres.edge] + internal.edge.length
    }
  }
  # remove edges to drop
  new.edges <- tree$edge[-edges.to.drop, ]
  # re-number nodes
  # remove 1 from all nodes greater than or equal to
  #  the internal node
  new.edges[new.edges[ ,1] >= internal.node,1] <-
    new.edges[new.edges[ ,1] >= internal.node,1] - 1
  new.edges[new.edges[ ,2] >= internal.node,2] <-
    new.edges[new.edges[ ,2] >= internal.node,2] - 1
  # remove 1 from all nodes above tip node
  new.edges[new.edges[ ,1] > tip.node,1] <-
    new.edges[new.edges[ ,1] > tip.node,1] - 1
  new.edges[new.edges[ ,2] > tip.node,2] <-
    new.edges[new.edges[ ,2] > tip.node,2] - 1
  # if there is a tripartition ...
  #  subtract 1 to the internal node that once was
  # (ordering has changed as a result of dropping)
  trip.bool <- sum (new.edges[ ,1] == internal.node - 2) > 2
  if (trip.bool) {
    if (corres.node > tip.node) {
      bool.edge <- new.edges[ ,2] == corres.node - 1
      new.edges[bool.edge, 1] <-
        new.edges[bool.edge, 1] - 1
    }
    else {
      bool.edge <- new.edges[ ,2] == corres.node
      new.edges[bool.edge, 1] <-
        new.edges[bool.edge, 1] - 1
    }
  }
  # replace old with new
  tree$edge <- new.edges
  tree$edge.length <- tree$edge.length[-edges.to.drop]
  tree$tip.label <- tree$tip.label[-tip.node]
  # update Nnode
  tree$Nnode <- tree$Nnode - 1
  # tidy up new tree
  if (!is.null (attr (tree, "order"))) 
    attr(tree, "order") <- NULL
  if (!is.null (tree$all.node.label)) {
    tree$all.node.label <- tree$all.node.label[c (-tip.node,
                                          -internal.node)]
  }
  if (!is.null (tree$node.ages)) {
    tree.node.ages <- tree$node.ages[c (-tip.node,
                                        -internal.node)]
  }
  tree
}