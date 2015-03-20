#' @name reLoad
#' @title Change loading of tree
#' @description # Return a tree with loading changed by factor.
#' @details The functions works by adding and removing branch length
#' between parent and children branches. The amount by which the branch
#' lengths are changed is determiend by factor.
#' @template base_template
#' @param factor proportion by which loading is changed, positive values
#' increase loading, negative values decrease it.
#' @export
#' @examples
#' # Create balanced tree, plot positively and negatively reloaded
#' n <- 16
#' par (mfrow = c (2, 2), mar = c (1, 1, 1, 1))
#' tree <- compute.brlen (stree (n, 'balanced'))
#' plot (reLoad (tree, factor = -0.9), show.tip.label=FALSE, edge.width=2)
#' plot (reLoad (tree, factor = 0.9), show.tip.label=FALSE, edge.width=2)

reLoad <- function (tree, factor) {
  # Return a tree with the loading changed by factor
  downLoad <- function (node) {
    # decrease tree loading
    # get the preceding edge for node
    preceding.edge <- tree$edge[ ,2] == node
    if (any (preceding.edge)) {
      # get the following edges
      following.edges <- tree$edge[ ,1] == node
      # calculate amount to change the preceding and following edges
      change.by <- tree$edge.length[preceding.edge] * factor
      # add and remove by change.by
      tree$edge.length[preceding.edge] <-
        tree$edge.length[preceding.edge] - change.by
      tree$edge.length[following.edges] <-
        tree$edge.length[following.edges] + change.by
    }
    tree <<- tree
  }
  upLoad <- function (node) {
    # increase tree loading
    # get the following edges for node
    following.edges <- tree$edge[ ,1] == node
    if (any (following.edges)) {
      # get the preceding edges
      preceding.edge <- tree$edge[ ,2] == node
      # calculate amount to change the preceding and following edges
      change.by <- tree$edge.length[following.edges][1] * factor
      # add and remove by change.by
      tree$edge.length[preceding.edge] <-
        tree$edge.length[preceding.edge] + change.by
      tree$edge.length[following.edges] <-
        tree$edge.length[following.edges] - change.by
    }
    tree <<- tree
  }
  # list internal nodes
  internal.nodes <- length (tree$tip.label)+1:
    (length (tree$tip.label) +tree$Nnode)
  # if factor is negative, downLoad
  if (factor < 0) {
    factor <- abs (factor)
    m_ply (.data = data.frame (node = internal.nodes),
           .fun = downLoad)
  } else {
    # else upLoad
    m_ply (.data = data.frame (node = internal.nodes),
           .fun = upLoad)
  }
  tree
}

#' @name reBalance
#' @title Change balance of tree
#' @description Return a tree with balance changed.
#' @details The functions switches tips in tree to make the tree
#' more or less balanced. If steps is positive, the tree will be
#' made more balanced else more unbalanced. Steps determines the
#' number of tip switches in a balanced or unbalanced direction.
#' @template base_template
#' @param steps number of tip switches to perform
#' @export
#' @examples
#' # Make an unbalanced tree balanced
#' n <- 16
#' tree <- compute.brlen (stree (n, 'left'))
#' plot (tree, show.tip.label=FALSE, edge.width=2, main='Unbalanced')
#' # +11 steps to make tree balanced
#' tree <- reBalanced (tree, 11)
#' plot (tree, show.tip.label=FALSE, edge.width=2, main='Rebalanced')

reBalance <- function (tree, steps) {
  run <- function (i) {
    rtts <- diag (vcv.phylo (tree))
    dists <- rtts - log2 (getSize (tree))
    if (balancing) {
      drop <- names (rtts) [dists == max (dists)][1]
      add <- names (rtts) [dists == min (dists)][1]
    } else {
      add <- names (rtts) [dists == max (dists)][1]
      tds <- cophenetic.phylo (tree)
      drop <- rownames (tds)[which.max (tds[add, ])]
    }
    # remove drop.tip
    tree <- drop.tip (tree, tip = drop)
    # add new tip
    target.edge <- which (tree$edge[ ,2] == which (tree$tip.label == add))
    tree <- addTip (tree = tree, edge = target.edge,
                    tip.name = drop, node.age = 0)
    # reset edge.lengths all to 1
    tree$edge.length <- rep (1, nrow (tree$edge))
    tree <<- tree
  }
  # tree must be bifurcating
  if (tree$Nnode != (getSize (tree) - 1)) {
    stop ('Tree must be bifurcating')
  }
  if (steps > 0) {
    # if steps positive, balanced is True
    balancing <- TRUE
  } else if (steps < 0) {
    # if steps negative, balanced is True
    balancing <- FALSE
  } else {
    # if steps is 0, return tree
    return (tree)
  }
  # set all edges to 1 -- so branch is comensurate with taxonomic distance
  tree$edge.length <- rep (1, nrow (tree$edge))
  m_ply (.data = data.frame (i = 1:abs (steps)), .fun = run)
  tree
}