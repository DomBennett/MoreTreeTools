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
    node.1.age <- getAge (tree, node=node.1)$age
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
#' @description Remove tips from a phylogenetic tree
#' @details Based on \code{ape}'s \code{drop.tip} but with
#' \code{preserve.age}. This will ensure
#' the root to tip distance of the tree is maintained by adding length
#' to next immediate internal edge whenever the removal of a tip
#' leads to a reduction in a tree's age. If root node is lost, then
#' lost edge length is added to \code{root.edge} slot.
#' @template base_template
#' @param tips vector of tip names
#' @param preserve.age maintain the root to tip distance of the tree
#' @export
#' @examples
#' tree <- rtree (10)
#' par (mfrow=c (1,3))
#' plot (tree)
#' plot (removeTip (tree, tips=c ('t1', 't2', 't3'), preserve.age=FALSE))
#' plot (removeTip (tree, tips=c ('t1', 't2', 't3'), preserve.age=TRUE), root.edge=TRUE)

removeTip <- function (tree, tips, preserve.age) {
  # internal
  .run <- function (tip.name) {
    # names description:
    # - tree.age, age of the tree from root to tip
    # - root node, the node number of the root of the tree
    # - max.node, the maximum internal node number of the tree
    # - edges.to.drop, edge numbers that will be dropped (max 2, target and parent edges)
    # - target.node, node number of tip.name
    # - target.edge, connecting edge of tip node
    # - node.1, first internal node connecting to target.node and sister.node
    # - node.2, second internal node connecting to tree to node.1
    # - sister.node, node of sister
    # - sister.edge, edge of sister
    # - parent.edge, edge of node.2 and node.1
    # - .ori, original data provided
    # - .mod, modified without missing elements
    # - .new, modified with missing elements 
    #print (tip.name)
    # base vars
    tree.age <- getSize (tree, 'rtt')
    edges.ori <- edges.mod <- tree$edge
    root.node <- length (tree$tip.label) + 1
    max.node <- tree$Nnode + length (tree$tip.label)
    target.node <- which (tree$tip.label == tip.name)
    target.edge <- which (tree$edge[ ,2] == target.node)
    edges.to.drop <- target.edge
    node.1 <- tree$edge[target.edge, 1]
    # get sister
    sister.edge <- which (tree$edge[, 1] == node.1)
    sister.edge <- sister.edge[sister.edge != target.edge]
    sister.node <- tree$edge[sister.edge, 2]
    if (node.1 != root.node) {
      # get parent
      parent.edge <- which (tree$edge[ ,2] == node.1)
      node.2 <- tree$edge[parent.edge,1]
      edges.to.drop <- c (edges.to.drop, parent.edge)
      # update length of sister.edge
      parent.edge.length <- tree$edge.length[parent.edge]
      sister.edge.length <- tree$edge.length[sister.edge]
      target.edge.length <- tree$edge.length[target.edge]
      sister.edge.length.new <- sister.edge.length +
        parent.edge.length
      if (target.edge.length > sister.edge.length &
            preserve.age) {
        sister.edge.length.new <- sister.edge.length.new +
          (target.edge.length - sister.edge.length)
      }
      tree$edge.length[sister.edge] <- sister.edge.length.new
      # replace node.1 with node.2 for sister.edge
      edges.mod[sister.edge, 1] <- node.2
      # remove edges to drop
      edges.new <- edges.mod[-edges.to.drop, ]
      # re-number nodes
      # remove 1 from all nodes greater than or equal to
      #  the node.1
      edges.new[edges.new[ ,1] >= node.1, 1] <-
        edges.new[edges.new[ ,1] >= node.1, 1] - 1
      edges.new[edges.new[ ,2] >= node.1, 2] <-
        edges.new[edges.new[ ,2] >= node.1, 2] - 1
      # remove 1 from all nodes above tip node
      edges.new[edges.new[ ,1] > target.node,1] <-
        edges.new[edges.new[ ,1] > target.node,1] - 1
      edges.new[edges.new[ ,2] > target.node,2] <-
        edges.new[edges.new[ ,2] > target.node,2] - 1
    } else {
      # if node is root, there is no parent edge, must remove sister
      edges.to.drop <- c (edges.to.drop, sister.edge)
      # remove edges to drop
      edges.new <- tree$edge[-edges.to.drop, ]
      # re-number nodes
      edges.new <- edges.new - 1
      edges.new[ ,1] <- edges.new[ ,1] - 1
      edges.new[edges.new[ ,2] >= node.1, 2] <-
        edges.new[edges.new[ ,2] >= node.1, 2] - 1
      if (preserve.age) {
        sister.node.age <- getAge (tree, node=sister.node)[ ,'age']
        tree$root.edge <- tree.age - sister.node.age
      }
    }
    # replace old with new
    tree$edge <- edges.new
    tree$edge.length <- tree$edge.length[-edges.to.drop]
    if (!is.null (tree$edge.label)) {
      tree$edge.label <- tree$edge.label[-edges.to.drop]  
    }
    tree$tip.label <- tree$tip.label[-target.node]
    # update Nnode
    tree$Nnode <- tree$Nnode - 1
    # tidy up new tree
    if (!is.null (attr (tree, "order"))) {
      attr(tree, "order") <- NULL
    }
    #storage.mode (tree$edge) <- "integer"
    #tree.c <- collapse.singles (tree)
    tree <<- tree
  }
  if (!all (tips %in% tree$tip.label)) {
    stop ('Not all tips are in the tree')
  }
  if (length (tips) >= (length (tree$tip.label) - 2)) {
    stop ('Removing too many tips, smallest resulting tree is 3 tips')
  }
  loop.data <- data.frame (tip.name=tips, stringsAsFactors=FALSE)
  plyr::m_ply (.data=loop.data, .fun=.run)
  tree
}

#' @name collapseTips
#' @title Remove all tips smaller than min.length
#' @description Return a phylo object with tip edges less than min.length
#' removed.
#' @details Tips are removed, but the tree age is preserved. Returned
#' trees are therefore not the same phylogenetic tree as that given.
#' The purpose of the function is to generate representations of
#' trees with fewer tips for plotting purposes.
#' 
#' 
#' Set iterative as TRUE to keep removing all tip edges less than min.length.
#' @template base_template
#' @param min.length smallest length for a tip edge
#' @param iterative run until no tip edges less than min.length, default TRUE
#' @export
#' @examples
#' tree <- rtree (100)
#' par (mfrow = c (1,2))
#' plot (tree, show.tip.label=FALSE, main='Full tree')
#' tree.age <- getSize (tree, 'rtt')
#' ctree <- collapseTips (tree, min.length=tree.age*0.1)
#' plot (ctree, show.tip.label=FALSE, main='Collapsed tree')

collapseTips <- function (tree, min.length, iterative=TRUE) {
  while (TRUE) {
    tip.edges <- tree$edge[ ,2] < getSize (tree)
    drop.tip.edges <- (tree$edge.length < min.length) &
      tip.edges
    to.drop <- tree$tip.label[tree$edge[drop.tip.edges,2]]
    if (length (to.drop) == 0) {
      break
    }
    tree <- try (removeTip (tree, tip=to.drop, preserve.age=TRUE),
                 silent=TRUE)
    if (grepl ('Removing too many tips', tree[[1]][1])) {
      stop ('min.length is too great, too many tips are being dropped')
    }
    if (!iterative) {
      break
    }
  }
  tree
}