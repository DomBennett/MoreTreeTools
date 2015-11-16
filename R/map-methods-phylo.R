# Hidden methods for pinning and mapping names
# Placed here to prevent double-coding and make testing and debugging easier



# phylo methods
.mnAddTipWAge <- function (tree, tip.is, new.name,
                           min.age, max.age) {
  print(new.name)
  # add new tip
  if (length (tip.is) == 1) {
    edge <- which (tip.is == tree$edge[ ,2])
    edge.age <- getAge (tree, edge=edge)
    pull <- edge.age$max.age >= max.age
    if (!pull) {
      edge <- NA
    }
  } else {
    # randomly map new tip to any edge in the clade
    # represented by the matching tip.is
    # find the parent node of all the matching tips
    children <- tree$tip.label[tip.is]
    parent.node <- getParent (tree, tips=children)
    # get all descending
    edges <- getEdges(tree, node=parent.node)
    # INCLUDING supporting edge!!!
    edges <- c (edges, which (tree$edge[ ,2] ==parent.node))
    # check ages
    edge.ages <- getAge (tree, edge=edges)
    pull <- edge.ages$max.age >= max.age
    edges <- edges[pull]
    # choose edge at random
    if (length (edges) > 1) {
      edge <- sample (edges, size=1)
    } else {
      edge <- edges[1]
    }
  }
  # exit if no suitable edge found
  if (!is.numeric (edge)) {
    return (tree)
  }
  # random tip and node age
  age.range <- getAge (tree, edge=edge)
  node.age <- runif (n=1, min=age.range[1, 'max.age'],
                     max=age.range[1, 'max.age'])
  tip.age <- runif (n=1, min=min.age, max=node.age)
  tree <- addTip (tree, edge=edge, tip.name=new.name, node.age=node.age,
                  tip.age=tip.age)
  return (tree)
}
.mnAddTip <- function (tree, tip.is, new.name) {
  # add new tip
  if (length (tip.is) == 1) {
    # add to the pendant edge at a random time point
    edge <- which (tip.is == tree$edge[ ,2])
  } else {
    # randomly map new tip to any edge in the clade
    # represented by the matching tip.is
    # find the parent node of all the matching tips
    children <- tree$tip.label[tip.is]
    parent.node <- getParent (tree, tips=children)
    # get descending and supporting edges
    edges <- which (tree$edge[ ,1] == parent.node)
    edges <- c (edges, which (tree$edge[ ,2] == parent.node))
    # choose edge at random, bias sample based on branch length
    edge <- sample (edges, size=1, prob=tree$edge.length[edges])
  }
  # random node age somewhere on edge
  age.range <- getAge (tree, edge=edge)
  node.age <- runif (n=1, min=age.range[1, 'min.age'],
                     max=age.range[1, 'max.age'])
  # random tip age if not ultrametric
  if (is.ultrametric (tree)) {
    tip.age <- 0
  } else {
    # use random tip age in shared clade
    possibles <- getAge (tree, node=tip.is)$age
    possibles <- possibles[possibles < node.age]
    tip.age <- runif (n=1, min=min(possibles), max=max(possibles))
  }
  tree <- addTip (tree, edge=edge, tip.name=new.name, node.age=node.age,
                  tip.age=tip.age)
  return (tree)
}
.mnExtract <- function (tree, names) {
  # return tree representing names
  if (class (tree) == 'multiPhylo') {
    tree <- tree[[sample (1:length (tree), 1)]]
  }
  if (sum (names %in% tree$tip.label) > 1) {
    res.tree <- drop.tip (tree, tip = tree$tip.label[!tree$tip.label %in% names])
    return (res.tree)
  } else {
    warning ('Too few names could be mapped to tree')
    return (NA)
  }
}
.mnEarlyReturn <- function (tree, names, iterations) {
  # if mapNames returns tree early, return an object
  #  that was expect i.e. phylo or multiphylo
  tree <- .mnExtract (tree, names)
  if (iterations == 1) {
    return (tree)
  } else {
    trees <- list ()
    for (i in 1:iterations) {
      trees <- c (trees, list (tree))
    }
    class (trees) <- 'multiPhylo'
    return (trees)
  }
}
.mnClean <- function (trees) {
  # drop _ in names of trees in a multiphylo
  # remove node labels
  # set branch lengths to 1 if no branch lengths
  if (class (trees) == 'multiPhylo') {
    .drop <- function (i) {
      tree <- trees[[i]]
      tree$tip.label <- gsub ('_', ' ', tree$tip.label)
      tree$node.label <- NULL
      if (is.null (tree$edge.length)) {
        tree$edge.length <- rep (1, nrow (tree$edge))
      }
      res <<- c (res, list (tree))
    }
    res <- list ()
    m_ply (.data=data.frame (i=1:length (trees)), .fun=.drop)
    class (res) <- 'multiPhylo'
    return (res)
  } else {
    trees$tip.label <- gsub ('_', ' ', trees$tip.label)
    trees$node.label <- NULL
    if (is.null (trees$edge.length)) {
      trees$edge.length <- rep (1, nrow (trees$edge))
    }
    return (trees)
  }
}
.mnGetNames <- function(trees) {
  # get tip names from a multiphylo
  if (class (trees) == 'multiPhylo') {
    .get <- function (i) {
      res <<- c (res, trees[[i]]$tip.label)
    }
    res <- NULL
    m_ply (.data=data.frame (i=1:length (trees)), .fun=.get)
    res <- unique (res)
  } else {
    res <- trees$tip.label
  }
  return (res)
}

# nodelist methods