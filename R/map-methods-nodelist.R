# TODO:
# adapt removeTip
# adapt addTip
# create methods for NodeList to rename
# test each function here

# All hidden functions for map-methods-visible.R using NodeList class
.addTip__mapNames <- function (paraenv, tips, qry) {
  new_name <- as.character (qry$search.name)
  max.age <- qry$max.age
  min.age <- qry$min.age
  parent <- paraenv[['getParent']](paraenv$grow_tree, tips=tips)
  edges <- paraenv[['getEdges']](paraenv$grow_tree, node=parent)
  ages <- paraenv[['getAge']](paraenv$grow_tree, edge=edges)
  pull <- ages$max.age >= max.age
  edges <- edges[pull]
  random <- sample (1:length (edges), size=1)
  age_range <- ages[random, ]
  node_age <- runif (n=1, min=age_range[1, 'max.age'],
                     max=age_range[1, 'max.age'])
  tip_age <- runif (n=1, min=min_age, max=node_age)
  tree <- addTip (tree, edge=edge, tip_name=new.name,
                  node_age=node_age,
                  tip_age=tip_age)
  return (tree)
}
.extract__mapNames <- function (tree, names) {
  # return tree representing names
  if (class (tree) == 'multiPhylo') {
    tree <- tree[[sample (1:length (tree), 1)]]
  }
  if (class (tree) == 'phylo') {
    tips <- tree$tip.label
  } else {
    tips <- tree@tips
  }
  if (sum (names %in% tips) > 1) {
    rm_tips <- tips[!tips %in% names]
    res_tree <- removeTip (tree, tip=rm_tips)
    return (res_tree)
  } else {
    warning ('Too few names could be mapped to tree')
    return (NA)
  }
}
.return__mapNames <- function (tree, names, iterations) {
  # if mapNames returns tree early, return an object
  #  that was expected i.e. NodeList, phylo or multiphylo
  tree <- .extract_mapNames (tree, names)
  if (iterations == 1) {
    return (tree)
  } else {
    trees <- list ()
    for (i in 1:iterations) {
      trees <- c (trees, list (tree))
    }
    if (class (tree) == 'phylo') {
      class (trees) <- 'multiPhylo'
    }
    return (trees)
  }
}
.clean__mapNames <- function (trees) {
  # drop _ in names of trees in a multiphylo
  # remove node labels
  # set branch lengths to 1 if no branch lengths
  .phylo <- function (tree) {
    tree$tip.label <- gsub ('_', ' ', tree$tip.label)
    tree$node.label <- NULL
    if (is.null (tree$edge.length)) {
      tree$edge.length <- rep (1, nrow (tree$edge))
    }
  }
  .drop <- function (i) {
    tree <- trees[[i]]
    .phylo (tree)
    res <<- c (res, list (tree))
  }
  if (class (trees) == 'multiPhylo') {
    res <- list ()
    m_ply (.data=data.frame (i=1:length (trees)), .fun=.drop)
    class (res) <- 'multiPhylo'
    return (res)
  } else if (class (tree) == 'phylo') {
    trees <- .phylo (trees)
  } else {
    if (any (grepl ('_', tree@tips))) {
      tips <- gsub ('_', ' ', tree@tips)
    }
    tree <- rename (tree, tree@tips, tips)
  }
  return (trees)
}
.getNames__mapNames <- function(trees) {
  # get tip names from a multiphylo
  if (class (trees) == 'multiPhylo') {
    .get <- function (i) {
      res <<- c (res, trees[[i]]$tip.label)
    }
    res <- NULL
    m_ply (.data=data.frame (i=1:length (trees)), .fun=.get)
    res <- unique (res)
  } else if (class (trees) == 'phylo') {
    res <- trees$tip.label
  } else {
    res <- trees@tips
  }
  return (res)
}