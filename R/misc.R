.findClade <- function (lineages) {
  # for a list of lineages, find the clade shared by all
  subj <- lineages[[1]]
  for (i in 2:length (lineages)) {
    query <- lineages[[i]]
    subj <- subj[subj %in% query]
  }
  subj[length (subj)]
}

cTrees <- function (treelist, tree) {
  ## Concatenate ape's phylo and multiPhylo classes
  ##  and return multiPhylo
  # check class types, and make sure they can be used
  ebool.multi <- any (class (tree) == 'multiPhylo')
  ebool.phylo <- any (class (tree) == 'phylo')
  if (!any (c (ebool.multi,ebool.phylo))) {
    stop ('Tree is neither phylo nor multiPhylo')
  }
  lbool.multi <- any (class (treelist) == 'multiPhylo')
  lbool.phylo <- any (class (treelist) == 'phylo')
  if (!any (c (lbool.multi,lbool.phylo))) {
    stop ('Treelist is neither phylo nor multiPhylo')
  }
  # container
  trees <- list ()
  # unpack each and add to container
  if (lbool.multi) {
    for (i in seq(length = length (treelist))) {
      trees <- c (trees, list (treelist[[i]]))
    }
  } else {
    trees <- c (trees, list (treelist))
  }
  if (ebool.multi) {
    for (i in seq(length = length (tree))) {
      trees <- c (trees, list (tree[[i]]))
    }
  } else {
    trees <- c (trees, list (tree))
  }
  # change class of container and return
  class (trees) <- 'multiPhylo'
  trees
}