#' @name mapNames
#' @title Map names to a tree
#' @description Map names to a tree with the option of searching online
#' to perform fuzzy-name matching.
#' @details This function firsts matches names to tip labels in tree, if not all names
#' are present in the tree and 'fuzzy' is true, it will then search the online taxonomic
#' names resolution service Global Names Resolver. Taxonomic IDs are then matched
#' between tip labels and names, if names are still not mapped to tree, names are mapped
#' to tips using random mapping. For example, if an unmapped name is a member of the same
#' genus as several species in the tree, the name will be mapped to the tree with any one
#' of these at random. In cases where an unmapped name can only be added to the tree by adding
#' an extra tip, then a new tip is added at random to the pendant edge of a suitable tip at
#' a random time point (it assumes the tree is time callibrated) -- see examples.
#' 
#' In cases where a large proportion of the tips are mapped using random placement, it may
#' be best to generate a distribution of trees that represent the range of possible phylogenetic
#' relationships for the names and tree given. If iterations is greater than 1, instead of a
#' single tree, a multiPhylo object of length iterations is returned representing this
#' range of possible mappings.
#' 
#' For non-time-callibrated or trees without branch lengths, the function using taxonomic
#' distance for random placement and returns a tree without branch lengths.
#' 
#' @template base_template
#' @param names vector of names of tips to be extracted
#' @param fuzzy boolean, if true will search Global Names Resolver online
#' @param datasource GNR datasource ID, default 4 for NCBI
#' @param iterations how many times to repeat?
#' @export
#' @examples
#' # bring in the catarrhines data
#' data ('catarrhines')
#' # we want to map these names to the catarrhines tree
#' names <- c ('Homo sapiens', 'Pongo pygmaeus', 'Gorilla gorila', 'Pan troglodytes',
#'             'Homo erectus', 'Homo neanderthalensis', 'Hylobates')
#' # 4 of the names are already in the tree (one has a spelling mistake) but the other 3
#' # do not exist and will be mapped randomly based on resolved taxonomic lineages of the
#' # names already in the tree and the names to be added.
#' hominid.map <- mapNames (tree=catarrhines, names=names, fuzzy=TRUE)
#' plot (hominid.map)

mapNames <- function (tree, names, fuzzy=TRUE, datasource=4,
                      iterations=1) {
  # SAFETY CHECK
  if (iterations < 1) {
    stop ('Iterations must be >=1')
  }
  if (!is.vector (names) && length (names) <= 1) {
    stop ('Names must be a vector of more than 1 character string')
  }
  # INIT
  tree$tip.label <- gsub ('_', ' ', tree$tip.label)
  names <- gsub ('_', ' ', names)
  # find matching and non-matching names
  matching.names <- names[names %in% tree$tip.label]
  nonmatching.names <- names[!names %in% tree$tip.label]
  # stop now if all names match names in tree or fuzzy is false
  if (!fuzzy || length (matching.names) == length (names)) {
    return (.mnEarlyReturn (tree, names, iterations))
  }
  # VARIABLES AND ENVIRONMENTS
  # place all parameters in an environment to save arguments
  paraenv <- new.env (parent=emptyenv ())
  paraenv$start.tree <- tree
  paraenv$grow.tree <- tree
  paraenv$datasource <- datasource
  paraenv$matching.names <- matching.names
  paraenv$names <- names
  # hold query name resolution results in a list
  qrylist <- .mnResolve (names=nonmatching.names, paraenv=paraenv)
  if (!nrow (qrylist$resolved) > 0) {  # stop now if none resolved
    return (.mnEarlyReturn (tree, names, iterations))
  }
  # hold subject name resolution results in a single env
  if (!exists ('sbjctenv')) {
    # allow user to specify sbjctenv outside of function
    # this allows the user to run the function several times
    # for different sets of names but with the same large tree
    # while reducing number of searches required
    sbjctenv <- new.env (parent=emptyenv ())
    .mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
  } else {
    paraenv$deja.vues <- sbjctenv$resolved$search.name
  }
  # add qrylist results to sbjctenv
  pull <- !as.vector (qrylist$resolved$search.name) %in%
    as.vector (sbjctenv$resolved$search.name)
  sbjctenv$resolved <- rbind (sbjctenv$resolved, qrylist$resolved[pull, ])
  sbjctenv$lineages <- c (sbjctenv$lineages, qrylist$lineages[pull])
  # hold all resulting trees and stats in a single env
  resenv <- new.env (parent=emptyenv ())
  resenv$trees <- list ()
  # ITERATE
  for (i in 1:iterations) {
    # loop through for each iteration, write results to resenv
    .mnMap (resenv=resenv, qrylist=qrylist, sbjctenv=sbjctenv,
            paraenv=paraenv)
  }
  # RETURN
  if (length (resenv$trees) > 1) {
    trees <- resenv$trees
    class (trees) <- 'multiPhylo'
  } else {
    trees <- resenv$trees[[1]]
  }
  return (trees)
}
# HIDDEN mapNames functions
.mnMap <- function (resenv, qrylist, sbjctenv, paraenv) {
  # INIT
  # randomise order to prevent order bias
  randomised <- sample (1:nrow (qrylist$resolved))
  qrylist$resolved <- qrylist$resolved[randomised, ]
  qrylist$lineages <- qrylist$lineages[randomised]
  paraenv$grow.tree <- paraenv$start.tree
  sbjctlist <- .mnTemporise(record=sbjctenv, tree=paraenv$grow.tree)  # convert to list
  # LOOP until qrylist in grow.tree or sbjct names exhausted
  while (TRUE) {
    if (nrow (sbjctlist$resolved) > 0) {
      not.in.tree <- which (!as.vector (qrylist$resolved$search.name) %in%
                              paraenv$grow.tree$tip.label)
      # loop through each in qrylist$resolved and match to sbjctlist$resolved
      for (i in not.in.tree) {
        # get lineage of the resovled qry name
        lineage <- qrylist$lineages[[i]]
        # match lineage to sbjct lineages
        lsrs <- .mnGetLSR (qry=lineage, sbjcts=sbjctlist$lineages)
        # of the most closely related lineages get the min dist -- calculated as
        # the smallest lsrs, this makes sure that they do not
        # form an in-group to the exclusion of the qry lineage
        possibles <- as.vector (which (lsrs == max (lsrs)))
        best.sbjcts <- sbjctlist$lineages[possibles]
        if (length (best.sbjcts) > 1) {
          min.lsr <- min (.mnGetLSR (qry=best.sbjcts[[1]],
                                     sbjcts=best.sbjcts[-1]))
        } else {
          # if there is only 1, it must be the best match in tree
          min.lsr <- max (lsrs)
        }
        if (max (lsrs) >= min.lsr) {
          paraenv$grow.tree <- .mnAddTip (tree=paraenv$grow.tree,
                                          tip.is=sbjctlist$resolved$tip.i[possibles],
                                          new.name=as.character (
                                            qrylist$resolved$search.name[i]))
          # update tip.i in sbjctlist
          sbjctlist <- .mnTemporise(record=sbjctlist, tree=paraenv$grow.tree)
        }
      }
    }
    # if all qry names in grow.tree or all names in tree sampled
    if (all (as.vector (qrylist$resolved$search.name) %in% paraenv$grow.tree$tip.label) |
          all (paraenv$start.tree$tip.label %in% paraenv$deja.vues)) {
      break
    }
    # if haven't broken out, update sbjctenv and extract new sbjctlist
    .mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
    sbjctlist <- .mnTemporise(record=sbjctenv, tree=paraenv$grow.tree)
  }
  # save results to resenv
  tree <- .mnExtract (tree=paraenv$grow.tree, names=paraenv$names)
  resenv$trees <- c (resenv$trees, list (tree))
}
.mnResolve <- function (names, paraenv) {
  # Resolve names using taxaResolve, return list of taxaResolve
  #  dataframe and a list of lineages for each name
  res <- list ()
  res['resolved'] <- list (taxaResolve (names=names,
                                        datasource=paraenv$datasource))
  # drop NAs
  res$resolved <- res$resolved[!is.na (res$resolved$name.string), ]
  # drop those w/o lineage
  res$resolved <- res$resolved[res$resolved$lineage != '', ]
  # separate lineages
  if (nrow (res$resolved) > 0) {
    res['lineages'] <- list (strsplit (as.vector (res$resolved$lineage),
                                       '\\|'))
  }
  return (res)
}
.mnGetLSR <- function (qry, sbjcts) {
  # get lowest shared ranks between qry lineage and subject lineages
  tds <- rep (NA, length (sbjcts))
  for (i in 1:length (sbjcts)) {
    tds[i] <- max (which (qry %in% sbjcts[[i]]))
  }
  tds
}
.mnSample <- function (paraenv) {
  # sample names from a tree in a way to reduce searching
  tree <- paraenv$grow.tree
  tip <- sample (paraenv$matching.names, 1)
  node <- tree$edge[tree$edge[ ,2] == which (tip == tree$tip.label), 1]
  while (TRUE) {
    children <- getChildren (tree, node=node)
    if (any (!children %in% paraenv$deja.vues)) {
      # exclude names already searched
      children <- children[!children %in% paraenv$deja.vues]
      return (children)
    }
    if (node == getSize (tree) + 1) {
      break
    }
    node <- getParent (tree, node=node)
  }
  vector ()
}
.mnResolveUpdate <- function (paraenv, sbjctenv) {
  # Update sbjenv -- by only searching when needed, reduce number of searches
  # get sample of names
  names <- .mnSample (paraenv)
  if (length (names) > 0) {
    # add names to deja.vues
    paraenv$deja.vues <- c (paraenv$deja.vues, names)
    # resolve these names
    res <- .mnResolve (names, paraenv)
    # stick to previous results
    if (!is.null (sbjctenv$resolved) & !is.null (sbjctenv$lineages)) {
      sbjctenv$resolved <- rbind (sbjctenv$resolved, res$resolved)
      sbjctenv$lineages <- c (sbjctenv$lineages, res$lineages)
    } else {
      sbjctenv$resolved <- res$resolved
      sbjctenv$lineages <- res$lineages
    }
  }
}
.mnTemporise <- function (record, tree) {
  # return a temporary record for name mapping
  # create a copy
  res <- list ()
  res$resolved <- record$resolved
  res$lineages <- record$lineages
  # add tip.i info
  tip.i <- match (res$resolved$search.name, tree$tip.label)
  not.in.tree <- is.na (tip.i)
  res$resolved$tip.i <- tip.i
  res$resolved <- res$resolved[!not.in.tree, ]
  res$lineages <- res$lineages[!not.in.tree]
  return (res)
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
  # return tree representing names, also stats
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