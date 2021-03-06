#' #' @name pinNames
#' #' @title Pin new tips to time-callibrated tree
#' #' @description Take a list of names with taxonomic lineages and age ranges, and pin them on
#' #' to existing time-callibrated tree.
#' #' @details This function ...
#' #' @template base_template
#' #' @export
#' #' @examples
#' #' # bring in the catarrhines data
#' 
#' pinNames <- function (tree, names, lineages, min.ages, max.ages,
#'                       fuzzy=TRUE, datasource=4,
#'                       iterations=1, resolve.list=NULL) {
#'   # SAFETY CHECK
#'   if (iterations < 1) {
#'     stop ('Iterations must be >=1')
#'   }
#'   if (!is.vector (names) && length (names) <= 1) {
#'     stop ('Names must be a vector of more than 1 character string')
#'   }
#'   if (!class (tree) %in%  c ('phylo', 'multiPhylo')) {
#'     stop ('Tree must be phylo or multiPhylo')
#'   }
#'   #TODO update addTip to handle node.label, for now drop node.label
#'   # INIT
#'   # drop underscores, check for branch lengths and drop node labels
#'   tree <- .mnClean (tree)
#'   names <- gsub ('_', ' ', names)
#'   # find matching and non-matching names
#'   tip.labels <- .mnGetNames (tree)
#'   matching.names <- names[names %in% tip.labels]
#'   nonmatching.names <- names[!names %in% tip.labels]
#'   # stop now if all names match names in tree or fuzzy is false
#'   if (!fuzzy || length (matching.names) == length (names)) {
#'     return (.mnEarlyReturn (tree, names, iterations))
#'   }
#'   # VARIABLES AND ENVIRONMENTS
#'   # place all parameters in an environment to save arguments
#'   paraenv <- new.env (parent=emptyenv ())
#'   if (class (tree) == 'multiPhylo') {
#'     paraenv$trees <- tree
#'     paraenv$start.tree <- tree[[sample (1:length (tree), 1)]]
#'   } else {
#'     paraenv$start.tree <- tree
#'   }
#'   paraenv$grow.tree <- paraenv$start.tree
#'   paraenv$datasource <- datasource
#'   paraenv$matching.names <- matching.names
#'   paraenv$names <- names
#'   paraenv$by.age <- TRUE
#'   # get query and subject resolved names
#'   qrylist <- list ()
#'   qrylist$lineages <- lineages
#'   qrylist$resolved <- data.frame (search.name=names,
#'                                   min.age=min.ages,
#'                                   max.age=max.ages)
#'   if (!is.null (resolve.list)) {
#'     # if resolve.list provided, no need to do any searches
#'     # unpack resolved names into sbjctenv
#'     sbjctenv <- new.env (parent=emptyenv ())
#'     sbjctenv$resolved <- resolve.list$resolved
#'     sbjctenv$lineages <- resolve.list$lineages
#'     paraenv$deja.vues <- tree$tip.label
#'   } else {
#'     # hold subject name resolution results in a single env
#'     sbjctenv <- new.env (parent=emptyenv ())
#'     .mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
#'   }
#'   # hold all resulting trees and stats in a single env
#'   resenv <- new.env (parent=emptyenv ())
#'   resenv$trees <- list ()
#'   # ITERATE
#'   # loop through for each iteration, write results to resenv
#'   m_ply (.data=data.frame(iteration=1:iterations), .fun=.mnMap,
#'          resenv=resenv, qrylist=qrylist, sbjctenv=sbjctenv,
#'          paraenv=paraenv)
#'   # RETURN
#'   if (length (resenv$trees) > 1) {
#'     trees <- resenv$trees
#'     class (trees) <- 'multiPhylo'
#'   } else {
#'     trees <- resenv$trees[[1]]
#'   }
#'   return (trees)
#' }
#' 
#' #' @name mapNames
#' #' @title Map names to a tree
#' #' @description Map names to a tree with the option of searching online
#' #' to perform fuzzy-name matching.
#' #' @details This function firsts matches names to tip labels in tree, if not all names
#' #' are present in the tree and 'fuzzy' is true, it will then search the online taxonomic
#' #' names resolution service Global Names Resolver. Taxonomic IDs are then matched
#' #' between tip labels and names, if names are still not mapped to tree, names are mapped
#' #' to tips using random mapping. For example, if an unmapped name is a member of the same
#' #' genus as several species in the tree, the name will be mapped to the tree with any one
#' #' of these at random. In cases where an unmapped name can only be added to the tree by adding
#' #' an extra tip, then a new tip is added at random to the pendant edge of a suitable tip at
#' #' a random time point (it assumes the tree is time callibrated) -- see examples.
#' #' 
#' #' In cases where a large proportion of the tips are mapped using random placement, it may
#' #' be best to generate a distribution of trees that represent the range of possible phylogenetic
#' #' relationships for the names and tree given. If iterations is greater than 1, instead of a
#' #' single tree, a multiPhylo object of length iterations is returned representing this
#' #' range of possible mappings.
#' #' 
#' #' For non-time-callibrated or trees without branch lengths, the function using taxonomic
#' #' distance for random placement and returns a tree without branch lengths.
#' #' 
#' #' Does not require trees to be bifurcating, but will randomly resolve polytomies for
#' #' every iteration if not.
#' #' 
#' #' @template base_template
#' #' @param names vector of names of tips to be extracted
#' #' @param fuzzy boolean, if true will search Global Names Resolver online
#' #' @param datasource GNR datasource ID, default 4 for NCBI
#' #' @param iterations how many times to repeat?
#' #' @param resolve.list pre-resolved names generated with mapNamesPreDownload,
#' #' default NULL
#' #' @export
#' #' @examples
#' #' # bring in the catarrhines data
#' #' data ('catarrhines')
#' #' # we want to map these names to the catarrhines tree
#' #' names <- c ('Homo sapiens', 'Pongo pygmaeus', 'Gorilla gorila', 'Pan troglodytes',
#' #'             'Homo erectus', 'Homo neanderthalensis', 'Hylobates')
#' #' # 4 of the names are already in the tree (one has a spelling mistake) but the other 3
#' #' # do not exist and will be mapped randomly based on resolved taxonomic lineages of the
#' #' # names already in the tree and the names to be added.
#' #' hominid.map <- mapNames (tree=catarrhines, names=names, fuzzy=TRUE)
#' #' plot (hominid.map)
#' 
#' mapNames <- function (tree, names, fuzzy=TRUE, datasource=4,
#'                       iterations=1, resolve.list=NULL) {
#'   # SAFETY CHECK
#'   if (iterations < 1) {
#'     stop ('Iterations must be >=1')
#'   }
#'   if (!is.vector (names) && length (names) <= 1) {
#'     stop ('Names must be a vector of more than 1 character string')
#'   }
#'   if (!class (tree) %in%  c ('phylo', 'multiPhylo')) {
#'     stop ('Tree must be phylo or multiPhylo')
#'   }
#'   #TODO update addTip to handle node.label, for now drop node.label
#'   # INIT
#'   # drop underscores, check for branch lengths and drop node labels
#'   tree <- .mnClean (tree)
#'   names <- gsub ('_', ' ', names)
#'   # find matching and non-matching names
#'   tip.labels <- .mnGetNames (tree)
#'   matching.names <- names[names %in% tip.labels]
#'   nonmatching.names <- names[!names %in% tip.labels]
#'   # stop now if all names match names in tree or fuzzy is false
#'   if (!fuzzy || length (matching.names) == length (names)) {
#'     return (.mnEarlyReturn (tree, names, iterations))
#'   }
#'   # VARIABLES AND ENVIRONMENTS
#'   # place all parameters in an environment to save arguments
#'   paraenv <- new.env (parent=emptyenv ())
#'   if (class (tree) == 'multiPhylo') {
#'     paraenv$trees <- tree
#'     paraenv$start.tree <- tree[[sample (1:length (tree), 1)]]
#'   } else {
#'     paraenv$start.tree <- tree
#'   }
#'   paraenv$grow.tree <- paraenv$start.tree
#'   paraenv$datasource <- datasource
#'   paraenv$matching.names <- matching.names
#'   paraenv$names <- names
#'   paraenv$by.age <- FALSE
#'   # get query and subject resolved names
#'   if (!is.null (resolve.list)) {
#'     # if resolve.list provided, no need to do any searches
#'     # unpack resolved names into qrylist and sbjctenv
#'     pull <- resolve.list$resolved$search.name %in% nonmatching.names
#'     if (sum (pull) < 1) {
#'       return (.mnEarlyReturn (tree, names, iterations))
#'     }
#'     qrylist <- list ()
#'     qrylist$resolved <- resolve.list$resolved[pull, ]
#'     qrylist$lineages <- resolve.list$lineages[pull]
#'     sbjctenv <- new.env (parent=emptyenv ())
#'     sbjctenv$resolved <- resolve.list$resolved
#'     sbjctenv$lineages <- resolve.list$lineages
#'     paraenv$deja.vues <- tree$tip.label
#'   } else {
#'     # hold query name resolution results in a list
#'     qrylist <- .mnResolve (names=nonmatching.names, paraenv=paraenv)
#'     if (!nrow (qrylist$resolved) > 0) {  # stop now if none resolved
#'       return (.mnEarlyReturn (tree, names, iterations))
#'     }
#'     # hold subject name resolution results in a single env
#'     sbjctenv <- new.env (parent=emptyenv ())
#'     .mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
#'     # add qrylist results to sbjctenv
#'     pull <- !as.vector (qrylist$resolved$search.name) %in%
#'       as.vector (sbjctenv$resolved$search.name)
#'     sbjctenv$resolved <- rbind (sbjctenv$resolved, qrylist$resolved[pull, ])
#'     sbjctenv$lineages <- c (sbjctenv$lineages, qrylist$lineages[pull])
#'   }
#'   # hold all resulting trees and stats in a single env
#'   resenv <- new.env (parent=emptyenv ())
#'   resenv$trees <- list ()
#'   # ITERATE
#'   # loop through for each iteration, write results to resenv
#'   m_ply (.data=data.frame(iteration=1:iterations), .fun=.mnMap,
#'          resenv=resenv, qrylist=qrylist, sbjctenv=sbjctenv,
#'          paraenv=paraenv)
#'   # RETURN
#'   if (length (resenv$trees) > 1) {
#'     trees <- resenv$trees
#'     class (trees) <- 'multiPhylo'
#'   } else {
#'     trees <- resenv$trees[[1]]
#'   }
#'   return (trees)
#' }
#' # HIDDEN mapNames functions
#' .mnMap <- function (iteration, resenv, qrylist, sbjctenv, paraenv) {
#'   # INIT
#'   # randomise order to prevent order bias
#'   randomised <- sample (1:nrow (qrylist$resolved))
#'   qrylist$resolved <- qrylist$resolved[randomised, ]
#'   qrylist$lineages <- qrylist$lineages[randomised]
#'   if (!is.binary.tree (paraenv$start.tree)) {
#'     paraenv$grow.tree <- multi2di (paraenv$start.tree)
#'   } else {
#'     paraenv$grow.tree <- paraenv$start.tree
#'   }
#'   sbjctlist <- .mnTemporise(record=sbjctenv, tree=paraenv$grow.tree)  # convert to list
#'   # LOOP until qrylist in grow.tree or sbjct names exhausted
#'   while (TRUE) {
#'     if (nrow (sbjctlist$resolved) > 0) {
#'       not.in.tree <- which (!as.vector (qrylist$resolved$search.name) %in%
#'                               paraenv$grow.tree$tip.label)
#'       # loop through each in qrylist$resolved and match to sbjctlist$resolved
#'       for (i in not.in.tree) {
#'         # get lineage of the resovled qry name
#'         lineage <- qrylist$lineages[[i]]
#'         # match lineage to sbjct lineages
#'         lsrs <- .mnGetLSR (qry=lineage, sbjcts=sbjctlist$lineages)
#'         # of the most closely related lineages get the min dist -- calculated as
#'         # the smallest lsrs, this makes sure that they do not
#'         # form an in-group to the exclusion of the qry lineage
#'         possibles <- as.vector (which (lsrs == max (lsrs)))
#'         best.sbjcts <- sbjctlist$lineages[possibles]
#'         if (length (best.sbjcts) < length (sbjctlist$lineages)) {
#'           # if the length of best is the length of all, then no resolution withing sbjcts
#'           if (paraenv$by.age) {
#'             paraenv$grow.tree <- .mnAddTipWAge (tree=paraenv$grow.tree,
#'                                                 tip.is=sbjctlist$resolved$tip.i[possibles],
#'                                                 new.name=as.character (
#'                                                   qrylist$resolved$search.name[i]),
#'                                                 max.age=qrylist$resolved$max.age[i],
#'                                                 min.age=qrylist$resolved$min.age[i])
#'           } else {
#'               paraenv$grow.tree <- .mnAddTip (tree=paraenv$grow.tree,
#'                                               tip.is=sbjctlist$resolved$tip.i[possibles],
#'                                               new.name=as.character (
#'                                                 qrylist$resolved$search.name[i]))
#'           }
#'           # update tip.i in sbjctlist
#'           sbjctlist <- .mnTemporise(record=sbjctlist, tree=paraenv$grow.tree)
#'         }
#'       }
#'     }
#'     # if all qry names in grow.tree or all names in tree sampled
#'     if (all (as.vector (qrylist$resolved$search.name) %in% paraenv$grow.tree$tip.label) |
#'           all (paraenv$start.tree$tip.label %in% paraenv$deja.vues)) {
#'       break
#'     }
#'     # if haven't broken out, update sbjctenv and extract new sbjctlist
#'     .mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
#'     sbjctlist <- .mnTemporise(record=sbjctenv, tree=paraenv$grow.tree)
#'   }
#'   # save results to resenv
#'   if (paraenv$by.age) {
#'     resenv$trees <- c (resenv$trees, list (paraenv$grow.tree))
#'   } else {
#'     tree <- .mnExtract (tree=paraenv$grow.tree, names=paraenv$names)
#'     resenv$trees <- c (resenv$trees, list (tree))
#'   }
#'   # choose new random tree
#'   if (!is.null (paraenv$trees)) {
#'     paraenv$start.tree <- paraenv$trees[[sample (1:length (paraenv$trees), 1)]]
#'   }
#' }
#' .mnResolve <- function (names, paraenv) {
#'   # Resolve names using taxaResolve, return list of taxaResolve
#'   #  dataframe and a list of lineages for each name
#'   res <- list ()
#'   res['resolved'] <- list (taxaResolve (names=names,
#'                                         datasource=paraenv$datasource))
#'   # drop NAs
#'   res$resolved <- res$resolved[!is.na (res$resolved$name.string), ]
#'   # drop those w/o lineage
#'   res$resolved <- res$resolved[res$resolved$lineage != '', ]
#'   # separate lineages
#'   if (nrow (res$resolved) > 0) {
#'     res['lineages'] <- list (strsplit (as.vector (res$resolved$lineage),
#'                                        '\\|'))
#'   }
#'   return (res)
#' }
#' .mnGetLSR <- function (qry, sbjcts) {
#'   # get lowest shared ranks between qry lineage and subject lineages
#'   tds <- rep (NA, length (sbjcts))
#'   for (i in 1:length (sbjcts)) {
#'     tds[i] <- max (which (qry %in% sbjcts[[i]]))
#'   }
#'   tds
#' }
#' .mnSample <- function (paraenv) {
#'   # sample names from a tree in a way to reduce searching
#'   tree <- paraenv$grow.tree
#'   if (length (paraenv$matching.names) > 1) {
#'     tip <- sample (paraenv$matching.names, 1)
#'   } else {
#'     tip <- sample (tree$tip.label, 1)
#'   }
#'   node <- tree$edge[tree$edge[ ,2] == which (tip == tree$tip.label), 1]
#'   while (TRUE) {
#'     children <- getChildren (tree, node=node)
#'     if (any (!children %in% paraenv$deja.vues)) {
#'       # exclude names already searched
#'       children <- children[!children %in% paraenv$deja.vues]
#'       return (children)
#'     }
#'     if (node == getSize (tree) + 1) {
#'       break
#'     }
#'     node <- getParent (tree, node=node)
#'   }
#'   vector ()
#' }
#' .mnResolveUpdate <- function (paraenv, sbjctenv) {
#'   # Update sbjenv -- by only searching when needed, reduce number of searches
#'   # get sample of names
#'   names <- .mnSample (paraenv)
#'   if (length (names) > 0) {
#'     # add names to deja.vues
#'     paraenv$deja.vues <- c (paraenv$deja.vues, names)
#'     # resolve these names
#'     res <- .mnResolve (names, paraenv)
#'     # stick to previous results
#'     if (!is.null (sbjctenv$resolved) & !is.null (sbjctenv$lineages)) {
#'       sbjctenv$resolved <- rbind (sbjctenv$resolved, res$resolved)
#'       sbjctenv$lineages <- c (sbjctenv$lineages, res$lineages)
#'     } else {
#'       sbjctenv$resolved <- res$resolved
#'       sbjctenv$lineages <- res$lineages
#'     }
#'   }
#' }
#' .mnTemporise <- function (record, tree) {
#'   # return a temporary record for name mapping
#'   # create a copy
#'   res <- list ()
#'   res$resolved <- record$resolved
#'   res$lineages <- record$lineages
#'   # add tip.i info
#'   tip.i <- match (res$resolved$search.name, tree$tip.label)
#'   not.in.tree <- is.na (tip.i)
#'   res$resolved$tip.i <- tip.i
#'   res$resolved <- res$resolved[!not.in.tree, ]
#'   res$lineages <- res$lineages[!not.in.tree]
#'   return (res)
#' }
#' .mnAddTipWAge <- function (tree, tip.is, new.name,
#'                            min.age, max.age) {
#'   print(new.name)
#'   # add new tip
#'   if (length (tip.is) == 1) {
#'     edge <- which (tip.is == tree$edge[ ,2])
#'     edge.age <- getAge (tree, edge=edge)
#'     pull <- edge.age$max.age >= max.age
#'     if (!pull) {
#'       edge <- NA
#'     }
#'   } else {
#'     # randomly map new tip to any edge in the clade
#'     # represented by the matching tip.is
#'     # find the parent node of all the matching tips
#'     children <- tree$tip.label[tip.is]
#'     parent.node <- getParent (tree, tips=children)
#'     # get all descending
#'     edges <- getEdges(tree, node=parent.node)
#'     # INCLUDING supporting edge!!!
#'     edges <- c (edges, which (tree$edge[ ,2] ==parent.node))
#'     # check ages
#'     edge.ages <- getAge (tree, edge=edges)
#'     pull <- edge.ages$max.age >= max.age
#'     edges <- edges[pull]
#'     # choose edge at random
#'     if (length (edges) > 1) {
#'       edge <- sample (edges, size=1)
#'     } else {
#'       edge <- edges[1]
#'     }
#'   }
#'   # exit if no suitable edge found
#'   if (!is.numeric (edge)) {
#'     return (tree)
#'   }
#'   # random tip and node age
#'   age.range <- getAge (tree, edge=edge)
#'   node.age <- runif (n=1, min=age.range[1, 'max.age'],
#'                      max=age.range[1, 'max.age'])
#'   tip.age <- runif (n=1, min=min.age, max=node.age)
#'   tree <- addTip (tree, edge=edge, tip.name=new.name, node.age=node.age,
#'                   tip.age=tip.age)
#'   return (tree)
#' }
#' .mnAddTip <- function (tree, tip.is, new.name) {
#'   # add new tip
#'   if (length (tip.is) == 1) {
#'     # add to the pendant edge at a random time point
#'     edge <- which (tip.is == tree$edge[ ,2])
#'   } else {
#'     # randomly map new tip to any edge in the clade
#'     # represented by the matching tip.is
#'     # find the parent node of all the matching tips
#'     children <- tree$tip.label[tip.is]
#'     parent.node <- getParent (tree, tips=children)
#'     # get descending and supporting edges
#'     edges <- which (tree$edge[ ,1] == parent.node)
#'     edges <- c (edges, which (tree$edge[ ,2] == parent.node))
#'     # choose edge at random, bias sample based on branch length
#'     edge <- sample (edges, size=1, prob=tree$edge.length[edges])
#'   }
#'   # random node age somewhere on edge
#'   age.range <- getAge (tree, edge=edge)
#'   node.age <- runif (n=1, min=age.range[1, 'min.age'],
#'                      max=age.range[1, 'max.age'])
#'   # random tip age if not ultrametric
#'   if (is.ultrametric (tree)) {
#'     tip.age <- 0
#'   } else {
#'     # use random tip age in shared clade
#'     possibles <- getAge (tree, node=tip.is)$age
#'     possibles <- possibles[possibles < node.age]
#'     tip.age <- runif (n=1, min=min(possibles), max=max(possibles))
#'   }
#'   tree <- addTip (tree, edge=edge, tip.name=new.name, node.age=node.age,
#'                   tip.age=tip.age)
#'   return (tree)
#' }
#' .mnExtract <- function (tree, names) {
#'   # return tree representing names
#'   if (class (tree) == 'multiPhylo') {
#'     tree <- tree[[sample (1:length (tree), 1)]]
#'   }
#'   if (sum (names %in% tree$tip.label) > 1) {
#'     res.tree <- drop.tip (tree, tip = tree$tip.label[!tree$tip.label %in% names])
#'     return (res.tree)
#'   } else {
#'     warning ('Too few names could be mapped to tree')
#'     return (NA)
#'   }
#' }
#' .mnEarlyReturn <- function (tree, names, iterations) {
#'   # if mapNames returns tree early, return an object
#'   #  that was expect i.e. phylo or multiphylo
#'   tree <- .mnExtract (tree, names)
#'   if (iterations == 1) {
#'     return (tree)
#'   } else {
#'     trees <- list ()
#'     for (i in 1:iterations) {
#'       trees <- c (trees, list (tree))
#'     }
#'     class (trees) <- 'multiPhylo'
#'     return (trees)
#'   }
#' }
#' .mnClean <- function (trees) {
#'   # drop _ in names of trees in a multiphylo
#'   # remove node labels
#'   # set branch lengths to 1 if no branch lengths
#'   if (class (trees) == 'multiPhylo') {
#'     .drop <- function (i) {
#'       tree <- trees[[i]]
#'       tree$tip.label <- gsub ('_', ' ', tree$tip.label)
#'       tree$node.label <- NULL
#'       if (is.null (tree$edge.length)) {
#'         tree$edge.length <- rep (1, nrow (tree$edge))
#'       }
#'       res <<- c (res, list (tree))
#'     }
#'     res <- list ()
#'     m_ply (.data=data.frame (i=1:length (trees)), .fun=.drop)
#'     class (res) <- 'multiPhylo'
#'     return (res)
#'   } else {
#'     trees$tip.label <- gsub ('_', ' ', trees$tip.label)
#'     trees$node.label <- NULL
#'     if (is.null (trees$edge.length)) {
#'       trees$edge.length <- rep (1, nrow (trees$edge))
#'     }
#'     return (trees)
#'   }
#' }
#' .mnGetNames <- function(trees) {
#'   # get tip names from a multiphylo
#'   if (class (trees) == 'multiPhylo') {
#'     .get <- function (i) {
#'       res <<- c (res, trees[[i]]$tip.label)
#'     }
#'     res <- NULL
#'     m_ply (.data=data.frame (i=1:length (trees)), .fun=.get)
#'     res <- unique (res)
#'   } else {
#'     res <- trees$tip.label
#'   }
#'   return (res)
#' }
#' 
#' #' @name mapNamesPreDownload
#' #' @title Download names for mapNames()
#' #' @description Run this function before running mapNames() when
#' #' mapping multiple sets of names.
#' #' @details mapNames() will try to minimise the number of searches it
#' #' needs to do to map names given to a phylogenetic tree. If you're running
#' #' mapNames() multiple times for different vectors of names, it may be faster
#' #' to resolve names via this function and then pass these to mapNames()
#' #' 
#' #' @template base_template
#' #' @param names vector of names to be resolved and all names in tree to be mapped
#' #' @param datasource GNR datasource ID, default 4 for NCBI
#' #' @examples
#' #' # bring in the catarrhines data
#' #' data ('catarrhines')
#' #' names1 <- c ('Homo sapiens', 'Pongo pygmaeus', 'Gorilla gorila', 'Pan troglodytes',
#' #'              'Homo erectus', 'Homo neanderthalensis', 'Hylobates')
#' #' names2 <- c ('Hylobates', 'Homo', 'Macca')
#' #' resolve.list <- mapNamesPreDownload (names=c (names1, names2, catarrhines$tip.label),
#' #'                                      datasource=4)
#' #' map1 <- mapNames (tree=catarrhines, names=names1, resolve.list=resolve.list)
#' #' map2 <- mapNames (tree=catarrhines, names=names2, resolve.list=resolve.list)
#' #' plot (map1)
#' #' plot (map2)
#' 
#' mapNamesPreDownload <- function (names, datasource=4) {
#'   paraenv <- list ('datasource'=datasource)
#'   # parse names
#'   names <- unique (names)
#'   names <- gsub ('_', ' ', names)
#'   # search names
#'   resolved.list <- .mnResolve (names, paraenv)
#'   return (resolved.list)
#' }