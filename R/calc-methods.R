#' @name reduceTree
#' @title Reduce tree by taxonomic rank
#' @description Reduce tips by identifying tips in the same taxonomic
#' group at a specifed taxonomic rank
#' @details Search Global Names Resolver to find shared taxonomic groups of
#' all tips. Reduce tree by condensing members of the same taxonomic group
#' into a single tip. The user can specify at what rank the tree should be reduced.
#' N.B. unresolved names are dropped from the tree and are not represented in tip counts.
#' @template base_template
#' @param level rank by which tips will be reduce (e.g. kingdom)
#' @param datasource a number indicating the GNR datasource to search against
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

reduceTree <- function (tree, level, datasource = 4) {
  .match <- function (level, rank) {
    # a series of patterns to find matching level
    #  to control for subs, supers etc...
    matched <- grep (level, rank)
    pattern <- paste0 ('(super|sub|parv|infra)', level)
    avoid <- grep (pattern, rank)
    matched <- matched[!matched %in% avoid]
    if (length (matched) > 1) {
      # if still greater than 1, choose at random
      cat ('Random! In .match\n')
      matched <- sample (matched, 1)
    }
    matched
  }
  .pull <- function (lineage, rank) {
    # extract matching level name
    lineage <- strsplit (as.character (lineage), split = '\\|')[[1]]
    rank <- strsplit (as.character (rank), split = '\\|')[[1]]
    resolved.name <- lineage[.match (level, rank)]
    if (length (resolved.name) > 0) {
      return (resolved.name)
    } else {
      NA
    }
  }
  # resolve names
  names <- gsub ("_", " ", tree$tip.label)
  gnr.obj <- taxaResolve (names, datasource = datasource)
  data <- data.frame (lineage = gnr.obj['lineage'], rank = gnr.obj['rank'])
  # match level to rank, and extract lineage
  higher.names <- mdply (.data = data, .fun = .pull)[ ,3]
  # drop all unresolved tips
  if (any (is.na (higher.names))) {
    cat (paste0 ("[", sum (is.na (higher.names)),
                 "] unresovled names have been dropped.\n"))
    tree <- drop.tip (tree, tree$tip.label[is.na (higher.names)])
    higher.names <- higher.names[!is.na (higher.names)]
    if (length (higher.names) < 2) {
      stop ('Too few tips in tree -- try a different level.')
    }
  }
  # record new names in names.df
  higher.names.table <- table (higher.names)
  tip.labels <- paste0 (names (higher.names.table), " -- ", higher.names.table,
                        " tip(s)")
  tip.labels <- tip.labels[match (higher.names, names (higher.names.table))]
  name.df <- data.frame (names = tree$tip.label, higher.names, tip.labels)
  # drop all non-unique tips
  to.drop <- -1 * match (unique (tip.labels), tip.labels)
  tree <- drop.tip (tree, tree$tip.label[to.drop])
  # rename based on names.df
  bool <- name.df[['names']] %in% tree$tip.label
  tree$tip.label <- as.character(name.df[['tip.labels']][bool])
  tree
}

#' @name calcComPhyMets
#' @title Calculate Community Phylogenetic Metrics
#' @description One-stop function for calculating a range of community phylogenetic metrics
#' from tree and community matrix.
#' @details Metrics calculated: PD (types 1, 2 and 3), PSV, PSD, PSE, PSC, and PSR.
#' The PD metric in this function has 3 different types:
#'  PD1: considers the case if all other taxa from the phylogeny are dropped. In this
#'  way, type 1 is context independent but has a minimum number of taxa of 2. (This is
#'  the normal PD.)
#'  PD2: sums the lengths of all the edges connecting the specified taxa to the terminal node.
#'  PD3: sums the lengths of all the edges represented uniquely by the specified taxa.
#' @template base_template
#' @param cmatrix community/trait data matrix (cols taxa, rows sites)
#' @param metric what metric to use as vector, default all metrics
#' @param min.spp minimum number of species at site, default 2
#' @references Faith, D. (1992). Conservation evaluation and phylogenetic diversity.
#'  Biological Conservation, 61, 1–10.
#'
#' Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007) Phylogenetic measures of
#'  biodiversity. American Naturalist, 169, E68-E83
#'  @examples
#'  data (catarrhines)
#'  cmatrix <- randCommData (catarrhines, nsites=20, nspp=10)
#'  res <- calcComPhyMets (cmatrix, catarrhines)
#'  print (res)

calcComPhyMets <- function(cmatrix, tree, min.spp = 2,
                           metrics = c ('PD1', 'PD2', 'PD3', 'PSV',
                                     'PSD', 'PSE', 'PSC')) {
  calcPD <- function (i, type) {
    # get tips present in this site
    tips <- colnames (cmatrix)[cmatrix[i, ] > 0]
    if (length (tips) <= min.spp) {
      return (NA)
    }
    # get edges given type
    edges <- getEdges (tree, tips=tips, type=type)
    sum (tree$edge.length[edges])
  }
  if ('PD1' %in% metrics & min.spp < 2) {
    stop("Cannot compute type 1 PD for fewer than 2 species")
  }
  # ensure tip names and colnames are shared
  if (!any (colnames (cmatrix) %in% tree$tip.label)) {
    stop ('No tip labels in community matrix')
  }
  # convert to matrix, if not
  if (!is.matrix (cmatrix)) {
    cmatrix.dimnames <- dimnames (cmatrix)
    cmatrix <- as.matrix (cmatrix)
    dimnames (cmatrix) <- cmatrix.dimnames
  }
  # add site names if none
  if (is.null (rownames (cmatrix))) {
    rownames (cmatrix) <- paste0 ('site_', 1:nrow (cmatrix))
  }
  # init res data.frame
  res <- data.frame ('Site'=rownames (cmatrix))
  # calc each metric and add to res
  if ('PD1' %in% metrics) {
    temp <- mdply (.data = data.frame (i = 1:nrow (cmatrix)),
                   .fun = calcPD, type=1)[ ,2]
    res <- data.frame (res, 'PD1'=temp)
  }
  if ('PD2' %in% metrics) {
    temp <- mdply (.data = data.frame (i = 1:nrow (cmatrix)),
                   .fun = calcPD, type=2)[ ,2]
    res <- data.frame (res, 'PD2'=temp)
  }
  if ('PD3' %in% metrics) {
    temp <- mdply (.data = data.frame (i = 1:nrow (cmatrix)),
                   .fun = calcPD, type=3)[ ,2]
    res <- data.frame (res, 'PD3'=temp)
  }
  if ('PSV' %in% metrics) {
    temp <- psv (samp=cmatrix, tree=tree)[ ,1]
    res <- data.frame (res, 'PSV'=temp)
  }
  if ('PSR' %in% metrics) {
    temp <- psr (samp=cmatrix, tree=tree)[ ,1]
    res <- data.frame (res, 'PSR'=temp)
  }
  if ('PSE' %in% metrics) {
    temp <- pse (samp=cmatrix, tree=tree)[ ,1]
    res <- data.frame (res, 'PSE'=temp)
  }
  if ('PSC' %in% metrics) {
    temp <- psc (samp=cmatrix, tree=tree)[ ,1]
    res <- data.frame (res, 'PSC'=temp)
  }
  if ('PSD' %in% metrics) {
    temp <- psd (samp=cmatrix, tree=tree)[ ,1]
    res <- data.frame (res, 'PSD'=temp)
  }
  res
}

#' @name calcED
#' @title Calculate Evolutionary Distinctiveness
#' @description Calculate species evolutionary distinctiveness (ED) using one of three methods:
#'  Fair Proportion (FP), Equal Splits (ES) or Pendant Edge (PE).
#' @details Evolutionary distinctiveness is a measure of how much independent evolution a species
#'  represents. Multiple methods exist for its calculation all of which require an ultrametric
#'  phylogenetic tree. The methods used here are Pendant Edge (PE) the length of a species branch
#'  that connects it to the tree (Altschul and Lipman, 1990), Fair Proportion (FP) the total
#'  proportion of the phylogenetic tree that a species represents where each branch is equally
#'  divided between all descdendants (Isaac et al. 2007) and Equal Splits (ES) where branch lengths
#'  are equally divided between descendents at every node in the tree (Redding and Mooers 2006)
#'  
#'  N.B. \code{picante} already has a function for doing this. \code{calcED}, however, uses \code{plyr}
#'  vectorisation and is much faster.
#' @template base_template
#' @param tips vector of tips for which ED is calculated, else 'all'
#' @param type method of ED calculation, either 'all', 'FP', 'ES' or 'PE' (default FP)
#' @references Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). 
#'  Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#'  
#'  Altschul, S., and Lipman, D. (1990). Equal animals. Nature.
#'  
#'  Redding, D. W., and Mooers, A. Ø. (2006). Incorporating evolutionary measures into conservation 
#'  prioritization. Conservation Biology : The Journal of the Society for Conservation Biology, 
#'  20(6), 1670–8.
#' @export
#' @examples
#' #

calcED <- function (tree, tips = 'all', type = 'FP') {
  calcFairProportion <- function (tips) {
    # first create an edge.matrix that will collect proportions of branch
    #  per species
    .calc <- function (node) {
      # find all of node's children
      children <- getChildren (tree, node)
      # work out proportion of branch between all children
      edge.length <- tree$edge.length[which (tree$edge[ ,2] == node)]
      edge.length.share <- edge.length/length (children)
      # fill in edge.matrix
      edge.matrix[node, children] <<- edge.length.share
    }
    edge.matrix <- matrix (0, ncol = getSize (tree),
                           nrow = tree$Nnode + getSize (tree))
    colnames (edge.matrix) <- tree$tip.label
    nodes <- 1:(getSize (tree) + tree$Nnode)
    # remove root node
    nodes <- nodes[-1 * (getSize (tree) + 1)]
    m_ply (.data = data.frame (node = nodes), .fun = .calc)
    res <- data.frame (FP = colSums (edge.matrix[ ,tips]))
  }
  calcEqualSplits <- function (tips) {
    .calc <- function (sp) {
      # find all edges that connect sp to root
      edges <- getEdges (tree, tips = as.character (sp), type = 2)
      # get all lengths and divide by the number of splits, from that node
      res <- 0
      for (i in length (edges):1) {
        res <- res + tree$edge.length[edges[i]]/i
      }
      res
    }
    res <- mdply (.data = data.frame (sp = tips), .fun = .calc)
    res <- data.frame (ES = res[, 2])
    rownames (res) <- tips
    res
  }
  calcPendantEdge <- function (tips) {
    getEdgeLength <- function (sp) {
      tip.index <- which (tree$tip.label == sp)
      tip.edge <- which (tree$edge[, 2] == tip.index)
      tree$edge.length[tip.edge]
    }
    res <- mdply (.data = data.frame (sp = tips), .fun = getEdgeLength)
    res <- data.frame (PE = res[, 2])
    rownames (res) <- tips
    res
  }
  if (tips[1] == 'all') {
    tips <- tree$tip.label
  }
  if (type == 'all') {
    FPs <- calcFairProportion (tips)
    ESs <- calcEqualSplits (tips)
    PEs <- calcPendantEdge (tips)
    return (cbind (FPs, ESs, PEs))
  }
  if (type == 'FP') {
    EDs <- calcFairProportion (tips)
  } else if (type == 'ES') {
    EDs <- calcEqualSplits (tips)
  } else if (type == 'PE') {
    EDs <- calcPendantEdge (tips)
  } else {
    stop (paste0 ('Type [', type, '] not known.
                  Type must be either all, FP, ES or PE'))
  }
  colnames (EDs) <- 'ED'
  EDs
}

#' @name randCommData
#' @title Generate Random Community Data Based on Community Phylogeny
#' @description Generate a community matrix based on a phylogeny without phylogenetic signal.
#' @details The default is to generate incidences, but if \code{lam} is numeric will
#' generate using lam as lambda.
#' @template base_template
#' @param nsites the number of sites
#' @param the total number of species for each site
#' @param lambda for Poisson distribution (default NULL for incidences)

randCommData <- function (tree, nsites, nspp, lam = NULL) {
  ## TODO: write test
  ntips <- length (tree$tip.label)
  output <- matrix (rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames (output) <- tree$tip.label
  if (is.null (lam)) {
    for (i in 1:nsites) {
      output[i, ] <- sample (c (rep (1, nspp), rep (0, ntips - nspp)))
    }
  } else {
    for (i in 1:nsites) {
      output[i, ] <- sample (rpois (ntips, lam))
    }
  }
  output
}

#' @name genCommData
#' @title Generate Clustered/Overdispersed Data Based on Community Phylogeny
#' @description Generate a community matrix based on a phylogeny with phylogenetic signal.
#' @details Abundances/incidences for a commuity are generated clustering around a focal
#'  point in the phylogeny (clustering \code{psi} > 0) or diserpsing around a focal point
#'  (overdispersion \code{psi} < 0). Return either incidence or abundance data.
#' @template base_template
#' @param focal numeric index, indicating which tip to perform cluster/dispersion
#' @param psi scaling coefficient: larger the number the more pronounced the
#' effect by a power law (hence 0 not allowed)
#' @param mean.incid the mean incidence of species in the community
#' @param mean.abun the mean abundance per site, if given output will be abundances
#' @param nsites number of sites

genCommData <- function(tree, focal, psi = 1, mean.incid, mean.abun = FALSE,
                        nsites = 1) {
  ##TODO: write test
  invertVector <- function(dists) {
    # for reversing the probs for overdispersion
    u <- sort(unique(dists), TRUE)
    s <- sort(u)
    probs <- rep(NA, length(dists))
    for (i in 1:length(u)) {
      probs[u[i] == dists] <- s[i]
    }
    return (probs)
  }
  genAbuns <- function(row) {
    # for generating abundances for each row
    out.row <- rep(0, ntips)
    temp.probs <- probs
    temp.probs[row < 1] <- 0
    abundance <- abs(ceiling(rnorm(1, mean = mean.abun)))
    if (abundance == 0) {
      return (out.row)
    } else {
      abuns <- sample(1:ntips, abundance, prob = temp.probs, replace = TRUE)
      abuns <- table(abuns)
      out.row[as.numeric(names(abuns))] <- abuns
      return (out.row)
    }
  }
  ntips <- length(tree$tip.label)
  output <- matrix(rep(NA, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- tree$tip.label
  pd.dists <- cophenetic.phylo(tree)
  focal.dists <- pd.dists[ , focal] + 1 # avoid 0
  if (psi > 0) focal.dists <- invertVector(focal.dists)
  probs <- focal.dists^abs(psi)
  probs <- probs/sum(probs)
  for (i in 1:nsites) {
    incidence <- abs(ceiling(rnorm(1, mean = mean.incid))) # avoid negative numbers
    output[i, ] <- ifelse(1:ntips %in% sample(1:ntips, incidence,
                                              prob = probs), 1, 0)
  }
  if (mean.abun != FALSE) {
    for (i in 1:nsites) {
      output[i, ] <- genAbuns(output[i, ])
    }
  }
  return (output)
}

#' @name evenCommData
#' @title Generate Evenly Distributed Data Based on Community Phylogeny
#' @description Generate a community matrix based on a phylogeny with phylogenetic signal.
#' @details Generate evenly distributed community phylogenetic data based on a
#' phylogeny. Takes the phylogeny, calculates cummulative PD from random
#' taxon, cuts at even intervals based on nspp to maximise distance.
#' @template base_template
#' @param nsites number of sites
#' @param nspp number of species

evenCommData <- function(tree, nsites, nspp) {
  ## TODO: write test
  ntips <- length(tree$tip.label)
  output <- matrix(rep(0, ntips * nsites),
                   ncol = ntips, nrow = nsites)
  colnames(output) <- tree$tip.label
  pd.dists <- cophenetic.phylo(tree)
  for (i in 1:nsites) {
    focal.pd.dists <- pd.dists[ , sample(1:ntips, 1)]
    # select focal taxon at random
    splits <- .split0(1:sum(focal.pd.dists), nspp)
    cumsum.pd.dists <- cumsum(focal.pd.dists)
    for (j in 1:nspp) {
      pull <- which(abs(cumsum.pd.dists - splits[j])
                    == min(abs(cumsum.pd.dists - splits[j])))
      output[i, pull] <- 1
    }
  }
  return(output)
}

#' @name calcDist
#' @title Calculate the distance between trees using different methods
#' @description Calculate the normalised tree distance (topological and branch length)
#' between trees
#' @details This functions uses three different methods for calculating the distance
#' between trees: the 'PH85' and 'score' methods of ape's topo.dist and/or 1 - Pearson's r
#' calculated from the cophenetic matrix of the tip distances ('dmat'). The function returns
#' the normalised distances by default. Trees of different numbers of tips are shrunk
#' to the same size. Function will also rescale branch distances to sum to 1 before
#' calculating distances.
#' @param tree1 first tree of comparison
#' @param tree2 second tree of comparison
#' @param method 'all', 'PH85', 'score' or 'dmat'. Default 'all'
#' @param normalised boolean, return distances of 0-1

calcDist <- function (tree1, tree2, method = c ('all', 'PH85', 'score', 'dmat'),
                      normalised = TRUE) {
  # internals
  calcPH85 <- function (tree1, tree2) {
    # dist.topo assumes trees are unrooted
    tree1 <- unroot (tree1)
    tree2 <- unroot (tree2)
    # get distance
    d <- dist.topo (tree1, tree2, method = 'PH85')
    if (normalised) {
      # expected max, the sum of total internal branches
      max.d <- (tree1$Nnode + tree2$Nnode) - 2
      d <- d / max.d
    }
    d
  }
  calcScore <- function (tree1, tree2) {
    # If either have no lengths, return NA
    if (is.null (tree1$edge.length) | is.null (tree2$edge.length)) {
      return (NA)
    }
    # dist.topo assumes trees are unrooted
    tree1 <- unroot (tree1)
    tree2 <- unroot (tree2)
    # rescale branch lengths
    tree1$edge.length <- tree1$edge.length/sum (tree1$edge.length)
    tree2$edge.length <- tree2$edge.length/sum (tree2$edge.length)
    # get distance
    d <- dist.topo (tree1, tree2, method = 'score')
    if (normalised) {
      # get internal branches for tree 1 and 2
      ibs.t1 <- which (!tree1$edge[ ,2] %in% 1:length (tree1$tip.label))
      ibs.t2 <- which (!tree2$edge[ ,2] %in% 1:length (tree2$tip.label))
      # get internal branch lengths
      els.1 <- tree1$edge.length[ibs.t1]
      els.2 <- tree2$edge.length[ibs.t2]
      # calc max possible distance
      max.d <- sqrt (sum (els.1^2, els.2^2))
      d <- d / max.d
    }
    d
  }
  calcDmat <- function (tree1, tree2) {
    # If either have no lengths, return NA
    if (is.null (tree1$edge.length) | is.null (tree2$edge.length)) {
      return (NA)
    }
    # get tip distances
    dist1 <- cophenetic.phylo (tree1)
    dist2 <- cophenetic.phylo (tree2)
    # ensure tips are in the same order
    dist1 <- dist1[order (colnames(dist1)), order (colnames(dist1))]
    dist2 <- dist2[order (colnames(dist2)), order (colnames(dist2))]
    model <- cor.test (x = dist1, y = dist2)
    as.numeric (1 - model$estimate)  # 1 - Pearson's r
  }
  method <- match.arg (method)
  # safety first
  if (sum (tree1$tip.label %in% tree2$tip.label) < 3) {
    stop ('Trees must share more than 3 tips.')
  }
  # tip handling
  tree1.drop <- tree1$tip.label[!tree1$tip.label %in% tree2$tip.label]
  if (length (tree1.drop) > 0) {
    tree1 <- drop.tip (tree1, tip = tree1.drop)
  }
  tree2.drop <- tree2$tip.label[!tree2$tip.label %in% tree1$tip.label]
  if (length (tree2.drop) > 0) {
    tree2 <- drop.tip (tree2, tip = tree2.drop)
  }
  if (method == 'PH85') {
    return (calcPH85 (tree1, tree2))
  }
  if (method == 'score') {
    return (calcScore (tree1, tree2))
  }
  if (method == 'dmat') {
    return (calcDmat (tree1, tree2))
  }
  data.frame ('PH85' = calcPH85 (tree1, tree2),
              'score' = calcScore (tree1, tree2),
              'dmat' = calcDmat (tree1, tree2))
}

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
#' relationships for the names and tree given.
#' 
#' If stats is true, the function will return a list of the extracted tree and a data frame
#' giving the proportions of names mapped to the tree, names mapped by exact string matching,
#' fuzzy matching and random taxonomic group placement. Of those names that are mapped using
#' random taxonomic placement, the mean rank number is returned. For the NCBI taxonomy, these
#' rank numbers are roughly equivalent to:
#' 
#'     Rank         Number
#' 
#'     superkingdom    2
#' 
#'     kingdom         4
#' 
#'     phylum          9
#' 
#'     subphylum       10
#' 
#'     superclass      12
#' 
#'     class           18
#' 
#'     superorder      21
#' 
#'     order           22
#' 
#'     suborder        23
#'     
#'     infraorder      24
#'     
#'     parvorder       25
#'     
#'     superfamily     26
#' 
#'     family          27
#' 
#'     subfamily       28
#' 
#'     genus           29
#' 
#'     species         30
#' 
#' These numbers depend on datasource and taxonomic group; some groups have more ranks than others.
#' 
#' If the mean rank number is very high, it is likely that the tree on which the names have been
#' mapped does not represent them very well.
#' 
#' @template base_template
#' @param names vector of names of tips to be extracted
#' @param fuzzy boolean, if true will search Global Names Resolver online
#' @param stats boolean, if true returns list object of extracted tree and extraction statistics
#' @param datasource GNR datasource ID, default 4 for NCBI
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
#' hominid.map <- mapNames (tree=catarrhines, names=names, stats=TRUE, fuzzy=TRUE)
#' plot (hominid.map$tree)
#' print (hominid.map$stats)  # all names were mapped, 4 fuzzily on average at the genus level

mapNames <- function (tree, names, fuzzy=TRUE, stats=FALSE, datasource=4) {
  # early checks
  if (!is.rooted (tree)) {
    stop ('Tree must be rooted.')
  }
  if (!is.vector (names) && length (names) <= 1) {
    stop ('Names must be a vector of more than 1 character string')
  }
  if (is.null (tree$edge.length)) {
    warning ('No branch lengths in tree, branch lengths added with `compute.brlen()`')
    tree <- compute.brlen (tree)
  }
  # randomise names to prevent biasing random tree generated
  names <- sample (names)
  # stats
  ranks <- vector ()  # taxonomic rank where tips have been randomly added
  # drop _
  tree$tip.label <- gsub ('_', ' ', tree$tip.label)
  names <- gsub ('_', ' ', names)
  # finding matching and non-matching names
  matching.names <- names[names %in% tree$tip.label]
  nonmatching.names <- names[!names %in% tree$tip.label]
  # stop now if all names match names in tree or fuzzy is false
  if (!fuzzy || length (matching.names) == length (names)) {
    return (.extract (tree, names, stats, matching.names, ranks))
  }
  # resolve names that didn't string match
  resolved <- taxaResolve (names=nonmatching.names, datasource=datasource)
  # drop NAs
  resolved <- resolved[!is.na (resolved$name.string), ]
  # separate lineages
  resolved.lineages <- strsplit (as.vector (resolved$lineage), '\\|')
  # stop now if none resolved
  if (!nrow (resolved) > 0) {
    return (.extract (tree, names, stats, matching.names, ranks))
  }
  # get parent of matched names to reduce the number of searches needed
  if (length (matching.names) > 1) {
    parent.node <- getParent (tree = tree, tips = matching.names)
  } else {
    # else whole tree
    parent.node <- getSize (tree) + 1
  }
  res <- NULL
  # keep looping and increasing the sample size until no longer possible
  while (nrow (resolved) > 0) {
    # get resolved names for all tips in tree
    res <- .searchTreeNames (tree=tree, parent.node=parent.node,
                                       previous=res, datasource=datasource)
    if (nrow (res$resolved) > 0) {
      # calculate min rank
      bool <- rep (TRUE, length (res$lineages[[1]]))
      for (lineage in res$lineages[2:length (res$lineages)]) {
        bool <- bool & (res$lineages[[1]] %in% lineage)
      }
      min.rank <- which (bool)  # this ensures the sampled tree is not too narrow
      min.rank <- min.rank[length (min.rank)]
      # loop through each in resolved and match to tree subset
      drop.vector <- rep (FALSE, nrow (resolved))  # bool
      for (i in 1:nrow (resolved)) {
        lineage <- resolved.lineages[[i]]
        # loop through each in tree lineages to find the best matching tip
        matches <- rep (NA, length (res$lineages))
        for (j in 1:length (res$lineages)) {
          matches[j] <- max (which (lineage %in% res$lineages[[j]]))
        }
        if (max (matches, na.rm=TRUE)[1] >= min.rank) {
          possibles <- as.vector (which (matches == max (matches, na.rm=TRUE)))
          if (length (possibles) > 1) {
            # choose one at random
            match <- sample (possibles, 1)
          } else {
            match <- possibles
          }
          tree <- .addTip (tree=tree, tip.i=res$resolved$tip.i[match],
                               new.name=as.character (resolved$search.name[i]))
          # get parent node + 1
          parent.node <- parent.node + 1
          # update tip.i as new tip is added to tree
          res <- .searchTreeNames (tree=tree, parent.node=parent.node,
                                   previous=res, datasource=datasource)
          drop.vector[i] <- TRUE
          ranks <- append (ranks, max (matches, na.rm=TRUE))
        }
      }
      # drop resolved names now accounted for
      resolved <- resolved[!drop.vector, ]
    }
    # if parent node is root, break out of while loop
    if (parent.node == getSize (tree) + 1) {
      break
    }
    parent.node <- getParent (tree, node=parent.node)
  }
  return (.extract (tree, names, stats, matching.names, ranks))
}

# HIDDEN mapNames functions
.searchTreeNames <- function (tree, parent.node, previous=NULL,
                              datasource=datasource) {
  # search and return resolved names for a tree using a subset of names
  # based on parent.node
  # Returns a taxaResolve() data.frame
  if (!is.null (previous)) {
    deja.vues <- as.vector (previous$resolved$search.name)
  } else {
    deja.vues <- vector ()
  }
  # get subset of tree and resolve
  subset.tree <- extract.clade (tree, node=parent.node)
  # exclude names already searched
  subset.names <- subset.tree$tip.label[!subset.tree$tip.label %in% deja.vues]
  if (length (subset.names) > 1) {
    res <- taxaResolve (names=subset.names, datasource=datasource)
    # remove NA results
    res <- res[!is.na (res$name.string), ]
    # stick to previous results
    if (!is.null (previous)) {
      res$tip.i <- rep (NA, nrow (res))
      res <- rbind (previous$resolved, res)
    }
  } else {
    res <- previous$resolved
  }
  # add tip.is
  res$tip.i <- match (res$search.name, tree$tip.label)
  lineages <- strsplit (as.vector (res$lineage), '\\|')
  list ('resolved'=res, 'lineages'=lineages)
}

.addTip <- function (tree, tip.i, new.name) {
  # add new tip to new tree on pendant edge of matching tip
  edge <- which (tip.i == tree$edge[ ,2])
  # random node age somewhere on pendant edge
  node.age <- runif (1, 0, tree$edge.length[edge])
  tree <- addTip (tree, edge=edge, tip.name=new.name, node.age=node.age)
  return (tree)
}

.extract <- function (tree, names, stats, matching.names, ranks) {
  # return tree representing names, also stats
  if (sum (names %in% tree$tip.label) > 1) {
    res.tree <- drop.tip (tree, tip = tree$tip.label[!tree$tip.label %in% names])
    if (!stats) {
      return (res.tree)
    } else {
      p.names <- sum (names %in% tree$tip.label)/length (names)
      n.matched <- length (matching.names)
      n.fuzzy <- length (ranks)
      mean.rank <- mean (ranks)
      res.stats <- data.frame (p.names, n.matched, n.fuzzy, mean.rank)
      return (list ('tree'=res.tree, 'stats'=res.stats))
    }
  } else {
    warning ('Too few names could be mapped to tree')
    return (NA)
  }
}