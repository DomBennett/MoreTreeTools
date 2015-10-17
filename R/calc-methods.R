#' @name calcEdgeDiversity
#' @title Calculate the diversity of edges within set time intervals
#' @description Return a dataframe of edge number and a count
#' based on the number descendents within n.intervals uniformly
#' spaced time intervals.
#' @details Intervals are unifromely spaced from the root to the tip
#' of the tree. Output can be provided to treeplot for colouring edges,
#' see example.
#' @template base_template
#' @param n.intervals number of intervals
#' @export
#' @examples
#' tree <- rtree (100)
#' edge.diversity <- calcEdgeDiversity (tree, n.intervals=4)
#' # convert to logged z-score to increase colour spectrum in plot
#' edge.diversity$count <- log (edge.diversity$count)
#' edge.diversity$col <- (edge.diversity$count - mean (edge.diversity$count)) / sd (edge.diversity$count)
#' treeplot (tree, edge.cols=edge.diversity, legend.title='Diversity')

calcEdgeDiversity <- function (tree, n.intervals) {
  # internal
  .count <- function (i) {
    .findIgnoreTips <- function (tip) {
      node <- which (tip == tree$tip.label)
      tip.age <- getAge (tree, node=node)
      if (tip.age$age < t1) {
        ignore.tips <<- c (tip, ignore.tips)
      }
    }
    .accountForCutters <- function (cutter) {
      connecting.nodes <-
        getNodes (tree, node=edge.stats[cutter, 'node.2'])
      n.children <<- n.children +
        as.numeric (edge.stats$node.2 %in% connecting.nodes)
    }
    t0 <- ts[i-1]
    t1 <- ts[i]
    # find all edges within age range
    edges <- which (edge.stats.ref$age.1 <= t0 &
                      edge.stats.ref$age.1 >= t1)
    if (length (edges) > 1) {
      # find all tip nodes younger than t1
      ignore.tips <- NULL
      m_ply (.data=data.frame (tip=tree$tip.label, stringsAsFactors=FALSE),
             .fun=.findIgnoreTips)
      # calc edge stats ignoring ignore.tips
      edge.stats <- getEdgeStats (tree, edges=edges,
                                  ignore.tips=ignore.tips)
      # identify all edges that cut through time interval
      # add 1 to all connecting nodes to the cutters
      cutters <- edge.stats$age.1 > t1 & edge.stats$age.2 < t1
      n.children <- edge.stats$n.children
      m_ply (.data=data.frame (cutter=which (cutters)),
             .fun=.accountForCutters)
      n.children <- n.children + as.numeric (cutters)
      # keep edge and new count + edge label if available
      temp <- edge.stats[ ,!colnames (edge.stats) %in% 
                         c ('node.1', 'node.2', 'age.1', 'age.2', 'n.children')]
      temp$count <- n.children
      return (temp)
    }
  }
  if (n.intervals < 2) {
    stop ('n.intervals must be greater than or equal to 2')
  }
  edge.stats.ref <- getEdgeStats (tree)
  tree.age <- getSize (tree, 'rtt')
  ts <- seq (from=tree.age, to=0, length.out=n.intervals+1)
  mdply (.data=data.frame (i=2:length (ts)),
         .fun=.count)[ ,-1]
}

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
  if (any (is.na (cmatrix))) {
    stop ('NAs in cmatrix')
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
#' @title Calculate Evolutionary Distinctness
#' @description Calculate species evolutionary distinctness (ED) using one of three methods:
#'  Fair Proportion (FP), Equal Splits (ES) or Pendant Edge (PE).
#' @details Evolutionary distinctness is a measure of how much independent evolution a species
#'  represents. Multiple methods exist for its calculation all of which require an ultrametric
#'  phylogenetic tree. The methods used here are Pendant Edge (PE) the length of a species branch
#'  that connects it to the tree (Altschul and Lipman, 1990), Fair Proportion (FP) the total
#'  proportion of the phylogenetic tree that a species represents where each branch is equally
#'  divided between all descdendants (Isaac et al. 2007) and Equal Splits (ES) where branch lengths
#'  are equally divided between descendents at every node in the tree (Redding and Mooers 2006)
#' @template base_template
#' @param tips vector of tips for which ED is calculated, else 'all'
#' @param type method of ED calculation, either 'all', 'FP', 'ES' or 'PE' (default FP)
#' @seealso
#' \code{\link[picante]{evol.distinct}}, \code{\link{catarrhines}}
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
#' # Load catarrhine tree
#' data(catarrhines)
#' ed.vals <- calcED (catarrhines, type='all')
#' # Which Old World monkeys are the most distinct?
#' ed.vals[ed.vals$FP == min (ed.vals$FP), ]
#' ed.vals[ed.vals$ES == min (ed.vals$ES), ]
#' ed.vals[ed.vals$PE == min (ed.vals$PE), ]

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
                    ncol=ntips, nrow=nsites)
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
#' @details This functions uses four different methods for calculating the distance
#' between trees: the 'PH85' and 'score' methods of ape's topo.dist, 1 - Pearson's r
#' calculated from the cophenetic matrix of the tip distances ('dmat') or the triplet
#' distance of Critchlow et al. (1996). The function returns the normalised distances
#' by default. Trees of different numbers of tips are pruned to the same size.
#' Function will also rescale branch distances to sum to 1 before
#' calculating distances. Note, triplet metric will only work for bifurcating and
#' rooted trees otherwise it will return NA.
#' @param tree1 first tree of comparison
#' @param tree2 second tree of comparison
#' @param method 'all', 'PH85', 'score', 'trip' or 'dmat'. Default 'all'
#' @param normalised boolean, return distances of 0-1
#' @references
#' Critchlow DE, Pearl DK, Qian C. (1996) The Triples Distance for rooted bifurcating
#' phylogenetic trees. Systematic Biologly, 45, 323–34.
#' 
#' Billera, L. J., Holmes, S. P. and Vogtmann, K. (2001) Geometry of the space
#' of phylogenetic trees. Advances in Applied Mathematics, 27, 733–767.
#' 
#' Kuhner, M. K. and Felsenstein, J. (1994) Simulation comparison of phylogeny
#' algorithms under equal and unequal evolutionary rates. Molecular Biology and
#' Evolution, 11, 459–468.
#' 
#' Nei, M. and Kumar, S. (2000) Molecular Evolution and Phylogenetics. Oxford:
#' Oxford University Press.
#' 
#' Penny, D. and Hendy, M. D. (1985) The use of tree comparison metrics.
#' Systemetic Zoology, 34, 75–82.
#' 
#' Rzhetsky, A. and Nei, M. (1992) A simple method for estimating and testing
#' minimum-evolution trees. Molecular Biology and Evolution, 9, 945–967.

calcDist <- function (tree1, tree2, method=c ('all', 'PH85', 'score', 'dmat', 'trip'),
                      normalised=TRUE) {
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
    return (calcPH85 (tree1, tree2, normalised))
  }
  if (method == 'score') {
    return (calcScore (tree1, tree2, normalised))
  }
  if (method == 'dmat') {
    return (calcDmat (tree1, tree2, normalised))
  }
  if (method == 'trip') {
    return (calcTrip (tree1, tree2, normalised))
  }
  data.frame ('PH85'=calcPH85 (tree1, tree2, normalised),
              'score'=calcScore (tree1, tree2, normalised),
              'triplet'=calcTrip (tree1, tree2, normalised),
              'dmat'=calcDmat (tree1, tree2, normalised))
}

calcPH85 <- function (tree1, tree2, normalised) {
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

calcScore <- function (tree1, tree2, normalised) {
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
calcDmat <- function (tree1, tree2, normalised) {
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

calcTrip <- function (tree1, tree2, normalised) {
  # internal
  .count <- function (i) {
    tips <- all.tips[ ,i]
    tree1.outgroup <- getOutgroup (tree1, tips)
    tree2.outgroup <- getOutgroup (tree2, tips)
    # test if the two are the same
    if (length (tree1.outgroup) != length (tree2.outgroup) ||
          tree1.outgroup != tree2.outgroup) {
      counter <<- counter + 1
    }
  }
  # requires bifurcating and rooted trees
  if (!is.rooted (tree1) && !is.rooted (tree2)) {
    return (NA)
  }
  if (!is.binary.tree (tree1) && !is.binary.tree (tree2)) {
    return (NA)
  }
  # shared tips
  pull.1 <- tree1$tip.label %in% tree2$tip.label
  pull.2 <- tree2$tip.label %in% tree1$tip.label
  # get all combinations of 3 tips
  all.tips <- combn (tree1$tip.label[pull.1], 3)
  # loop through all combinations,
  # count all triplets where the outgroups differ
  counter <- 0
  m_ply (.data=data.frame(i=1:ncol (all.tips)),
         .fun=.count)
  if (normalised) {
    return (counter / ncol (all.tips))
  } else {
    return (counter)
  }
}

#' @name calcBalance
#' @title Calculate balance of tree
#' @description Return data.frame of nodal imbalance
#' @template base_template
#' @param node which node to calculate balance from, default 'all'
#' @export

calcBalance <- function (tree, node = 'all') {
  # Return for a given node:
  #  the absolute difference of an even split [absdiff]
  #  the normalised balance (0 - 1) [nabsdiff]
  #  the number of transitions to maximal balance [b.steps]
  #  the number of transitions to imbalance [i.steps]
  .calc <- function (node) {
    following.edges <- tree$edge[ ,1] == node
    sister.nodes <- tree$edge[following.edges, 2]
    ntot <- length (getChildren (tree, node = node))
    n1 <- length (getChildren (tree, node = sister.nodes[1]))
    max.absdiff <- (ntot / 2) - 1  # max absdiff possible
    absdiff <- abs ((ntot / 2) - n1)
    b.steps <- absdiff %/% 1
    i.steps <- (max.absdiff - b.steps) %/% 1
    nabsdiff <- b.steps / (b.steps + i.steps)
    data.frame (n=ntot, absdiff, nabsdiff, b.steps, i.steps)
  }
  if (node == 'all') {
    node <- (getSize (tree) + 1):
      (getSize (tree) + (tree$Nnode))
  }
  by.node <- mdply (.data = data.frame (node=node), .fun = .calc)
  # total number of steps to balance in tree
  b.steps <- sum (by.node$b.steps[1:(nrow (by.node) - 1)] !=
                    by.node$b.steps[2:nrow (by.node)])
  # total number of steps to imbalance in tree
  i.steps <- sum (by.node$i.steps[1:(nrow (by.node) - 1)] !=
                    by.node$i.steps[2:nrow (by.node)])
  by.tree <- data.frame (b.steps, i.steps)
  list (by.node=by.node, by.tree=by.tree)
}