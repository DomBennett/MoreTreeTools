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
  gnr.obj <- .taxaResolve (names, datasource = datasource)
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