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
#' @title Calculate Evolutionary Distinctiveness (Fair Proportion or Pendant Edge)
#' @description Calculate species evolutionary distinctiveness (ED) using either of two methods:
#'  Fair Proportion (FP) or Pendant Edge (PE).
#' @details Evolutionary distinctiveness is a measure of how much independent evolution a species
#'  represents. Multiple methods exist for its calculation all of which require an ultrametric
#'  phylogenetic tree. The methods used here are Pendant Edge (PE) the length of a species branch
#'  that connects it to the tree (Altschul and Lipman, 1990) or Fair Proportion (FP) the total
#'  proportion of the phylogenetic tree that a species represents where each branch is equally
#'  divided between all descdendants (Isaac et al. 2007).
#'  
#'  N.B. \code{picante} already has a function for doing this. \code{calcED}, however, uses \code{plyr}
#'  vectorisation and should be faster.
#' @template base_template
#' @param tips vector of tips for which ED is calculated, else 'all'
#' @param type method of ED calculation, either 'FP' or 'PE'
#' @param progress bar showing the progress of the calculation. Either
#' time, text, tk or none. Uses \code{plyr}'s \code{create_progress_bar}.
#' @references Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). 
#'  Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#'  
#'  Altschul, S., and Lipman, D. (1990). Equal animals. Nature.
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

calcED <- function (tree, tips = 'all', type = 'FP', progress = 'none') {
  calcFairProportion <- function (tips) {
    countChildren <- function (node) {
      length (getChildren (tree, node))
    }
    calcSpecies <- function (sp) {
      edges <- getEdges (tree, tips = as.character (sp), type = 2)
      n.children = mdply (.data = data.frame (node = tree$edge[edges, 2]),
                          .fun = countChildren)[ ,2]
      sum (tree$edge.length[edges]/n.children)
    }
    EDs <- mdply (.data = data.frame (sp = tips),
                  .fun = calcSpecies,
                  .progress = create_progress_bar (name = progress))
    EDs
  }
  calcPendantEdge <- function (tips) {
    getEdgeLength <- function (sp) {
      tip.index <- which (tree$tip.label == sp)
      tip.edge <- which (tree$edge[, 2] == tip.index)
      tree$edge.length[tip.edge]
    }
    mdply (.data = data.frame (sp = tips), .fun = getEdgeLength,
           .progress = create_progress_bar (name = progress))
  }
  if (tips[1] == 'all') {
    tips <- tree$tip.label
  }
  if (type == 'FP') {
    EDs <- calcFairProportion (tips)
  } else if (type == 'PE') {
    EDs <- calcPendantEdge (tips)
  } else {
    stop (paste0 ('Type [', type, '] not known.
                  Type must be either FP or PE'))
  }
  colnames (EDs)[2] <- 'ED'
  EDs
}