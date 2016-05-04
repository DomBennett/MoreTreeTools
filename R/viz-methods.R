#' @name commplot
#' @title Plot a community or traits on a phylogenetic tree
#' @description Plot community or trait data on phylogeny to visualise
#' differences in structure between traits or sites. Use colours to distinguish
#' groups and alpha to distinguish abundances (works like rainbow())
#' @details No details
#' @template base_template
#' @param cmatrix community/trait data matrix (cols taxa, rows sites)
#' @param groups trait or site groups
#' @export
#' @examples
#' # Generate a balanced tree
#' tree <- compute.brlen (stree (16, 'balanced'))
#' # Generate random community data, 2 of 3 different communities
#' cmatrix <- matrix (round (runif (16*6, 0, 2)), nrow = 6, byrow = TRUE)
#' colnames (cmatrix) <- tree$tip.label
#' rownames (cmatrix) <- paste0 ('site_', 1:6)
#' # Plot community, specify group identitiy of each of the 6 sites
#' commplot (cmatrix, tree, groups = c (1,1,2,2,3,3))

commplot <- function(cmatrix, tree, groups = rep(1, nrow(cmatrix)),
                     no.margin = TRUE, ...){
  # ultrametricize tree FIRST
  tree <- compute.brlen(tree, method="Grafen")
  # plot phylogeny, allow space for points
  tree$edge.length <- tree$edge.length/getSize (tree, 'rtt')
  # make all trees the same length before plotting i.e. all branches from terminal
  # node to tip equal 1
  # for some weird reason the rules of plotting are dyanmic!
  if (nrow(cmatrix) < 20) {
    variable.max <- 1 + nrow(cmatrix)/20
    spacing.start <- 0.55
    spacing.i <- 0.05
  } else {
    variable.max <- nrow(cmatrix)/10
    spacing.i <- 0.1 - 1/nrow(cmatrix)
    spacing.start <- 0.5 + spacing.i
  }
  plot(tree, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max), ...)
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- cmatrix/max(cmatrix)
  # loop init
  ntips <- length(tree$tip.label)
  spacing <- spacing.start
  group <- groups[1]
  j <- 1
  # loop through sites and plot points for present species
  for(i in 1:nrow(cmatrix)){
    j <- ifelse(group == groups[i], j, j + 1)
    pull <- as.logical(cmatrix[i,])
    taxa <- tree$tip.label[pull]
    abunds <- alphas[i, pull]
    tiplabels(tip = match(taxa, tree$tip.label),
              pch = 19, adj = spacing, col = hsv(rep(hs[j], ntips), 1, 1, abunds))
    spacing <- spacing + spacing.i
    group <- groups[i]
  }
}

#' @name chromatophylo
#' @title Plot tree through time with coloured edges
#' @description Return a geom_object() of a phylogenetic tree through time.
#' Requires tree given to be rooted, time callibrated and bifurcating.
#' @details No details -- see details.
#' @param tree
#' @param edge.cols data.frame of values used to colour edges (see example)
#' @param edge.sizes data.frame of sizes used to determine the size of edges
#' @param legend.title text for legend if edge.cols provided
#' @param reduce.overlap boolean, TRUE will prevent edges from overlapping in plot, may take longer to plot
#' @export
#' @examples
#' # Quick example without colours
#' tree <- rtree (50)  # random tree with fossils
#' plot (tree)  # ape plot
#' chromatophylo (tree)  # chromatophylo
#' 
#' 
#' # Example with collapseTips() + calcEdgeDiversity()
#' # generate random tree
#' tree <- rtree (100)
#' # add edge labels to track edges after collapse
#' tree$edge.label <- paste0 ('edge_', 1:nrow (tree$edge))
#' # generate edge diversity using two intervals
#' ed <- calcEdgeDiversity (tree, n.intervals=2)
#' # iteratively drop all tip edges with > 10% tree age
#' tree.age <- getSize (tree, 'rtt')
#' ctree <- collapseTips (tree, min.length=tree.age*0.1, iterative=TRUE)
#' # reconstruct edge diversity data frame without dropped edges
#' new.ed <- data.frame (edge.label=ctree$edge.label)
#' new.ed$edge <- 1:nrow (new.ed)
#' match.counts <- match (ctree$edge.label, ed$edge.label)
#' new.ed$count <- ed$count[match.counts]
#' # convert count to logged Z-score for maximum colour separation
#' new.ed$col <- (log (new.ed$count) - mean (log (new.ed$count))) /
#'   sd (log (new.ed$count))
#' ed$col <- (log (ed$count) - mean (log (ed$count))) /
#'   sd (log (ed$count))
#' # plot before
#' chromatophylo (tree, edge.cols=ed, legend.title='Diversity') + ggtitle ('Before collapse')
#' # plot after
#' chromatophylo (ctree, edge.cols=new.ed, legend.title='Diversity') + ggtitle ('After collapse')

# TODO:
# 1. Take tree data from Newick file efficiently
# 3. add tips option

chromatophylo <- function (tree, edge.cols=NULL, edge.sizes=NULL, legend.title='',
                           reduce.overlap=TRUE) {
  if (!is.binary.tree (tree)) {
    stop ('Tree must be bifurcating')
  }
  # init plot data
  N <- length (tree$tip.label)
  tree.age <- getSize (tree, 'rtt')
  x.spacer <- tree.age/100
  y.spacer <- 1/N
  x1 <- tree.age
  node <- length (tree$tip.label) + 1
  edge <- which (tree$edge[ ,1] == node)
  # init p.env
  p.env <- new.env (parent=emptyenv ())
  p.env$N <- N
  p.env$all.edges <- 1:nrow (tree$edge)
  p.env$tree.stats <- getTreeStats (tree)
  p.env$p.data <- data.frame (x=c(x1, x1), y=c(0, 1),
                              edge=c (0,0))
  p.env$y.spacer <- y.spacer
  p.env$x.spacer <- x.spacer
  p.env$tree <- tree
  p.env$root.node <- length (tree$tip.label) + 1
  #t.data <- data.frame (x=NA, y=NA, label=NA)
  # generate p.data
  .cpMkPData (x1, c (0, 1), edge, p.env)
  # find all lines
  p.env$p.data$line <- rep (1:(nrow (p.env$p.data)/2), each=2)
  p.env$lines <- unique (p.env$p.data$line)[-1]  # ignore root
  edges <- p.env$p.data$edge[match (p.env$lines, p.env$p.data$line)]
  # get max and min
  p.env$maxmin <- data.frame (line=p.env$lines, edge=edges,
                              max.x=NA, min.x=NA,
                              max.y=NA, min.y=NA)
  .cpGetMaxMin (p.env$lines, p.env)
  while (reduce.overlap) {
    # loop through p.data rearranging y coords to avoid overlap
    p.env$overlap <- FALSE
    # loop through ignoring root
    plyr::m_ply (.data=data.frame (l=p.env$lines), .fun=.cpCheckLine, p.env=p.env)
    #return (p.env$p.data)
    if (!p.env$overlap) {
      break
    }
  }
  # make geom_object
  p.data <- p.env$p.data
  if (!is.null (edge.cols)) {
    matched.is <- match (p.data$edge, edge.cols$edge)
    p.data$col <- edge.cols$col[matched.is]
  }
  if (!is.null (edge.sizes)) {
    matched.is <- match (p.data$edge, edge.sizes$edge)
    p.data$size <- edge.sizes$col[matched.is]
  }
  p <- ggplot (p.data, aes (x=x, y=y, group=edge))
  colscale <- scale_colour_gradient2 (low="blue", mid='white',
                                      high='red', name=legend.title)
  if (all (c ('col', 'size') %in% names (p.data))) {
    p <- p + geom_line (aes (colour=col, size=size)) + colscale
  } else if ('size' %in% names (p.data)) {
    p <- p + geom_line (aes (size=size), colour='white')
  } else if ('col' %in% names (p.data)) {
    p <- p + geom_line (aes (colour=col), size=2)  + colscale
  } else {
    p <- p + geom_line (colour='white', size=2)
  }
  # set ggplot environment parameters
  theme.opts <- theme(panel.background = element_rect(fill="gray15"),
                      axis.title.y = element_blank(),
                      axis.ticks.y =element_blank(),
                      axis.text.y = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.major.x = element_line(colour="black"),
                      panel.grid.minor.x = element_line(colour="gray5"))
  # add time and reverse to end at 0
  p <- p + scale_x_reverse() + xlab ('Time') + theme.opts
  p
}
.cpGetMaxMin <- function (lines, p.env) {
  # Write max min for x and y for lines given
  .maxmin <- function (l) {
    # get line coords
    which.line <- p.env$p.data$line == l
    line <- p.env$p.data[which.line, ]
    max.x <- max (line$x) + p.env$x.spacer/4
    min.x <- min (line$x) - p.env$x.spacer/4
    min.y <- min (line$y) - p.env$y.spacer/4
    max.y <- max (line$y) + p.env$y.spacer/4
    res <- data.frame ('max.x'=max.x, 'min.x'=min.x,
                       'max.y'=max.y, 'min.y'=min.y)
    which.line <- p.env$maxmin$line == l
    p.env$maxmin[which.line, 3:6] <- c (max.x, min.x, max.y, min.y)
  }
  l.data <- data.frame (l=lines)
  plyr::m_ply (.data=l.data, .fun=.maxmin)
}
.cpCorrectOverlap <- function (edge.1, edge.2, p.env) {
  # set overlap to TRUE
  p.env$overlap <- TRUE
  # find all edges descending from higher edge and add y.spacer/2
  # find all edges descending from lower edge and remove y.spacer/2
  # find parent vertical edge and add and remove y.spacer/2
  parent.node <- .cpGetParentNode (edge.1, edge.2, p.env)
  if (parent.node == p.env$root.node) {
    parent.edge <- 0
  } else {
    parent.edge <- which (p.env$tree$edge[ ,2] == parent.node)
  }
  # get next edges and their ys
  edges <- p.env$tree.stats[[parent.node]][['next.edges']]
  edge.1 <- p.env$p.data[p.env$p.data$edge == edges[1], ]
  edge.2 <- p.env$p.data[p.env$p.data$edge == edges[2], ]
  if (mean (edge.1$y) >= mean (edge.2$y)) {
    higher <- edge.1$edge[1]
    lower <- edge.2$edge[1]
  } else {
    higher <- edge.2$edge[1]
    lower <- edge.1$edge[1]
  }
  # higher
  higher.node <- p.env$tree$edge[higher, 2]
  higher.edges <- c (p.env$tree.stats[[higher.node]][['descend.edges']],
                     higher)
  p.env$p.data[p.env$p.data$edge %in% higher.edges, 'y'] <-
    p.env$p.data[p.env$p.data$edge %in% higher.edges, 'y'] + p.env$y.spacer/2
  # lower
  lower.node <- p.env$tree$edge[lower, 2]
  lower.edges <- c (p.env$tree.stats[[lower.node]][['descend.edges']],
                    lower)
  p.env$p.data[p.env$p.data$edge %in% lower.edges, 'y'] <-
    p.env$p.data[p.env$p.data$edge %in% lower.edges, 'y'] - p.env$y.spacer/2
  # vertical
  vertical.bool <- which (p.env$p.data$edge == parent.edge)
  vertical.bool <- vertical.bool[(length (vertical.bool) - 1):
                                   length (vertical.bool)]
  p.env$p.data[vertical.bool[2], 'y'] <- p.env$p.data[vertical.bool[2], 'y'] + p.env$y.spacer/2
  p.env$p.data[vertical.bool[1], 'y'] <- p.env$p.data[vertical.bool[1], 'y'] - p.env$y.spacer/2
  # get all lines where maxmin have changed
  p.env$maxmin$line[p.env$maxmin$edge %in% c (higher.edges, lower.edges)]
}
.cpCheckLine <- function (l, p.env) {
  # Check if line overlaps with other lines
  which.line <- p.env$p.data$line == l
  line <- p.env$p.data[which.line, ]
  # find connecting edges
  edge <- line$edge[1]
  node <- p.env$tree$edge[edge, 2]
  connecting.edges <- c (p.env$tree.stats[[node]][['ascend.edges']],
                         p.env$tree.stats[[node]][['descend.edges']])
  # check and correct overlap
  potential.edges <- p.env$all.edges[!p.env$all.edges %in% connecting.edges]
  if (length (potential.edges) == 0) {
    overlap <- FALSE
  } else {
    overlap <- TRUE
  }
  sbjct <- p.env$maxmin[p.env$maxmin$line == l, ]
  while (TRUE) {
    # check for overlap, correct, check again
    qry <- p.env$maxmin[p.env$maxmin$edge %in% potential.edges, ]
    overlap <- .cpCheckOverlap (sbjct, qry)
    if (any (overlap)) {
      overlapping.edge <- qry$edge[overlap][1]
      changed.lines <- .cpCorrectOverlap (edge.1=edge, edge.2=overlapping.edge, p.env)
      # recalc maxmin for lines that have moved
      .cpGetMaxMin (changed.lines, p.env)
    } else {
      break
    }
  }
}
.cpCheckOverlap <- function (sbjct, qry) {
  # Check if sbjct (line) and qry (multiple potential lines)
  #  overlap using p.env$maxmin
  # unpack
  maxsx <- sbjct[ ,'max.x']
  minsx <- sbjct[ ,'min.x']
  maxsy <- sbjct[ ,'max.y']
  minsy <- sbjct[ ,'min.y']
  maxqx <- qry[ ,'max.x']
  minqx <- qry[ ,'min.x']
  maxqy <- qry[ ,'max.y']
  minqy <- qry[ ,'min.y']
  # run bools
  xspace <- maxsx >= minqx & minsx <= maxqx
  yspace <- maxsy >= minqy & minsy <= maxqy
  ycross <- maxsy >= maxqy & minsy < minqy
  xcross <- maxsx >= maxqx & minsx <= minqx
  xspace & ycross | yspace & xcross |
    yspace & xspace
}
.cpMkPData <- function (x1, y, edge, p.env) {
  if (length (edge) == 2) {
    # run recursively for both
    .cpMkPData (x1, y[1], edge[1], p.env)
    .cpMkPData (x1, y[2], edge[2], p.env)
  } else if (length (edge) == 1){
    # find the next node on the edge
    next.node <- p.env$tree$edge[edge,2]
    # get n.children
    n.children <- p.env$tree.stats[[next.node]][['n.children']]
    # get the length of the edge
    x2 <- x1 - p.env$tree$edge.length[edge]
    l.data <- data.frame (x=c(x1, x2), y=c(y, y), edge)
    p.env$p.data <- rbind (p.env$p.data, l.data)
    # run again for next edge
    next.edge <- p.env$tree.stats[[next.node]][['next.edges']]
    if (next.node <= length (p.env$tree$tip.label)) {
      # EDIT THIS TO ADD TIP LABELS
      #t.data <<- data.frame (x=(x2+1/N), y=y,
      #                       label=tree$tip.label[next.node])
    } else {
      # use n.children to determine y and reduce overlap
      y.length <- n.children / p.env$N
      ys <- c (y-(y.length/2), y+(y.length/2))
      # add vertical connector
      l.data <- data.frame (x=c(x2, x2), y=ys, edge)
      p.env$p.data <- rbind (p.env$p.data, l.data)
      # move on to next edge
      .cpMkPData (x2, ys, next.edge, p.env)
    }
  } else {
    stop ('Invalid tree')
  }
}
.cpGetParentNode <- function (edge.1, edge.2, p.env) {
  # get nodes
  node.1 <- p.env$tree$edge[edge.1, 2]
  node.2 <- p.env$tree$edge[edge.2, 2]
  # find ascending nodes
  nodes.1 <- p.env$tree.stats[[node.1]][['ascend.nodes']]
  nodes.2 <- p.env$tree.stats[[node.2]][['ascend.nodes']]
  # match to find highest shared
  maxi <- which.max (match (nodes.1, nodes.2))
  nodes.1[maxi]
}