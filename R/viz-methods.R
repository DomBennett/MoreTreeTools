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
#' @details No details
#' @param tree
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
# 2. reverse x axis
# 3. add tips option
# 4. test
# 5. speed up

chromatophylo <- function (tree, edge.cols=NULL, edge.sizes=NULL, legend.title='') {
  if (!is.binary.tree (tree)) {
    stop ('Tree must be bifurcating')
  }
  # init plot data
  N <- length (tree$tip.label)
  tree.age <- getSize (tree, 'rtt')
  x.length <- tree.age/100
  y.length <- 1
  x1 <- 0
  node <- length (tree$tip.label) + 1
  edge <- which (tree$edge[ ,1] == node)
  # init p.env
  p.env <- new.env (parent=emptyenv ())
  p.env$tree.stats <- getTreeStats (tree)
  p.env$p.data <- data.frame (x=c(0, 0), y=c(0, y.length),
                              edge=c (0,0))
  p.env$y.length <- y.length
  p.env$x.length <- x.length
  p.env$tree <- tree
  p.env$root.node <- length (tree$tip.label) + 1
  #t.data <- data.frame (x=NA, y=NA, label=NA)
  # generate p.data
  .cpMkLine (x1, c (0, y.length), edge, p.env)
  while (TRUE) {
    # loop through p.data rearranging y coords to avoid overlap
    p.env$overlap <- FALSE
    p.env$lines <- seq (from=3, to=nrow (p.env$p.data), by=2)
    m_ply (.data=data.frame (sl=lines), .fun=.cpEachSbjct, p.env=p.env)
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
  
  p <- p + xlab ('Time') + theme.opts
  p
}
.cpEachQry <- function (ql, p.env) {
  # get qry
  qry <- p.env$p.data[ql:(ql+1), ]
  ql.edge <- unique (qry$edge)
  if (!ql.edge %in% p.env$avoid.edges) {
    overlap <- .cpCheckOverlap (p.env$sbjct, qry,
                                p.env$y.length, p.env$x.length)
    if (overlap) {
      # set overlap to TRUE
      p.env$overlap <- TRUE
      # shift ys
      p.env$p.data <- .cpShiftY (p.env)
    }
  }
}
.cpEachSbjct <- function (sl, p.env) {
  p.env$sbjct <- p.env$p.data[sl:(sl+1), ]
  # find connecting edges
  p.env$sl.edge <- unique (sbjct$edge)
  nodes <- p.env$tree$edge[sl.edge, ]
  p.env$avoid.edges <- c (which (p.env$tree$edge[ ,2] == nodes[1]),
                          which (p.env$tree$edge[ ,1] == nodes[2]),
                          sl.edge)
  m_ply (.data=data.frame (ql=p.env$lines), .cpEachQry, p.env=p.env)
}
.cpCheckOverlap <- function (sbjct, qry, y.length, x.length) {
  # check if sbjct and qry straight lines overlap (within y.length)
  maxqx <- max (qry$x) + x.length/4
  minqx <- min (qry$x) - x.length/4
  minqy <- min (qry$y) - y.length/4
  maxqy <- max (qry$y) + y.length/4
  maxsx <- max (sbjct$x) + x.length/4
  minsx <- min (sbjct$x) - x.length/4
  minsy <- min (sbjct$y) - y.length/4
  maxsy <- max (sbjct$y) + y.length/4
  xspace <- maxsx > minqx & minsx < maxqx
  yspace <- maxsy > minqy & minsy < maxqy
  ycross <- maxsy > maxqy & minsy < minqy
  xcross <- maxsx > maxqx & minsx < minqx
  xspace & ycross | yspace & xcross |
    yspace & xspace
}
.cpMkLine <- function (x1, y, edge, p.env) {
  if (length (edge) == 2) {
    # run recursively for both
    .cpMkLine (x1, y[1], edge[1], p.env)
    .cpMkLine (x1, y[2], edge[2], p.env)
  } else if (length (edge) == 1){
    # find the next node on the edge
    next.node <- p.env$tree.stats[[node]][['next.nodes']]
    # get n.children
    n.children <- p.env$tree.stats[[node]][['n.children']]
    # get the length of the edge
    x2 <- x1 + p.env$tree$edge.length[edge]
    l.data <- data.frame (x=c(x1, x2), y=c(y, y), edge)
    p.env$p.data <- rbind (p.env$p.data, l.data)
    # run again for next edge
    next.edge <- p.env$tree.stats[[node]][['next.edges']]
    if (next.node <= length (p.env$tree$tip.label)) {
      # EDIT THIS TO ADD TIP LABELS
      #t.data <<- data.frame (x=(x2+1/N), y=y,
      #                       label=tree$tip.label[next.node])
    } else {
      ys <- c (y-(p.env$y.length/2), y+(p.env$y.length/2))
      # add vertical connector
      l.data <- data.frame (x=c(x2, x2), y=ys, edge)
      p.env$p.data <- rbind (p.env$p.data, l.data)
      # move on to next edge
      .cpMkLine (x2, ys, next.edge, p.env)
    }
  } else {
    stop ('Invalid tree')
  }
}
.cpGetParentNode <- function (p.env) {
  # get nodes
  ql.node <- p.env$tree$edge[p.env$ql.edge, 2]
  sl.node <- p.env$tree$edge[p.env$sl.edge, 2]
  # find ascending nodes
  ql.nodes <- p.env$tree.stats[[ql.node]][['ascend.nodes']]
  sl.nodes <- p.env$tree.stats[[sl.node]][['ascend.nodes']]
  # match to find highest shared
  maxi <- which.max (match (ql.nodes, sl.nodes))
  ql.nodes[maxi]
}
.cpShiftY <- function (p.env) {
  # find all edges descending from higher edge and add y.length/2
  # find all edges descending from lower edge and remove y.length/2
  # find parent vertical edge and add and remove y.length/2
  parent.node <- .cpGetParent (p.env)
  if (parent.node == p.env$root.node) {
    parent.edge <- 0
  } else {
    parent.edge <- which (p.env$tree$edge[ ,2] == parent.node)
  }
  # get next edges and their ys
  edges <- p.env$tree.stats[[parent.node]][['next.edges']]
  edge.1 <- p.data[p.data$edge == edges[1], ]
  edge.2 <- p.data[p.data$edge == edges[2], ]
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
  p.data[p.data$edge %in% higher.edges, 'y'] <-
    p.data[p.data$edge %in% higher.edges, 'y'] + y.length/2
  # lower
  lower.node <- p.env$tree$edge[lower, 2]
  lower.edges <- c (p.env$tree.stats[[lower.node]][['descend.edges']],
                    lower)
  p.data[p.data$edge %in% lower.edges, 'y'] <-
    p.data[p.data$edge %in% lower.edges, 'y'] - y.length/2
  # vertical
  vertical.bool <- which (p.data$edge == parent.edge)
  vertical.bool <- vertical.bool[(length (vertical.bool) - 1):
                                   length (vertical.bool)]
  p.data[vertical.bool[2], 'y'] <- p.data[vertical.bool[2], 'y'] + y.length/2
  p.data[vertical.bool[1], 'y'] <- p.data[vertical.bool[1], 'y'] - y.length/2
  p.data
}