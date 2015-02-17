#' @name blockplot
#' @title Plot Blocks for Detecting Phylogenetic Signal
#' @description A plot function for visualing phylogenetic signal in
#'  a tree using blocks to represent branches.
#' @details No details
#' @template base_template
#' @param trait vector or matrix of traits for each tip
#' @param title title of plot (default NULL)
#' @export
#' @examples
#' # Generate a balanced tree
#' tree <- compute.brlen (stree (16, 'balanced'))
#' # Generate perfect phylogenetic signal traits
#' # (Note: trait vector should be named)
#' trait <- c (rep (1, 8), rep (0, 8))
#' names (trait) <- tree$tip.label
#' # Plot blocks ...
#' blockplot (tree, trait, title = 'Perfect signal', gradient = TRUE)
#' # Repeat with evenly distributed traits ...
#' trait <- rep (c (0, 1), 8)
#' names (trait) <- tree$tip.label
#' blockplot (tree, trait, gradient = TRUE, title = 'Maximally overdispered')
#' # The more evenly distributed, the more symetrical
#' # the pattern
#' # A random distribution of traits ...
#' trait <- randCommData (tree, nsites = 1, nspp = 8)[1, ]
#' blockplot (tree, trait, title = 'Random trait', gradient = TRUE)
#' # ... isn't as symetrical

blockplot <- function (tree, trait, gradient = TRUE,
                       title = NULL) {
  .byTraitType <- function (i) {
    # calculate proportion of trait per edge for categorical
    #  traits
    .byEdge <- function (i, trait.type.names) {
      sum (trait.type.names %in% children.by.edge[[i]]) /
        length (children.by.edge[[i]])
    }
    trait.type.names <- names (trait)[
      trait == unique (trait)[i]]
    res <- mdply (.data = data.frame (i = 1:length (
      children.by.edge)), .fun = .byEdge,
      trait.type.names)[ ,2]
    res * tree$edge.length
  }
  .byTraitValue <- function (i) {
    # get mean value for each edge based on children
    mean (trait[names (trait) %in% 
                  children.by.edge[[i]]]) * tree$edge.length[i]
  }
  .byTip <- function (i) {
    # calculate blocks for each tip
    .eachTip <- function (j) {
      ord <- order (temp[ ,j], decreasing = TRUE)
      data.frame (trait = rownames (temp)[ord],
                  p.trait = temp[ord, j], tip,
                  edges = edges[j])
    }
    tip <- tree$tip.label[i]
    edges <- getEdges (tree, tips = tip, type = 2)
    temp <- edge.mappings[tree$edge[, 2] %in%
                            tree$edge[edges, 2]]
    res <- mdply (.data = data.frame (
      j = 1:length (temp)), .fun = .eachTip)[ ,-1]
    res$x2 <- cumsum (res$p.trait)
    res$x1 <- c (0, res$x2[-nrow (res)])
    res$y1 <- y1[i]
    res$y2 <- y2[i]
    res
  }
  .forGradient <- function (i) {
    # merge separate blocks of different traits, into one
    .byEdge <- function (j) {
      edge.data <- tip.data[tip.data$edges ==
                              unique (tip.data$edges)[j], ]
      # p.trait is now the proportion of what trait it reps
      #  between -1 and +1
      p.trait <- (edge.data$p.trait[edge.data$trait == 1] -
        edge.data$p.trait[edge.data$trait == 2]) /
        sum (edge.data$p.trait)
      x1 <- edge.data$x1[1]
      x2 <- edge.data$x2[2]
      y1 <- edge.data$y1[1]
      y2 <- edge.data$y2[2]
      edges <- edge.data$edges[1]
      tip <- edge.data$tip[1]
      data.frame (p.trait, edges, y1, y2, x1, x2)
    }
    tip.data <- rectcoords[rectcoords$tip == unique (
      rectcoords$tip)[i], ]
    mdply (.data = data.frame (j = 1:length (
      unique (tip.data$edges))), .fun = .byEdge)[ ,-1]
  }
  # 1/calc FP -- needed to determine width of bar for each tip
  # low FP, means larger width
  y2 <- cumsum (1/calcED (tree)[ ,1])
  y1 <- c (0, y2[-length (y2)])
  # get tips descending from each edge (or more precisely the
  #  terminal node of each edge)
  children.by.edge <- mlply (.data = data.frame (
    node = tree$edge[ ,2]), .fun = getChildren, tree)
  # get proportion of edge.length of each trait based on FP
  if (is.factor (trait)) {
    edge.mappings <- mdply (.data = data.frame (
      i = 1:length (trait.types)), .fun = .byTraitType)[ ,-1]
  } else {
    edge.mappings <- t (mdply (.data = data.frame (
      i = 1:nrow (tree$edge)), .fun = .byTraitValue)[ ,-1, FALSE])
  }
  colnames (edge.mappings) <- tree$edge[, 2]
  # for each tip, find all its edges, stick together
  #  into a single data.frame for geom_rect ()
  # geom_rect () requires xy coordinates
  rectcoords <- mdply (.data = data.frame (
    i = 1:getSize (tree)), .fun = .byTip)
  if (gradient & length (trait.types) == 2) {
    # if gradient, merge p.trait for each edge for
    #  each tip
    rectcoords <- mdply (.data = data.frame (
      i = 1:length (unique (rectcoords$tip))),
      .fun = .forGradient)[ ,-1]
    p <- ggplot (rectcoords, aes (xmin = x1, xmax = x2,
                                  ymin = y1, ymax = y2)) +
                   geom_rect (aes (fill = p.trait)) +
      scale_fill_gradient2 (low = "red", high = "blue")
  } else {
    rectcoords <- rectcoords[rectcoords$p.trait != 0, ]
    p <- ggplot (rectcoords, aes (xmin = x1, xmax = x2,
                                  ymin = y1, ymax = y2,
                                  fill = trait)) +
      geom_rect ()
  }
  # plot
  p <- p + theme (line = element_blank (),
                  line = element_blank (),
                  axis.text = element_blank (),
                  panel.background = element_blank ())
  if (!is.null (title)) {
    print (p + ggtitle (title))
  } else {
    print (p)
  }
}

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