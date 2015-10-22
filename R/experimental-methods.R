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
  # internal
  mkline <- function (x1, y, edge) {
    if (length (edge) == 2) {
      # run recursively for both
      mkline (x1, y[1], edge[1])
      mkline (x1, y[2], edge[2])
    } else if (length (edge) == 1){
      # find the next node on the edge
      next.node <- tree$edge[edge,2]
      # get n.children
      n.children <- length (getChildren (tree, next.node))
      # get the length of the edge
      x2 <- x1 + tree$edge.length[edge]
      l.data <- data.frame (x=c(x1, x2), y=c(y, y), edge)
      p.data <<- rbind (p.data, l.data)
      # run again for next edge
      next.edge <- which (tree$edge[ ,1] == next.node)
      if (next.node <= length (tree$tip.label)) {
        #t.data <<- data.frame (x=(x2+1/N), y=y,
        #                       label=tree$tip.label[next.node])
      } else {
        ys <- c (y-(y.length/2), y+(y.length/2))
        # add vertical connector
        l.data <- data.frame (x=c(x2, x2), y=ys, edge)
        p.data <<- rbind (p.data, l.data)
        # move on to next edge
        mkline (x2, ys, next.edge)
      }
    } else {
      stop ('Invalid tree')
    }
  }
  if (!is.rooted (tree)) {
    stop ('Tree must be rooted')
  }
  if (!is.binary.tree (tree)) {
    stop ('Tree must be binary')
  }
  if (is.null (tree$edge.length)) {
    stop ('Tree must have branch lengths')
  }
  # init plot data
  N <- length (tree$tip.label)
  y.length <- 1/N
  x1 <- 0
  node <- length (tree$tip.label) + 1
  edge <- which (tree$edge[ ,1] == node)
  p.data <- data.frame (x=c(0, 0), y=c(0, y.length),
                        edge=c (0,0))
  #t.data <- data.frame (x=NA, y=NA, label=NA)
  mkline (x1, c (0, y.length), edge)
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

#' @name blockplot
#' @title Plot Blocks for Detecting Phylogenetic Signal (experimental)
#' @description Experimental plot function for visualing phylogenetic signal in
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
#' blockplot (tree, factor(trait), title = 'Perfect signal', gradient = FALSE)
#' # Or gradient ....
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

# TODO -- fix example
blockplot <- function (tree, trait, gradient = TRUE, title = NULL) {
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
    temp <- edge.mappings[,tree$edge[, 2] %in%
                            tree$edge[edges, 2], drop=FALSE]
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
  trait.types <- unique (trait)
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

#' @name flowplot
#' @title Plot tree as turbulence flow image (experimental)
#' @description Experimental plot where a tree is broken down into
#' X and Y pixels. Edges are mapped to this Cartesian space with
#' colour representing number of descendents. The idea is to better
#' visualise very large trees, particular ones with fossils. Currently only
#' works with ultrametric trees.
#' @details No details
#' @template base_template
#' @param n number of pixels
#' @export
#' @examples
#' # Generate and plot a random tree
#' tree <- compute.brlen (rtree (100))
#' flowplot (tree)

flowplot <- function (tree, n=getSize(tree)) {
  .getTimeSlice <- function (tree, xmin, xmax, ys) {
    .get <- function (edge) {
      # get Ds and Ys for all edges that pass through xmin and xmax
      # TODO: make this y calculation more general
      children <- getChildren (tree, tree$edge[edge, 2])
      res <- c (mean (children), length (children))
      names (res) <- c ('Y', 'D')
      res
    }
    .fill <- function (ymax) {
      # fill in y data frame
      vues <- which (res[ ,1] <= ymax & res[ ,1] > ymin)
      ymin <<- ymax
      data.frame ('x'=xmax, 'y'=ymax, 'd'=mean (res[vues, 2]))
    }
    edges <- which ((edge.ages[ ,1] < xmax | edge.ages[ ,2] < xmax) &
                      (edge.ages[ ,1] >= xmin | edge.ages[ ,2] >= xmin))
    res <- mdply (.data=data.frame (edge=edges), .fun=.get)[ ,-1]
    ymin <- 0
    mdply (.data=data.frame(ymax=ys), .fun=.fill)[ ,-1]
  }
  .getPixels <- function (tree, n) {
    xs <- seq (0, 1, length.out=n+1)
    ys <- seq (getSize (tree) / n, getSize (tree), length.out=n)
    pixels <- data.frame (x=-1, y=NA, d=NA)  # init df
    for (i in 1:n) {
      pixels <- rbind (pixels,
                       .getTimeSlice (tree, xmin=xs[i], xmax=xs[i+1], ys=ys))
    }
    pixels <- pixels[-1, ]
    pixels[is.na (pixels[, 'd']), 'd'] <- 0
    pixels
  }
  # convert tips to Y coords
  tree$tip.label <- 1:getSize (tree)
  # get edge ages
  node.ages <- getAge (tree)
  edge.ages <- tree$edge
  edge.ages[ ,1] <- node.ages[match (tree$edge[ ,1], node.ages[ ,1]), 2]
  edge.ages[ ,2] <- node.ages[match (tree$edge[ ,2], node.ages[ ,1]), 2]
  # get pixels
  pixels <- .getPixels (tree, n)
  # plot
  # TODO: workout how to get more colour range
  fill.grad <- scale_fill_continuous (limits=c(0, getSize (tree)), low='black',
                                    high='red')
  ggplot (pixels, aes (x=x, y=y)) +
    geom_tile (aes (fill=d)) + fill.grad +
    theme (line=element_blank (), line=element_blank (), axis.text=element_blank (),
           panel.background=element_blank (), legend.position="none")
}