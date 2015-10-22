

library (MoreTreeTools)

## Function
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
    lines <- seq (from=3, to=nrow (p.env$p.data), by=2)
    m_ply (.data=data.frame (sl=lines), .fun=.cpEachSbjct, lines=lines,
           p.env=p.env)
    #return (p.env$p.data)
    print ('here')
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
.cpEachQry <- function (ql, p.env, avoid.edges, sbjct, sl.edge) {
  # unpack p.env
  p.data <- p.env$p.data
  y.length <- p.env$y.length
  x.length <- p.env$x.length
  # get qry
  qry <- p.data[ql:(ql+1), ]
  ql.edge <- unique (qry$edge)
  if (!ql.edge %in% avoid.edges) {
    overlap <- .cpCheckOverlap (sbjct, qry, y.length, x.length)
    if (overlap) {
      # set overlap to TRUE
      p.env$overlap <- TRUE
      # find all edges descdning from higher edge and add y.length/2
      # find all edges descending from lower edge and remove y.length/2
      # find parent vertical edge and add and remove y.length/2
      parent.node <- getParent (p.env$tree, edges=c (sl.edge, ql.edge))
      if (parent.node == p.env$root.node) {
        parent.edge <- 0
      } else {
        parent.edge <- which (p.env$tree$edge[ ,2] == parent.node)
      }
      edges <- which (p.env$tree$edge[ ,1] == parent.node)
      edge.1 <- p.data[p.data$edge == edges[1], ]
      edge.2 <- p.data[p.data$edge == edges[2], ]
      if (mean (edge.1$y) >= mean (edge.2$y)) {
        higher <- edge.1$edge[1]
        lower <- edge.2$edge[1]
      } else {
        higher <- edge.2$edge[1]
        lower <- edge.1$edge[1]
      }
      higher.edges <- c (getEdges (p.env$tree, node=p.env$tree$edge[higher, 2]),
                         higher)
      p.data[p.data$edge %in% higher.edges, 'y'] <-
        p.data[p.data$edge %in% higher.edges, 'y'] + y.length/2
      lower.edges <- c (getEdges (p.env$tree, node=p.env$tree$edge[lower, 2]),
                        lower)
      p.data[p.data$edge %in% lower.edges, 'y'] <-
        p.data[p.data$edge %in% lower.edges, 'y'] - y.length/2
      vertical.bool <- which (p.data$edge == parent.edge)
      vertical.bool <- vertical.bool[(length (vertical.bool) - 1):
                                       length (vertical.bool)]
      p.data[vertical.bool[2], 'y'] <- p.data[vertical.bool[2], 'y'] + y.length/2
      p.data[vertical.bool[1], 'y'] <- p.data[vertical.bool[1], 'y'] - y.length/2
      # update p.data in p.env
      p.env$p.data <- p.data
    }
  }
}
.cpCheckOverlap <- function (sbjct, qry, y.length, x.length) {
  # check if sbjct and qry straight lines overlap (within y.length)
  maxqx <- max (qry$x) - x.length/4
  minqx <- min (qry$x) + x.length/4
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
.cpEachSbjct <- function (sl, lines, p.env) {
  sbjct <- p.env$p.data[sl:(sl+1), ]
  # find connecting edges
  sl.edge <- unique (sbjct$edge)
  nodes <- p.env$tree$edge[sl.edge, ]
  avoid.edges <- c (which (p.env$tree$edge[ ,2] == nodes[1]),
                    which (p.env$tree$edge[ ,1] == nodes[2]),
                    sl.edge)
  m_ply (.data=data.frame (ql=lines), .cpEachQry, p.env=p.env,
         avoid.edges=avoid.edges, sbjct=sbjct, sl.edge=sl.edge)
}
.cpMkLine <- function (x1, y, edge, p.env) {
  if (length (edge) == 2) {
    # run recursively for both
    .cpMkLine (x1, y[1], edge[1], p.env)
    .cpMkLine (x1, y[2], edge[2], p.env)
  } else if (length (edge) == 1){
    # find the next node on the edge
    next.node <- p.env$tree$edge[edge,2]
    # get n.children
    n.children <- length (getChildren (p.env$tree, next.node))
    # get the length of the edge
    x2 <- x1 + p.env$tree$edge.length[edge]
    l.data <- data.frame (x=c(x1, x2), y=c(y, y), edge)
    p.env$p.data <- rbind (p.env$p.data, l.data)
    # run again for next edge
    next.edge <- which (p.env$tree$edge[ ,1] == next.node)
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

# explore
tree <- rtree (20)
Rprof(tmp <- tempfile())
chromatophylo (tree)
Rprof()
summaryRprof(tmp)
