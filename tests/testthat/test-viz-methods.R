# Test viz-methods
# D.J. Bennett
# 26/10/2015

# USEFUL
#plot (p.env$tree);nodelabels();edgelabels();

# LIBS
library (MoreTreeTools)
library (testthat)

# PARAMETERS
ntips <- 10  # number of tips in the tree

# FUNCTIONS
init.1 <- function () {
  # push random tree and key parameters
  tree <- rtree (ntips)
  N <- length (tree$tip.label)
  tree.age <- getSize (tree, 'rtt')
  x.spacer <- tree.age/100
  y.spacer <<- 1
  x1 <<- tree.age
  node <<- length (tree$tip.label) + 1
  edge <<- which (tree$edge[ ,1] == node)
  # init p.env
  p.env <<- new.env (parent=emptyenv ())
  p.env$N <- N
  p.env$all.edges <- 1:nrow (tree$edge)
  p.env$p.data <- data.frame (x=c(0, 0), y=c(0, y.spacer),
                              edge=c (0,0))
  p.env$y.spacer <- y.spacer
  p.env$x.spacer <- x.spacer
  p.env$tree <- tree
  p.env$root.node <- length (tree$tip.label) + 1
  p.env$tree.stats <- getTreeStats (tree)
}

init.2 <- function () {
  # run subsequent init for advanced function testing
  MoreTreeTools:::.cpMkPData (x1, c (0, y.spacer), edge, p.env)
  p.env$p.data$line <- rep (1:(nrow (p.env$p.data)/2), each=2)
  p.env$lines <- unique (p.env$p.data$line)[-1]  # ignore root
  edges <- p.env$p.data$edge[match (p.env$lines, p.env$p.data$line)]  
  # get max and min
  p.env$maxmin <- data.frame (line=p.env$lines, edge=edges,
                              max.x=NA, min.x=NA,
                              max.y=NA, min.y=NA)
  MoreTreeTools:::.cpGetMaxMin (p.env$lines, p.env)
}


# RUNNING
context ('Testing \'viz-methods\'')
test_that ('.cpMkPData([basic]) works', {
  # init
  init.1 ()
  # run
  MoreTreeTools:::.cpMkPData (x1, c (0, y.spacer), edge, p.env)
  # sanity check
  res <- max (p.env$p.data$edge) == nrow (p.env$tree$edge)
  expect_that (res, is_true ())
  res <- max (p.env$p.data$x) == getSize (p.env$tree, 'rtt')
  expect_that (res, is_true ())
})
test_that ('.cpGetMaxMin([basic]) works', {
  # init
  init.1 ()
  init.2 ()
  # sanity checks
  expect_that (x1 < max (p.env$maxmin$max.x), is_true ())
  expect_that (nrow (p.env$maxmin), equals (length (p.env$lines)))
})
test_that ('.cpCheckOverlap([basic]) works', {
  # 1 spacer in x and y (lines must be without 0.5 from each other)
  # parallel lines (overlap is false)
  sbjct <- data.frame (max.x=1, min.x=0,
                       max.y=1, min.y=1)
  parallel <- data.frame (max.x=1, min.x=0,
                          max.y=2, min.y=2)
  parallel.res <- MoreTreeTools:::.cpCheckOverlap (sbjct, parallel)
  expect_that (parallel.res, is_false ())
  # crossing lines (overlap is true)
  crossing <- data.frame (max.x=0.5, min.x=0.5,
                          max.y=1.5, min.y=0.5)
  crossing.res <- MoreTreeTools:::.cpCheckOverlap (sbjct, crossing)
  expect_that (crossing.res, is_true ())
  # successive lines (overlap is true)
  # (this scenario would usually be ignored as this is expected
  #  if lines are descending)
  successive <- data.frame (max.x=1.5, min.x=0.5,
                            max.y=1, min.y=1)
  successive.res <- MoreTreeTools:::.cpCheckOverlap (sbjct, successive)
  expect_that (successive.res, is_true ())
  # all together
  all.together <- rbind (parallel, crossing, successive)
  all.together.res <- MoreTreeTools:::.cpCheckOverlap (sbjct, all.together)
  expect_that (all.together.res, equals (c (FALSE, TRUE, TRUE)))
})
test_that ('.cpGetParentNode([basic]) works', {
  # init
  init.1 ()
  # get random edges and make sure they match with list obj
  edges <- 1:nrow (p.env$tree$edge)
  edge.1 <- sample (edges, 1)
  edge.2 <- sample (edges, 1)
  p.node <- getParent (p.env$tree, edges=c (edge.1, edge.2))
  res <- MoreTreeTools:::.cpGetParentNode (edge.1, edge.2, p.env)
  expect_that (res, equals (p.node))
})
test_that ('.cpCorrectOverlap([basic]) works', {
  # mini function to prevent rounding errors
  test <- function (after, before) {
    signif (abs (after[1] - before[1]), 1) == spacer |
      signif (abs (after[1] + before[1])) == spacer
  }
  # init
  init.1 ()
  init.2 ()
  # get random connected edges
  internal.nodes <- (ntips + 1):((ntips*2)-1)
  rand.node <- sample (internal.nodes, 1)
  edge.1 <- p.env$tree.stats[[rand.node]][['next.edges']][1]
  edge.2 <- p.env$tree.stats[[rand.node]][['next.edges']][2]
  yi.1 <- which (p.env$p.data$edge == edge.1)
  yi.2 <- which (p.env$p.data$edge == edge.2)
  lines.to.change <- c (p.env$p.data$line[yi.1][1],
                        p.env$p.data$line[yi.2][1])
  # record before
  y.before.1 <- p.env$p.data[yi.1, 'y']
  y.before.2 <- p.env$p.data[yi.2, 'y']
  # run shift Y
  changed.lines <- MoreTreeTools:::.cpCorrectOverlap (edge.1, edge.2, p.env)
  # check Y shifted  -- should be equal to spacer/2
  spacer <- y.spacer/2
  y.after.1 <- p.env$p.data[yi.1, 'y']
  y.after.2 <- p.env$p.data[yi.2, 'y']
  shifted.1 <- test (y.after.1, y.before.1)
  shifted.2 <- test (y.after.2, y.before.2)
  expect_that (shifted.1, is_true ())
  expect_that (shifted.2, is_true ())
  expect_that (all (lines.to.change %in% changed.lines), is_true ())
})
test_that ('.cpCheckLine([basic]) works', {
  # init
  init.1 ()
  p.env$y.spacer <- 0
  init.2 ()
  # larger y.spacer, to force overlap
  p.env$y.spacer <- 1
  # random line
  rand.line <- sample (p.env$lines, 1)
  MoreTreeTools:::.cpCheckLine (rand.line, p.env)
  # make sure y has changed at least n.desc.edge times
  which.line <- p.env$p.data$line == rand.line
  line <- p.env$p.data[which.line, ]
  # find connecting edges
  edge <- line$edge[1]
  node <- p.env$tree$edge[edge,2]
  n.desc.edges <- length (p.env$tree.stats[[node]][['descend.edges']])
  expect_that (sum (abs (p.env$p.data$y)) , is_more_than (n.desc.edges))
})
test_that ('chromatophylo([basic]) works', {
  tree <- rtree (ntips)
  p <- chromatophylo (tree, reduce.overlap=TRUE)
  expect_that (class (p)[1], equals ('gg'))
})