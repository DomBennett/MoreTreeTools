# Test viz-methods
# D.J. Bennett
# 26/10/2015

# LIBS
library (MoreTreeTools)
library (testthat)

# FUNCTIONS
init <- function () {
  # push random tree and key parameters
  tree <<- rtree (10)
  N <<- length (tree$tip.label)
  tree.age <- getSize (tree, 'rtt')
  x.length <<- tree.age/100
  y.length <<- 1
  x1 <<- 0
  node <<- length (tree$tip.label) + 1
  edge <<- which (tree$edge[ ,1] == node)
  # init p.env
  p.env <<- new.env (parent=emptyenv ())
  p.env$p.data <- data.frame (x=c(0, 0), y=c(0, y.length),
                              edge=c (0,0))
  p.env$y.length <- y.length
  p.env$x.length <- x.length
  p.env$tree <- tree
  p.env$root.node <- length (tree$tip.label) + 1
}


# RUNNING
context ('Testing \'viz-methods\'')
test_that ('.cpMkLine([basic]) works', {
  # init
  init ()
  # run
  .cpMkLine (x1, c (0, y.length), edge, p.env)
  # sanity check
  res <- max (p.env$p.data$edge) == nrow (tree$edge)
  expect_that (res, is_true ())
  res <- max (p.env$p.data$x) == getSize (tree, 'rtt')
  expect_that (res, is_true ())
})

test_that ('.cpCheckOverlap([basic]) works', {
  # 1 spacer in x and y (lines must be without 0.5 from each other)
  x.length <- y.length <- 1
  # parallel lines (overlap is false)
  sbjct <- data.frame (x = c (0, 1),
                       y = c (1, 1))
  qry <- data.frame (x = c (0, 1),
                     y = c (2, 2))
  res <- .cpCheckOverlap (sbjct, qry, y.length, x.length)
  expect_that (res, is_false ())
  # crossing lines (overlap is true)
  sbjct <- data.frame (x = c (0, 1),
                       y = c (1, 1))
  qry <- data.frame (x = c (0.5, 0.5),
                     y = c (0.5, 1.5))
  res <- .cpCheckOverlap (sbjct, qry, y.length, x.length)
  expect_that (res, is_true ())
  # parallel lines, too close (overlap is true)
  sbjct <- data.frame (x = c (0, 1),
                       y = c (1, 1))
  qry <- data.frame (x = c (0, 1),
                     y = c (1.49, 1.49))
  res <- .cpCheckOverlap (sbjct, qry, y.length, x.length)
  expect_that (res, is_true ())
  # parallel lines, at the right distance (overlap is false)
  sbjct <- data.frame (x = c (0, 1),
                       y = c (1, 1))
  qry <- data.frame (x = c (0, 1),
                     y = c (1.51, 1.51))
  res <- .cpCheckOverlap (sbjct, qry, y.length, x.length)
  expect_that (res, is_false ())
  # successive lines (overlap is true)
  # (this scenario would usually be ignored as this is expected
  #  if lines are descending)
  sbjct <- data.frame (x = c (0, 1),
                       y = c (1, 1))
  qry <- data.frame (x = c (1, 2),
                     y = c (1, 1))
  res <- .cpCheckOverlap (sbjct, qry, y.length, x.length)
  expect_that (res, is_true ())
})

test_that ('.cpEachQry([basic]) works', {
  # init
  init ()
  .cpMkLine (x1, c (0, y.length), edge, p.env)
  # run 1  -- edge is avoided
  ql <- 0
  sl.edge <- 1
  .cpEachQry (ql, p.env, avoid.edges, sbjct, sl.edge)
  expect_that (p.env$p.data, equals (p.data))
  # run 2 -- edge is not avoided
  avoid.edges <- c (4)
  .cpEachQry (ql, p.env, avoid.edges, sbjct, sl.edge)
  expect_that (p.env$p.data, equals (p.data))
  # run 3  -- edge no longer overlaps
  .cpEachQry (ql, p.env, avoid.edges, sbjct, sl.edge)
  expect_that (p.env$p.data, equals (p.data))
})

test_that ('.cpGetNChildren([basic]) works', {
  # init
  init ()
  # sanity checks
  .cpGetNChildren (p.env)
  expect_that (p.env$n.children[[11]], equals (10))
  expect_that (p.env$n.children[[1]], equals (1))
})