# Test build-methods
# D.J. Bennett
# 16/06/2014

# LIBS
library (MoreTreeTools)
library (testthat)

# TEST DATA
data ('hominoids')

# RUNNING
context ('Testing \'build-methods\'')
test_that ('addTip([basic]) works', {
  # simply show that an errors are thrown and handled
  expect_that (addTip (hominoids, edge = 1, tip.age = 10,
                        node.age = 9, tip.name = 'new.tip'),
               throws_error ())
  expect_that (addTip (hominoids, edge = 1, tip.age = 10,
                        node.age = 40, tip.name = 'new.tip'),
               throws_error ())
  # actual test of functionality
  res <- addTip (hominoids, edge = 1, tip.age = 10,
          node.age = 20, tip.name = 'new.tip')
  expect_that (length (res$tip.label),
               equals (length (hominoids$tip.label) + 1))
})

test_that ('removeTip([preserve.age=FALSE]) works', {
  tree <- rtree (10)
  res <- removeTip (tree, tips=c ('t1', 't2'), preserve.age=FALSE)
  mat.before <- cophenetic.phylo (tree)
  mat.after <- cophenetic.phylo (res)
  bool <- mat.before['t3', c ('t4', 't5')] ==
    mat.after['t3', c ('t4', 't5')]
  expect_that (all (bool), is_true ())
})

test_that ('removeTip([preserve.age=TRUE]) works', {
  tree <- rtree (10)
  res <- removeTip (tree, tips=c ('t1', 't2'), preserve.age=TRUE)
  age.before <- getSize (tree, 'rtt')
  age.after <- getSize (res, 'rtt')
})

test_that ('collapseTips([basic]) works', {
  tree <- rtree (50)
  min.length <- getSize (tree, 'rtt')*0.1
  res <- collapseTips (tree, min.length, iterative=TRUE)
  tip.edges <- tree$edge[ ,2] < getSize (tree)
  bool <- all ((tree$edge.length < min.length) & tip.edges)
  expect_that (bool, is_false())
})