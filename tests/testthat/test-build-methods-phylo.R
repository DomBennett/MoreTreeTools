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
test_that ('.addTip__phylo([basic]) works', {
  res <- MoreTreeTools:::.addTip__phylo (hominoids, edge=1,
                                         tip_age=10, node_age=20,
                                         tip_name='new_tip')
  expect_that (length (res$tip.label),
               equals (length (hominoids$tip.label) + 1))
})

test_that ('.removeTip__phylo([preserve.age=FALSE]) works', {
  tree <- rtree (5)
  res <- MoreTreeTools:::.removeTip__phylo (tree, tips=c ('t1', 't2'),
                                            preserve_age=FALSE)
  mat_before <- cophenetic.phylo (tree)
  mat_after <- cophenetic.phylo (res)
  bool <- mat_before['t3', c ('t4', 't5')] ==
    mat_after['t3', c ('t4', 't5')]
  expect_that (all (bool), is_true ())
})

test_that ('.removeTip__phylo([preserve.age=TRUE]) works', {
  # TODO: occasionally the root.edge is not added, explore why this might fail
  tree <- rtree (5)
  res <- MoreTreeTools:::.removeTip__phylo (tree, tips=c ('t1', 't2'),
                                            preserve_age=TRUE)
  age_before <- getSize (tree, 'rtt')
  age_after <- getSize (res, 'rtt')
  expect_that (age_before, equals (age_after))
})

test_that ('.collapseTips__phylo([basic]) works', {
  tree <- rtree (50)
  min_length <- getSize (tree, 'rtt')*0.1
  res <- MoreTreeTools:::.collapseTips__phylo (tree,
                                               min_length,
                                               iterative=TRUE)
  tip_edges <- tree$edge[ ,2] < getSize (tree)
  bool <- all ((tree$edge.length < min_length) & tip_edges)
  expect_that (bool, is_false())
})