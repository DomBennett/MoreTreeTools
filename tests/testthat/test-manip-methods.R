# Test manip-methods
# D.J. Bennett
# 19/05/2015

# LIBS
library (MoreTreeTools)
library (testthat)

# RUNNING
context ('Testing \'manip-methods\'')
test_that ('reLoad([basic]) works ....', {
  tree <- compute.brlen (rtree (50))
  tree.uploaded <- reLoad (tree, factor=1)
  gamma.up <- gammaStat (tree.uploaded)
  tree.downloaded <- reLoad (tree, factor=-1)
  gamma.down <- gammaStat (tree.downloaded)
  expect_more_than (gamma.up, gamma.down)
  # TODO -- fix reLoad so it can work with any tree
})
# test_that ('reBalance([basic]) works ....', {
#   TODO
# })