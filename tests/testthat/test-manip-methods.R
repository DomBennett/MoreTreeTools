# Test manip-methods
# D.J. Bennett
# 19/05/2015

# LIBS
library (MoreTreeTools)
library (testthat)

# RUNNING
context ('Testing \'manip-methods\'')
test_that ('reGravitise([basic]) works ....', {
  tree <- compute.brlen (rtree (100))
  tree.uploaded <- reGravitise (tree, factor=1)
  gamma.up <- gammaStat (tree.uploaded)
  tree.downloaded <- reGravitise (tree, factor=-1)
  gamma.down <- gammaStat (tree.downloaded)
  expect_more_than (gamma.up, gamma.down)
})