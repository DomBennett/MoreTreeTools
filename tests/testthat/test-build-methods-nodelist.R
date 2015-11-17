# Test build-methods
# D.J. Bennett
# 16/06/2014

# LIBS
library (MoreTreeTools)
library (testthat)

# TEST DATA
data ('hominoids')
hominoids <- as (hominoids, 'NodeList')

# RUNNING
context ('Testing \'build-methods\'')
test_that ('.addTip__NodeList([basic]) works', {
  # actual test of functionality
  res <- MoreTreeTools:::.addTip__NodeList (hominoids,
                                            edge='Macaca mulatta',
                                            tip_age=10,
                                            node_age=20,
                                            tip_name='new_tip',
                                            node_name='new_node')
  expect_that (nTips (res),
               equals (nTips (hominoids) + 1))
})

test_that ('.removeTip__NodeList([preserve.age=FALSE]) works', {
  
})

test_that ('.removeTip__NodeList([preserve.age=TRUE]) works', {
  
})

test_that ('.collapseTips__NodeList([basic]) works', {
  
})