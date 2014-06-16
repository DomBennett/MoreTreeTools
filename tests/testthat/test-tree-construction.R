## Test tree construction (which btw is really hard...)
## D.J. Bennett
## 16/06/2014

## Libraries
library (MoreTreeTools)
library (testthat)

## Test data
data ('hominoids')

## Running tests
context ('Testing \'tree-construction\'')
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