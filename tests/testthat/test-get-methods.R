## Test get methods
## D.J. Bennett
## 05/05/2014

## Libraries
library (MoreTreeTools)
library (testthat)

## Test data
data ('hominoids')
hominids <- c ("Pongo pygmaeus", "Gorilla gorilla",
     "Homo sapiens", "Pan paniscus",
     "Pan troglodytes")
hominins <- c ('Homo sapiens', 'Gorilla gorilla',
               'Pan troglodytes', 'Pan paniscus')
hominins.tree <- drop.tip (hominoids, tip = hominoids$tip.label
                      [!hominoids$tip.label %in% hominins])
clade.node <- c (5, 6, 7, 1, 2, 3, 4)
clade.children <- list ('5' = c ("Gorilla gorilla", "Homo sapiens",
                                 "Pan paniscus", "Pan troglodytes"),
                        '6' = c ("Homo sapiens", "Pan paniscus",
                                 "Pan troglodytes"),
                        '7' = c ("Pan paniscus", "Pan troglodytes"),
                        '1' = 'Gorilla gorilla', '2' = 'Homo sapiens',
                        '3' = 'Pan paniscus', '4' = 'Pan troglodytes')
hominins.clades <- list (clade.children = clade.children,
                         clade.node = clade.node)

## Running tests
context ('Testing \'get-methods\'')
test_that ('getChildren([basic]) works',{
  expect_that (getChildren (hominoids, 20),
               is_equivalent_to (hominids))
  expect_that (getChildren (hominoids, '20'),
               throws_error ())
  expect_that (getChildren (hominoids, 31),
               throws_error ())
})

test_that ('getAge([basic]) works', {
  root.node <- length (hominoids$tip.label) + 1
  root.age <- getAge (hominoids, root.node)
  stem.age <- getAge (hominoids, root.node + 1)
  # root age should ALWAYS be greater than stem age
  expect_that (root.age > stem.age, is_true ())
})

test_that ('getSister([basic]) works', {
  # 21 is hominins, their sister is the Orangutan
  orangutan <- which (hominoids$tip.label == 'Pongo pygmaeus')
  expect_that (getSister (hominoids, 21), is_equivalent_to (orangutan))
})

test_that ('getEdges([node specified]) works', {
  expect_that (getEdges (hominoids, 21),
               is_equivalent_to (c (5, 6, 7, 8, 9, 10)))
})

test_that ('getEdges([tips specified, type = 1]) works', {
  # a phylogeny containing only the tips
  expect_that (getEdges (hominoids, tips = hominins, type = 1),
               is_equivalent_to (c (7, 5, 10, 9, 8, 6)))
})

test_that ('getEdges([tips specified, type = 2]) works', {
  # a phylogeny containing only the tips, plus the branches to the root
  expect_that (getEdges (hominoids, tips = hominins, type = 2),
               is_equivalent_to (c (7, 5, 10, 9, 8, 6, 4, 3, 2)))
})

test_that ('getEdges([tips specified, type = 3]) works', {
  # the branches unique to the tips specified
  expect_that (getEdges (hominoids, tips = hominins, type = 3),
               is_equivalent_to (c (7, 5, 10, 9, 8, 6, 4)))
})

test_that ('getNodes([basic]) works',{
  expect_that (getNodes (hominoids, 21),
               is_equivalent_to (c (20, 19, 18)))
})

test_that ('getClades([basic]) works',{
  # Use hominins -- a more manageable size
  expect_that (getClades (hominins.tree),
               is_equivalent_to (hominins.clades))
})

test_that ('getNodeLabels([basic]) works', {
  # Use hominins -- a lot faster
  expect_that (getNodeLabels (hominins.tree),
               is_equivalent_to (c ("Homininae", "Homininae", "Pan")))
})