# Test get-methods
# D.J. Bennett
# 05/05/2014

# LIBS
library (MoreTreeTools)
library (testthat)

# TEST DATA
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

# RUNNING
context ('Testing \'get-methods\'')
test_that ('.getNodeStats__phylo([basic]) works', {
  tree <- rtree (25)
  res <- MoreTreeTools:::.getNodeStats__phylo (tree)
  expect_that (nrow (res), equals (tree$Nnode))
  expect_that (max (res$n.children), equals (25))
})

test_that ('.getEdgeStats__phylo([basic]) works', {
  tree <- rtree (25)
  res <- MoreTreeTools:::.getEdgeStats__phylo (tree)
  expect_that (nrow (res) == nrow (tree$edge), is_true ())
  expect_that (max (res$n.children), is_less_than (25))
})

test_that ('.getOutgroup__phylo([basic]) works',{
  tips <- hominoids$tip.label[-1]
  tips <- c (sample (tips, 3), "Macaca mulatta")
  outgroup <- MoreTreeTools:::.getOutgroup__phylo (hominoids, tips)
  expect_true (outgroup == "Macaca mulatta")
  tips <- c ("Hylobates muelleri", "Hylobates pileatus",
             "Hylobates lar")
  outgroup <- MoreTreeTools:::.getOutgroup__phylo (hominoids, tips)
  expect_true (length (outgroup) == length (tips))
})

test_that ('.getExtant__phylo([basic]) works',{
  tree <- stree (20)
  tree$edge.length <- rep (10, 20)
  tree$edge.length[1] <- 1  # this will be t1
  res <- MoreTreeTools:::.getExtant__phylo (tree, tol=1e-8)
  expect_false ('t1' %in% res)
  expect_true ('t2' %in% res)
})

test_that ('.getChildren__phylo([basic]) works',{
  expect_that (MoreTreeTools:::.getChildren__phylo (hominoids, 20),
               is_equivalent_to (hominids))
  expect_that (MoreTreeTools:::.getChildren__phylo (hominoids, '20'),
               throws_error ())
  expect_that (MoreTreeTools:::.getChildren__phylo (hominoids, 31),
               throws_error ())
})

test_that ('.getSize__phylo([basic]) works', {
  expect_that (MoreTreeTools:::.getSize__phylo (hominoids),
               equals (length (hominoids$tip.label)))
  expect_that (MoreTreeTools:::.getSize__phylo (hominoids, type = 'pd'),
               equals (sum (hominoids$edge.length)))
  expect_that (MoreTreeTools:::.getSize__phylo (hominoids, type = 'rtt'),
               equals (max (diag (vcv.phylo (hominoids)))))
})

test_that ('.getAge__phylo([basic]) works', {
  root.node <- length (hominoids$tip.label) + 1
  root.age <- MoreTreeTools:::.getAge__phylo (hominoids, root.node)[ ,2]
  stem.age <- MoreTreeTools:::.getAge__phylo (hominoids, root.node + 1)[ ,2]
  # root age should ALWAYS be greater than stem age
  expect_that (root.age > stem.age, is_true ())
})

test_that ('.getSister__phylo([basic]) works', {
  # 21 is hominins, their sister is the Orangutan
  orangutan <- which (hominoids$tip.label == 'Pongo pygmaeus')
  expect_that (MoreTreeTools:::.getSister__phylo (hominoids, 21),
               is_equivalent_to (orangutan))
})

test_that ('.getParent__phylo([node specified]) works', {
  expect_that (MoreTreeTools:::.getParent__phylo (hominoids, node = 21),
               equals (20))
})

test_that ('.getParent__phylo([tips specified]) works', {
  tips <- c ('Gorilla gorilla', 'Pan troglodytes', 'Pongo pygmaeus')
  expect_that (MoreTreeTools:::.getParent__phylo (hominoids, tips=tips),
               equals (20))
})

test_that ('.getParent__phylo([edges specified]) works', {
  res <- MoreTreeTools:::.getParent__phylo (hominoids, edges=c (6, 25))
  expect_that (res, equals (19))
  # test with edge labels
  hominoids$edge.label <- paste0 ('edge_', 1:nrow (hominoids$edge))
  res <- MoreTreeTools:::.getParent__phylo (hominoids,
                                            edges=c ('edge_6', 'edge_25'))
  expect_that (res, equals (19))
})

test_that ('.getEdges__phylo([node specified]) works', {
  expect_that (MoreTreeTools:::.getEdges__phylo (hominoids, 21),
               is_equivalent_to (c (5, 6, 7, 8, 9, 10)))
})

test_that ('.getEdges__phylo([tips specified, type = 1]) works', {
  # a phylogeny containing only the tips
  expect_that (MoreTreeTools:::.getEdges__phylo (hominoids,
                                                 tips = hominins, type = 1),
               is_equivalent_to (c (7, 5, 10, 9, 8, 6)))
})

test_that ('.getEdges__phylo([tips specified, type = 2]) works', {
  # a phylogeny containing only the tips, plus the branches to the root
  expect_that (MoreTreeTools:::.getEdges__phylo (hominoids,
                                                 tips = hominins, type = 2),
               is_equivalent_to (c (7, 5, 10, 9, 8, 6, 4, 3, 2)))
})

test_that ('.getEdges__phylo([tips specified, type = 3]) works', {
  # the branches unique to the tips specified
  expect_that (MoreTreeTools:::.getEdges__phylo (hominoids,
                                                 tips = hominins, type = 3),
               is_equivalent_to (c (7, 5, 10, 9, 8, 6, 4)))
})

test_that ('.getNodes__phylo([basic]) works',{
  expect_that (MoreTreeTools:::.getNodes__phylo (hominoids, 21),
               is_equivalent_to (c (20, 19, 18)))
})

test_that ('.getClades__phylo([basic]) works',{
  # Use hominins -- a more manageable size
  expect_that (MoreTreeTools:::.getClades__phylo (hominins.tree),
               is_equivalent_to (hominins.clades))
})

test_that ('.getNodeLabels__phylo([basic]) works', {
  # Use hominins -- a lot faster
  expect_that (MoreTreeTools:::.getNodeLabels__phylo (hominins.tree),
               is_equivalent_to (c ("Homininae", "Homininae", "Pan")))
})

test_that ('.getSubtrees__phylo([basic]) works', {
  clade.trees <- MoreTreeTools:::.getSubtrees__phylo (hominoids, 5, 10)
  # should return two trees of 5 and 8
  expect_that (getSize (clade.trees[[1]]), equals (5))
  expect_that (getSize (clade.trees[[2]]), equals (8))
  expect_that (length (clade.trees), equals (2))
})

test_that ('.getTreeStats__phylo([basic]) works', {
  tree <- rtree (100)
  nodes <- 1:(length (tree$tip.label) + tree$Nnode)
  tree.stats <- MoreTreeTools:::.getTreeStats__phylo (tree)
  # check children
  rand.node <- sample (nodes, 1)
  children <- getChildren (tree, node=rand.node)
  expect_that (tree.stats[[rand.node]][['children']], equals (children))
  # check edges
  rand.node <- sample (nodes, 1)
  edges <- getEdges (tree, node=rand.node)
  bool <- all (tree.stats[[rand.node]][['descend.edges']] %in% edges &
                 edges %in% tree.stats[[rand.node]][['descend.edges']])
  expect_that (bool, is_true ())
  # check ascending nodes
  rand.node <- sample (nodes, 1)
  nodes <- getNodes (tree, node=rand.node)
  bool <- all (tree.stats[[rand.node]][['ascend.nodes']] %in% nodes &
                 nodes %in% tree.stats[[rand.node]][['ascend.nodes']])
  expect_that (bool, is_true ())
  # check age
  rand.node <- sample (nodes, 1)
  age <- getAge (tree, node=rand.node)[ ,2]
  expect_that (tree.stats[[rand.node]][['age']], equals (age))
})