# TODO:
# -- convert these to vectorized functions, recursive limits the tree size
# -- rethink age, no need to calculate it on the fly, age is immutable unless node spans are extended

.Node_preDist <- function (tree, node, res=0) {
  i <- 1
  for (n in node$prenode) {
    res[i] <- res[i] + node$span
    res[i] <- .Node_preDist (tree, tree@nodelist[[n]], res[i])
    i <- 1 + length (res)
    res[i] <- 0
  }
  max (res)
}

.Node_postDist <- function (tree, node, res=0) {
  i <- 1
  for (n in node$postnode) {
    res <- c (res[-i], .Node_postDist (tree, tree@nodelist[[n]], res[i]))
    i <- 1 + length (res)
    res[i] <- 0
  }
  max (res + node$span)
}

.getAge__NodeList <- function (tree, node) {
  tree_age <- .Node_postDist (tree, tree@nodelist[[tree@root]])
  tree_age - .Node_preDist (tree, tree@nodelist[[node]])
}


