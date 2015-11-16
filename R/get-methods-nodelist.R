# All hidden functions for get-methods-visible.R using NodeList class
.getTreeStats__NodeList <- function (tree) {
  cat ('\nNot implemented for NodeList')
}
.getNodeStats__NodeList <- function (tree, nodes='all', ignore.tips=NULL) {
  cat ('\nNot implemented for NodeList')
}
.getEdgeStats__NodeList <- function (tree, edges='all', ignore.tips=NULL) {
  cat ('\nNot implemented for NodeList')
}
.getOutgroup__NodeList <- function (tree, tips) {
  cat ('\nNot implemented for NodeList')
}
.getExtant__NodeList <- function (tree, tol) {
  tree <- setTol (tree, tol)
  extant (tree)
}
.getChildren__NodeList <- function (tree, node, display=FALSE) {
  node <- tree[[node]]
  node$children ()
}
.getSize__NodeList <- function (tree, type = c ('ntips', 'pd', 'rtt')) {
  cat ('\nNot implemented for NodeList')
}
.getAge__NodeList <- function (tree, node='all', edge=NULL) {
  cat ('\nNot implemented for NodeList')
}
.getSister__NodeList <- function (tree, node = 'all') {
  cat ('\nNot implemented for NodeList')
}
.getParent__NodeList <- function (tree, node=NULL, tips=NULL, edges=NULL) {
  prenodes <- vector ()
  nodes <- c (node, tips, edges)  # all must be node IDs
  if (!all (nodes %in% c (tree@nodes, tree@tips))) {
    stop ('Not all IDs in NodeList')
  }
  refprenodes <- tree[[nodes[1]]]$prenodes()
  min_rank <- 0
  for (n in nodes) {
    prenodes <- tree[[n]]$prenodes ()
    rank <- min (match (refprenodes, prenodes))
    if (rank > min_rank) min_rank <- rank
  }
  refprenodes[min_rank]
}
.getEdges__NodeList <- function (tree, node = NULL, tips = NULL, type = 1) {
  cat ('\nNot implemented for NodeList')
}
.getNodes__NodeList <- function (tree, node) {
  cat ('\nNot implemented for NodeList')
}
.getClades__NodeList <- function (tree) {
  cat ('\nNot implemented for NodeList')
}
.getNodeLabels__NodeList <- function (tree, all = FALSE, datasource = 4) {
  cat ('\nNot implemented for NodeList')
}
.getSubtrees__NodeList <- function (tree, min.n, max.n, verbose = FALSE) {
  cat ('\nNot implemented for NodeList')
}