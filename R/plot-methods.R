#' @name blockplot
#' @title Plot Blocks for Detecting Phylogenetic Signal
#' @description A plot function for visualing phylogenetic signal in
#'  a tree using blocks to represent branches.
#' @details No details
#' @template base_template
#' @param trait vector or matrix of traits for each tip
#' @param title title of plot (default NULL)
#' @export
#' @examples
#' # Generate a balanced tree
#' tree <- compute.brlen (stree (16, 'balanced'))
#' # Generate perfect phylogenetic signal traits
#' # (Note: trait vector should be named)
#' trait <- c (rep (1, 8), rep (0, 8))
#' names (trait) <- tree$tip.label
#' # Plot blocks ...
#' blockplot (tree, trait, 'Perfect signal')
#' # Repeat with evenly distributed traits ...
#' trait <- rep (c (0, 1), 8)
#' names (trait) <- tree$tip.label
#' blockplot (tree, trait, 'Maximally overdispered')
#' # The more evenly distributed, the more symetrical
#' # the pattern
#' # A random distribution of traits ...
#' trait <- randCommData (tree, nsites = 1, nspp = 8)[1, ]
#' blockplot (tree, trait, 'Random trait')
#' # ... isn't as symetrical

blockplot <- function (tree, trait, title = NULL) {
  # get tips descending from each edge (or more precisely the
  #  terminal node of each edge)
  children.by.edge <- mlply (.data = data.frame (
    node = tree$edge[ ,2]), .fun = getChildren, tree)
  .byEdge <- function (i, trait.type.names) {
    sum (trait.type.names %in% children.by.edge[[i]])
  }
  .byTraitType <- function (i) {
    trait.type.names <- names (trait)[trait ==
                                        trait.types[i]]
    res <- mdply (.data = data.frame (i = 1:length (
      children.by.edge)), .fun = .byEdge, trait.type.names)[ ,2]
    res * tree$edge.length
  }
  # get proportion of edge.length of each trait based on FP
  trait.types <- unique (trait)
  edge.mappings <- mdply (.data = data.frame (
    i = 1:length (trait.types)), .fun = .byTraitType)[ ,-1]
  colnames (edge.mappings) <- tree$edge[, 2]
  .byTip <- function (x) {
    tip <- tree$tip.label[x]
    edges <- getEdges (tree, tips = tip, type = 2)
    temp <- edge.mappings[tree$edge[, 2] %in%
                            tree$edge[edges, 2]]
    .eachTip <- function (i) {
      ord <- order (temp[ ,i], decreasing = TRUE)
      data.frame (trait = rownames (temp)[ord],
                  p.trait = temp[ord, i], tip,
                  edges = rep (edges, each = 2))
    }
    res <- mdply (.data = data.frame (
      i = 1:length (temp)), .fun = .eachTip)[ ,-1]
    res <- res[res$p.trait != 0, ]
    res$y2 <- cumsum (res$p.trait)
    res$y1 <- c (0, res$y2[-nrow (res)])
    res
  }
  # for each tip, find all its edges, stick together
  #  into a single data.frame for geom_rect ()
  # geom_rect () requires xy coordinates
  res <- mdply (.data = data.frame (
    x = 1:getSize (tree)), .fun = .byTip)
  p <- ggplot (res, aes (xmin = y1, xmax = y2,
                        ymin = x, ymax = x + 1,
                        fill = trait)) +
    geom_rect () + theme (line = element_blank (),
                         line = element_blank (),
                         axis.text = element_blank (),
                         panel.background = element_blank (),
                         legend.position = 'none')
  if (!is.null (title)) {
    print (p + ggtitle (title))
  } else {
    print (p)
  }
}