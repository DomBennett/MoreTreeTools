## D.J Bennett
## Testing out new ideas

blockPlot <- function (tree, trait, title = NULL) {
  children.by.edge <- mlply (.data = data.frame (node = tree$edge[ ,2]),
                             .fun = getChildren, tree)
  .byEdge <- function (i, trait.type.names) {
    sum (trait.type.names %in% children.by.edge[[i]]) /
      length (children.by.edge[[i]])
  }
  .byTraitType <- function (i) {
    trait.type.names <- names (trait)[trait == trait.types[i]]
    res <- mdply (.data = data.frame (i = 1:length (children.by.edge)),
                  .fun = .byEdge, trait.type.names)[ ,2]
    res * tree$edge.length
  }
  trait.types <- unique (trait)
  edge.mappings <- mdply (.data = data.frame (i = 1:length (trait.types)),
                          .fun = .byTraitType)[ ,-1]
  colnames (edge.mappings) <- tree$edge[, 2]
  
  .byTip <- function (x) {
    tip <- tree$tip.label[x]
    edges <- getEdges (tree, tips = tip, type = 2)
    temp <- edge.mappings[tree$edge[, 2] %in%
                            tree$edge[edges, 2]]
    .eachTip <- function (i) {
      ord <- order (temp[ ,i], decreasing = TRUE)
      data.frame (trait = rownames (temp)[ord],
                  p.trait = temp[ord, i], tip)
    }
    res <- mdply (.data = data.frame (i = 1:length (temp)),
                  .fun = .eachTip)[ ,-1]
    res <- res[res$p.trait != 0, ]
    res$y2 <- cumsum (res$p.trait)
    res$y1 <- c (0, res$y2[-nrow (res)])
    res
  }
  res <- mdply (.data = data.frame (x = 1:getSize (tree)),
                .fun = .byTip)
  p <- ggplot (res, aes(xmin = x, xmax = x + 1,
                   ymin = y1, ymax = y2,
                   fill = trait)) +
    geom_rect () + theme(line = element_blank(),
                         line = element_blank(),
                         axis.text = element_blank(),
                         panel.background = element_blank(),
                         legend.position = 'none')
  if (!is.null (title)) {
    print (p + ggtitle (title))
  } else {
    print (p)
  }
}

## Plotting phylogenetic signal
tree <- compute.brlen (stree (128, 'balanced'))
rtrait <- randCommData (tree, nsites = 1, nspp = 64)[1, ]
blockPlot (tree, rtrait, 'Random')
etrait <- evenCommData (tree, nsites = 1, nspp = 63)[1, ]
blockPlot (tree, etrait, 'Even')
ctrait1 <- genCommData (tree, focal = sample (1:getSize (tree), 1),
                      psi = 100, mean.incid = 5)[1, ]
ctrait2 <- genCommData (tree, focal = sample (1:getSize (tree), 1),
                       psi = 100, mean.incid = 5)[1, ]
ctrait3 <- genCommData (tree, focal = sample (1:getSize (tree), 1),
                        psi = 100, mean.incid = 5)[1, ]
ctrait2 <- ifelse (ctrait2 == 1, 2, 0)
ctrait3 <- ifelse (ctrait3 == 1, 3, 0)
ctrait <- ctrait1 + ctrait2 + ctrait3
blockPlot (tree, ctrait, 'Clustered')
ctrait <- genCommData (tree, focal = sample (1:getSize (tree), 1),
                        psi = 100, mean.incid = 64)[1, ]
blockPlot (tree, ctrait, 'Clustered big')