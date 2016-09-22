

# map plot
library(MoreTreeTools)
library(ggplot2)



tree <- read.tree(file='/Users/djb208/Coding/Project-EDBMM/data/raw_trees/literature/bininda.tre')
plot(tree)
nodelabels()



addToCoords <- function(nd, p1, p2, p3) {
  coords$cnd <- coords$cnd + 1
  coords$nd[coords$cnd] <- nd
  coords$cnd <- coords$cnd + 1
  coords$nd[coords$cnd] <- nd
  coords$cnd <- coords$cnd + 1
  coords$nd[coords$cnd] <- nd
  coords$cx <- coords$cx + 1
  coords$x[coords$cx] <- p1[1]
  coords$cx <- coords$cx + 1
  coords$x[coords$cx] <- p2[1]
  coords$cx <- coords$cx + 1
  coords$x[coords$cx] <- p3[1]
  coords$cy <- coords$cy + 1
  coords$y[coords$cy] <- p1[2]
  coords$cy <- coords$cy + 1
  coords$y[coords$cy] <- p2[2]
  coords$cy <- coords$cy + 1
  coords$y[coords$cy] <- p3[2]
}
getNxtNds <- function(nd, pts) {
  nxt <- which(tree$edge[ ,1] == nd)
  if(length(nxt) == 2) {
    nds <- tree$edge[nxt, 2]
    h1 <- log(length(getEdges(tree, node=nd)) + 1)
    h2 <- log(length(getEdges(tree, node=nd)) + 1)
    if(h1 > h2) {
      t1 <- getTrngl(nds[1], pts[[1]], pts[[3]], h1)
      getNxtNds(nds[1], t1)
      t2 <- getTrngl(nds[2], pts[[3]], pts[[2]], h2)
      getNxtNds(nds[2], t2)
    } else {
      t1 <- getTrngl(nds[1], pts[[3]], pts[[2]], h1)
      getNxtNds(nds[1], t1)
      t2 <- getTrngl(nds[2], pts[[1]], pts[[3]], h2)
      getNxtNds(nds[2], t2)
    }
  }
}
getTrngl <- function(nd, p1, p2, h) {
  pbar <- (p1 + p2)/2
  nhat <- c(p1[2]-p2[2], p2[1]-p1[1])/sqrt((p1[2]-p2[2])^2+(p2[1]-p1[1])^2)
  p3 <- pbar + h*nhat
  addToCoords(nd, p1, p2, p3)
  list(p1, p2, p3)
}
# init
tree <- rtree(10)
tree <- compute.brlen(tree)
plot(tree, show.tip.label = FALSE)
# new coords env
coords <- new.env()
coords$nd <- coords$y <- coords$x <- vector(length=(nrow(tree$edge) + 1)*3)
coords$cnd <- coords$cy <- coords$cx <- 0
# ROOT
root <- length(tree$tip.label) + 1
h <- log(length(getEdges(tree, node=root)) + 1)
p1 <- c(0, 0)
p3 <- c((0+h)/2, h)
p2 <- c(0+h, 0)
addToCoords(nd, p1, p2, p3)
pts <- list(p1, p2, p3)
getNxtNds(root, pts)


ages <- getAge(tree)
positions <- data.frame(id=factor(coords$nd),
                        x=coords$x, y=coords$y,
                        age=ages$age[match(coords$nd, ages$node)])

p <- ggplot(positions, aes(x=x, y=y)) + geom_polygon(aes(fill=age, group=id))
p + scale_fill_gradient2() +
  theme(panel.background = element_rect(fill='#0000ff'))

