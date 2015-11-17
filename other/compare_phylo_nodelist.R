# 16/11/2015
# Simple test to see whether NodeList object is faster for tree building than phylo

# ----
# Test
# ----
#
# 1. Grow tree from 3 tips to 100
# 2. Time how long it takes for phylo and Nodelist

# --------------
# Phylo function
# --------------
phyloBuilder <- function (tree, iterations) {
  tree$node.label <- paste0 ('n', (length (tree$tip.label)+1):
                               (length (tree$tip.label)+2))
  for (i in (length (tree$tip.label)+1):
         (length (tree$tip.label)+iterations)) {
    # add tip
    tree <- MoreTreeTools:::.addTip__phylo (tree, edge=1,
                                            tip_age=0, node_age=0.5,
                                            tip_name=paste0('t', i),
                                            node_name=paste0('n', i+2))
    # grow tip edges
    #tip_edges <- which (tree$edge[ ,2] %in% 1:length (tree$tip.label))
    #tree$edge.length[tip_edges] <- tree$edge.length[tip_edges] + 0.5
  }
  tree
}

# -----------------
# NodeList function
# -----------------
nodeListBuilder <- function (tree, iterations) {
  for (i in (nNodes (tree)+1):(nTips (tree)+iterations)) {
    # add tip
    tree <- MoreTreeTools:::.addTip__NodeList (tree, edge='t1',
                                               tip_age=0, node_age=0.5,
                                               tip_name=paste0('t', i),
                                               node_name=paste0('n', i+2))
  }
  tree
}

# -------
# Testing
# -------
iterations <- 1000
start_tree <- compute.brlen (stree (3, 'left'))
start_time <- Sys.time ()
tree <- phyloBuilder (start_tree, iterations)
end_time <- Sys.time ()
phylo_time <- end_time - start_time
rm (tree)
# start_tree <- as (start_tree, 'NodeList')
# start_time <- Sys.time ()
# tree <- nodeListBuilder (start_tree, iterations)
# end_time <- Sys.time ()
# nodelist_time <- end_time - start_time
start_time <- Sys.time ()
tree <- list ()
for (i in 1:iterations) {
  tree[[i]] <- i
}
end_time <- Sys.time ()
nodelist_time <- end_time - start_time
cat ('\nphylo time: ', phylo_time)
cat ('\nPseudo-NodeList time: ', nodelist_time)

# Really slow...
# Rprof(tmp <- tempfile())
# tree <- nodeListBuilder (start_tree, iterations)
# Rprof()
# summaryRprof(tmp)

# New thinks:
# The profile suggests that it is the creation of Node objects
#  that take time. I think I shold convert the Node objects
#  to S4 class or even just to a simple list. Either way it should be
#  lightweight. The major issue with the RefClass is it takes a copy
#  of all the methods associated with the object as well as the data.



