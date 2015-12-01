# 16/11/2015
# Simple test to see whether NodeList object is faster for tree building than phylo

# ----
# Test
# ----
#
# 1. Grow tree from 3 tips to 100
# 2. Time how long it takes for phylo and Nodelist

# LIBS
library (MoreTreeTools)

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
  for (i in 4:(iterations+3)) {
    # add tip
    tree <- MoreTreeTools:::.addTip__NodeList (tree, edge='t1',
                               tip_age=0, node_age=0.5,
                               tip_name=paste0('t', i),
                               node_name=paste0('n', i+2))
  }
  tree
}

# ------------------------
# Pseudo-NodeList function
# ------------------------
pseudoNodeListBuilder <- function (iterations) {
  #setClass ('Node', representation=representation(i='numeric'))
  tree <- list ()
  for (i in 1:iterations) {
    # create a new S4 class object and add to list
    #node <- new ('Node', i=i)
    node <- list (i)
    tree[[i]] <- node
  }
}


# ---------
# Profiling
# ---------
# Really slow...
# Rprof(tmp <- tempfile())
# tree <- nodeListBuilder (start_tree, iterations)
# Rprof()
# summaryRprof(tmp)


# -------
# Testing
# -------
repeats <- seq (10, 1010, 500)
res <- data.frame (N_iterations=NA, Timing=NA, Build_method=NA)
start_phylo <- compute.brlen (stree (3, 'left'))
start_nodelist <- as (start_phylo, 'NodeList')
for (iterations in repeats) {
  cat ('\nphylo....', iterations)
  phylo_time <- system.time (expr = {
    tree <- phyloBuilder (start_phylo, iterations)
    })
  rm (tree)
  cat ('\nnodelist....', iterations)
  nodelist_time <- system.time (expr={
    tree <- nodeListBuilder (start_nodelist, iterations)
    })
  rm (tree)
  cat ('\npseudo-nodelist....', iterations)
  pnodelist_time <- system.time (expr={
    pseudoNodeListBuilder (iterations)
    })
  res <- rbind (res, data.frame (N_iterations=iterations, Timing=c (phylo_time[3], nodelist_time[3],
                                                                    pnodelist_time[3]),
                                 Build_method=c('phylo', 'nodelist', 'pseudo-nodelist')))
}
p <- ggplot (res, aes (y=Timing, x=N_iterations, colour=Build_method)) + geom_line()
print (p)
