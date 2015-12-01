# test new NodeList
tree_phylo <- rtree (5)
tree <- as (tree_phylo, 'NodeList')

getSize (tree_phylo, 'rtt')
sum (tree_phylo$edge.length)

plot (tree_phylo);nodelabels();axisPhylo()
getAge (tree_phylo, 5)
age (tree)
tree[['t2']]$predist


# adding new tip
edge <- 'n7'
tip_name <- 't6'
node_name <- 'n10'
node_age <- 1.8
tip_age <- -1  # test whether negative numbers increases age of tree
tree <- MoreTreeTools:::.addTip__NodeList (tree, edge, tip_name, node_name,
                           node_age, tip_age)
tree
