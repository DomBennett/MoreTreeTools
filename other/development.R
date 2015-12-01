# test new NodeList
tree_phylo <- rtree (5)
tree <- as (tree_phylo, 'NodeList')

getSize (tree_phylo, 'rtt')
sum (tree_phylo$edge.length)

plot (tree_phylo);nodelabels()
getAge (tree_phylo, 5)
tree[['t5']]$age
