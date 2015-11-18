# register phylo as a pseudo S4
setOldClass ('phylo')

# phylo --> NodeList
setAs (from="phylo", to="NodeList",
       def=function (from, to) {
         # unique ids for each node (internal and tip)
         ids <- paste0 ('n', (length (from$tip.label)+1):
                          (length (from$tip.label) + from$Nnode))
         ids <- c (from$tip.label, ids)
         # init nodelist and add to list
         nodelist <- list ()
         for (i in 1:length (ids)) {
           span <- from$edge.length[from$edge[ ,2] == i]
           node <- list ('span'=span, 'id'=ids[i])
           nodelist[[ids[i]]] <- node
         }
         # add pre and post nodes
         for (i in 1:length (ids)) {
           # get prenodes
           prenodes_ids <- ids[from$edge[from$edge[ ,2] == i, 1]]
           prenodes <- vector ()
           for (id in prenodes_ids) {
             prenodes <- append (prenodes, nodelist[[id]]$id)
           }
           nodelist[[ids[i]]]$prenode <- prenodes
           # get postnodes
           postnodes_ids <- ids[from$edge[from$edge[ ,1] == i, 2]]
           postnodes <- vector ()
           for (id in postnodes_ids) {
             postnodes <- append (postnodes, nodelist[[id]]$id)
           }
           nodelist[[ids[i]]]$postnode <- postnodes
         }
         if (is.rooted (from)) {
           root_node <- ids[length (from$tip.label) + 1]
           nodelist[[root_node]]$span <- 0
           to <- new (to, nodelist=nodelist, root=root_node)
         } else {
           to <- new (to, nodelist=nodelist)
         }
         to
       })

# NodeList --> phylo
# TODO
# setAs (from='NodeLise', to='phylo',
#        def=function(from, to) {
#          
#        })