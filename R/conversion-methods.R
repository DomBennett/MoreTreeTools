# register phylo as a pseudo S4
setOldClass ('phylo')

# phylo --> NodeList
# TODO: handle no branchlengths
setAs (from="phylo", to="NodeList",
       def=function (from, to) {
         # key stats
         tree_age <- getSize (from, 'rtt')
         # unique ids for each node (internal and tip)
         ids <- paste0 ('n', (length (from$tip.label)+1):
                          (length (from$tip.label) + from$Nnode))
         ids <- c (from$tip.label, ids)
         # init nodelist and add to list
         nodelist <- list ()
         for (i in 1:length (ids)) {
           pre_is <- from$edge[from$edge[ ,2] == i, 1]
           post_is <- from$edge[from$edge[ ,1] == i, 2]
           age <- getAge (from, node=i)[ ,2]
           predist <- tree_age - age
           pd <- sum (from$edge.length[getEdges (from, node=i)])
           children <- getChildren (from, node=i)
           span <- from$edge.length[from$edge[ ,2] == i]
           node <- list ('span'=span, 'id'=ids[i], 'pd'=pd,
                         'children'=children,
                         'prenode'=ids[pre_is],
                         'postnode'=ids[post_is],
                         'predist'=predist)
           nodelist[[ids[i]]] <- node
         }
         if (is.rooted (from)) {
           root_node <- ids[length (from$tip.label) + 1]
           nodelist[[root_node]]$span <- 0
           to <- new (to, nodelist=nodelist, root=root_node)
         } else {
           to <- new (to, nodelist=nodelist)
         }
         .update (to)
       })

# NodeList --> phylo
# TODO
# setAs (from='NodeLise', to='phylo',
#        def=function(from, to) {
#          
#        })