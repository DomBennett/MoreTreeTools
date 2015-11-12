library (ape)
tree <- rtree (3)

# unique ids for each node (internal and tip)
ids <- paste0 ('n', (length (tree$tip.label)+1):
  (length (tree$tip.label) + tree$Nnode))
ids <- c (tree$tip.label, ids)
# init nodes and add to list
nodes <- list ()
for (i in 1:length (ids)) {
  span <- tree$edge.length[tree$edge[ ,2] == i]
  node <- new ('Node', span=span, id=ids[i])
  nodes[[ids[i]]] <- node
}
# add pre and post nodes
for (i in 1:length (ids)) {
  # get prenodes
  prenodes_ids <- ids[tree$edge[tree$edge[ ,2] == i, 1]]
  prenodes <- vector ()
  for (id in prenodes_ids) {
    prenodes <- append (prenodes, nodes[[id]])
  }
  nodes[[ids[i]]]$prenode <- prenodes
  # get postnodes
  postnodes_ids <- ids[tree$edge[tree$edge[ ,1] == i, 2]]
  postnodes <- vector ()
  for (id in postnodes_ids) {
    postnodes <- append (postnodes, nodes[[id]])
  }
  nodes[[ids[i]]]$postnode <- postnodes
}
nodelist <- new('NodeList', nodes=nodes)

# test
getAge(tree, node=5)
nodelist[['n5']]$age()
node <- which (tree$tip.label == 't1')
tree$tip.label[getSister (tree, node)]
nodelist[['t1']]$sister()
getChildren (tree, node=5)
nodelist[['n5']]$children()

# setGeneric ('addNode', signature= c('x', 'parent_node', 'min_age',
#                                     'max_age', 'tip'),
#             function(x, parent_node, min_age, max_age, tip) {
#               genericFunction ('addNode')
#             })
# setMethod ('addNode', c('x'='NodeList', 'parent_node'='integer',
#                         'min_age'='numeric', 'max_age'='numeric',
#                         'tip'=ANY),
#            function(x, parent_node, min_age, max_age, tip=NA) {
#              node_id <- as.integer (length (x@nodes) + 1)
#              children <- c (x@nodes[[parent_node]]@children, tip)
#              node <- new ('Node', min_age=min_age, max_age=max_age,
#                           down_node=parent_node, node_id=node_id,
#                           children=children)
#              x@nodes[[node_id]] <- node
#              return (x)
#            })

# setMethod ('[', c ('NodeList', 'numeric', 'missing', 'ANY'),
#            function(x, i, j, ..., drop=TRUE) {
#              i <- as.integer (i)
#              initialize(x, nodes=x@nodes[i])
#            })
# setMethod ('[<-', c ('NodeList', 'numeric', 'missing', 'NodeList'),
#            function(x, i, j, ..., value) {
#              # stop this for now
#              stop ('Single nodes only can be replaced, use [[<-')
#              #i <- as.integer (i)
#              #x@nodes[i] <- initialize(value)
#              #initialize(x)
#            })
# setMethod ('[[', c ('NodeList', 'numeric', 'missing'),
#            function(x, i, j, ...) {
#              i <- as.integer (i)
#              initialize(x@nodes[[i]])
#            })
# setMethod ('[[<-', c ('NodeList', 'numeric', 'missing', 'Node'),
#            function(x, i, j, ..., value) {
#              i <- as.integer (i)
#              x@nodes[[i]] <- initialize (value)
#              # update elements given new node
#              x <- update(x)
#              initialize(x)
#            })
# setGeneric ('get', signature= c('x', 'element'),
#             function(x, element) {
#               genericFunction ('get')
#             })
# setMethod ('get', c ('NodeList', 'character'),
#            function (x, element) {
#              if (element == 'ntips') {
#                return (length (x@tip_indexes))
#              } else if (element == 'tip_labels') {
#                return (unlist (lapply (nodelist@nodes, slot, 'label')))
#              } else {
#                stop ('Invalid element')
#              }
#            })