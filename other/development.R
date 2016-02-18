# Thinking about creating a class for faster tree manipulation

# TODO:
# -- create addNode
# -- create removeNode
# -- create updateNodeList (updates values given a change)
# -- create converter from Phylo to NodeList

#S4
# Process for creating a class
# 1. Use setClass to define the class
# 2. Use setGeneric to define a method for the class
# 3. Use setMethod to define th function for the method

# Node declaration
node_fields <- list (node_id='integer',  # id of current node
                     label='character',  # name of node
                     children='vector',  # children descending
                     age='numeric',  # age of node in tree
                     pd='numeric',  # phylogenetic diversity from node
                     span='numeric',  # span of preceding edge
                     down_nodes='vector',  # id of nodes from root to node
                     up_nodes='vector')  # id of nodes from node to tips
Node <- setRefClass ('Node', fields=node_fields)
setGeneric ('as.character')
setMethod ('as.character', c('x'='Node'),
           function(x) {
             paste0 ('Node Object (ID=[', x@node_id,'])')
           })
setGeneric ('print')
setMethod ('print', c('x'='Node'),
           function(x){
             field_names <- names (Node$fields())
             nslots <- length (field_names)
             msg <- paste0('Node Object with fields:\n')
             for (n in field_names) {
               e <- x[[n]]
               if (length (e) > 0) {
                 msg <- paste0 (msg, '  $', n, ' = ', e, '\n')
               }
             }
             cat (msg)
           })
setMethod ('show', 'Node',
           function(object){
             print (object)
           })
# NodeList declaration
nodelist_representation <- representation (nodes='list',  # list of Node classes
                                           tip_indexes='vector',  # index of tip nodes
                                           root_index='integer')  # index of root node
prototype_root_node <- new ('Node', node_id=1L)
nodelist_prototype <- prototype (nodes=list(prototype_root_node),
                                 elements=c('children', 'age'),
                                 ntips=1)
setClass ('NodeList', representation=nodelist_representation,
          prototype=nodelist_prototype)
setGeneric ('print')
setMethod ('print', c('x'='NodeList'),
           function(x){
             msg <- 'NodeList Object representing:\n'
             msg <- paste0 (msg, '  - [', length (x@nodes), '] nodes\n')
             msg <- paste0 (msg, '  - [', length (x@tip_indexes), '] tips\n')
             if (length (x@root_index) == 1) {
               msg <- paste0 (msg, '  - [', x@root_index, '] root node index\n')
             } else {
               msg <- paste0 (msg, '  - unrooted')
             }
             cat (msg)
           })
setMethod ('show', 'NodeList',
           function(object){
             print (object)
           })
setMethod ('[', c ('NodeList', 'numeric', 'missing', 'ANY'),
           function(x, i, j, ..., drop=TRUE) {
             i <- as.integer (i)
             # downdate elements and drop hanging nodes
             x <- downdate (x)
             initialize(x, nodes=x@nodes[i])
           })
setMethod ('[<-', c ('NodeList', 'numeric', 'missing', 'NodeList'),
           function(x, i, j, ..., value) {
             # stop this for now
             stop ('Single nodes only can be replaced, use [[<-')
             #i <- as.integer (i)
             #x@nodes[i] <- initialize(value)
             #initialize(x)
           })
setMethod ('[[', c ('NodeList', 'numeric', 'missing'),
           function(x, i, j, ...) {
             i <- as.integer (i)
             initialize(x@nodes[[i]])
           })
setMethod ('[[<-', c ('NodeList', 'numeric', 'missing', 'Node'),
           function(x, i, j, ..., value) {
             i <- as.integer (i)
             x@nodes[[i]] <- initialize (value)
             # update elements given new node
             x <- update(x)
             initialize(x)
           })
setGeneric ('get', signature= c('x', 'element'),
            function(x, element) {
              genericFunction ('get')
            })
setMethod ('get', c ('NodeList', 'character'),
           function (x, element) {
             if (element == 'ntips') {
               return (length (x@tip_indexes))
             } else if (element == 'tip_labels') {
               return (unlist (lapply (nodelist@nodes, slot, 'label')))
             } else {
               stop ('Invalid element')
             }
           })
children <- function (node, up_node) {
  temp_children <- c (node@children, up_node@children)
  if (temp_children != node@children) {
    node@children <- temp_children
  }
  node
}
age <- function (node, up_node) {
  temp_age <- (node@age + up_node@age)
  if (temp_age > node@age) {
    node@age <- temp_age
  }
  node
}
update <- function (nodelist) {
  # Recursively update from tip nodes to root
  .update <- function (nodelist, node_id, up_node_id) {
    node <- nodelist[[node_id]]
    up_node <- nodelist[[up_node_id]]
    node <- children (node, up_node)
    node <- age (node, up_node)
    if (temp_children != node@children) {
      nodelist@nodes[[node_id]] <- node
      if (length (nodelist[[node_id]]@down_node) > 0) {
        nodelist <- .update (nodelist, node@down_node, node_id)
      }
    }
    nodelist
  }
  node <- nodelist[[node_id]]
  if (length (node@down_node) > 0) {
    nodelist <- .update (nodelist, node@down_node, node_id)
  }
  nodelist
}

# process for adding
# 1. create Node
# 2. add Node to nodelist
# 3. update connected nodes

addNode <- function(nodelist, node_id, label, parent_node_id,
                    span, age, tip) {
  node <- new ('Node', span=span, age=age, down_node=parent_node_id,
               node_id=node_id, children=tip, label=label)
  # add new node to list
  nodelist@nodes[[node_id]] <- node
  # update
  #update (nodelist, node_id=i)
  return (nodelist)
}

library (ape)
tree <- rtree (3)

node <- new ('Node', span=0, age=0, node_id=1L, label='n1')
node <- new ('Node', span=0, age=0, node_id=2L, label='t1', down_node=node)
nodelist <- new('NodeList')
nodelist[[1]] <- node
nodelist[[2]] <- node

# Get each Node from tree and put into NodeList
nodelist <- new('NodeList')
phylo_nodes <- 1:(length (tree$tip.label) + tree$Nnode)
# init nodelist
for (i in phylo_nodes) {
  nodelist <- addNode(nodelist, node_id=i, label=as.character (i),
                      parent_node_id=integer(),
                      age=0, span=0, tip=vector())
}
# fill in details
for (i in phylo_nodes) {
  edge <- which (tree$edge[ ,2] == i)
  span <- tree$edge.length[edge]
  if (length (age) == 0) {
    age <- 0
  }
  parent_node_id <- tree$edge[tree$edge[ ,2] == i, 1]
  tip <- tree$tip.label[i]
  if (is.na (tip)) {
    tip <- vector ()
  }
  nodelist <- addNode(nodelist, node_id=i,
                      parent_node_id=parent_node_id,
                      span=span, tip=tip)
}


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