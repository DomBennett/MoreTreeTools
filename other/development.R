# Thinking about creating a class for faster tree manipulation

# TODO:
# -- create addNode
# -- create removeNode
# -- create updateNodeList (updates values given a change)
# -- create converter from Phylo to NodeList

# Process for creating a class
# 1. Use setClass to define the class
# 2. Use setGeneric to define a method for the class
# 3. Use setMethod to define th function for the method

# Node declaration
node_representation <- representation (children='vector',
                                       min_age='numeric',
                                       max_age='numeric',
                                       prev_node='integer',
                                       next_node='integer',
                                       node_id='integer')
setClass ('Node', representation=node_representation)
setGeneric ('print')
setMethod ('print', c('x'='Node'),
           function(x){
             cat ('Node: ID=[', x@node_id, '], children=[',
                  length (x@children), '], age=[',
                  x@max_age, ']', sep='')
           })
setMethod ('show', 'Node',
           function(object){
             print (object)
           })
# NodeList declaration
nodelist_representation <- representation (nodes='list',
                                           elements='vector',
                                           ntips='numeric')
prototype_root_node <- new ('Node', node_id=1L)
nodelist_prototype <- prototype (nodes=list(prototype_root_node),
                                 elements=c('children', 'age'),
                                 ntips=1)
setClass ('NodeList', representation=nodelist_representation,
          prototype=nodelist_prototype)
setGeneric ('print')
setMethod ('print', c('x'='NodeList'),
           function(x){
             cat ('NodeList: [',
                  length (x@nodes), '] nodes', sep='')
           })
setMethod ('show', 'NodeList',
           function(object){
             print (object)
           })
setMethod ('[', c ('NodeList', 'numeric', 'missing', 'ANY'),
           function(x, i, j, ..., drop=TRUE) {
             i <- as.integer (i)
             initialize(x, nodes=x@nodes[i])
           })
setMethod ('[<-', c ('NodeList', 'numeric', 'missing', 'NodeList'),
           function(x, i, j, ..., value) {
             i <- as.integer (i)
             x@nodes[i] <- initialize(value)
             initialize(x)
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
             initialize(x)
           })
update <- function (nodelist, node_id) {
  # Recursively update each nodelist element from node_id to root
  .update <- function (nodelist, node_id, next_node_id) {
    node <- nodelist[[node_id]]
    next_node <- nodelist[[next_node_id]]
    # for (element in nodelist@elements)
    node@children <- c (node@children, next_node@children)
    # TODO: fix age update
    test_age <- (node@age + next_node@age)
    node@age <- ifelse (test_age > node@age, next_node@age,
                        test_age)
    nodelist@nodes[[node_id]] <- node
    if (length (nodelist[[node_id]]@prev_node) > 0) {
      nodelist <- .update (nodelist, node@prev_node, node_id)
    }
    nodelist
  }
  node <- nodelist[[node_id]]
  if (length (node@prev_node) > 0) {
    nodelist <- .update (nodelist, node@prev_node, node_id)
  }
  nodelist
}

nodelist <- new('NodeList')
nodelist <- addNode(nodelist, parent_node=1L, min_age=0, max_age=1, tip='t1')
nodelist <- addNode(nodelist, parent_node=1L, min_age=0, max_age=0.5, tip='t2')

addNode <- function(nodelist, parent_node, min_age, max_age, tip=NA) {
  node_id <- as.integer (length (nodelist@nodes) + 1)
  children <- c (nodelist@nodes[[parent_node]]@children, tip)
  node <- new ('Node', min_age=min_age, max_age=max_age,
              prev_node=parent_node, node_id=node_id,
              children=children)
  nodelist@nodes[[node_id]] <- node
  nodelist <- update (nodelist, node_id)
  return (nodelist)
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
#                           prev_node=parent_node, node_id=node_id,
#                           children=children)
#              x@nodes[[node_id]] <- node
#              return (x)
#            })