# Node methods -- see test-node-methods.R
.Node_sister <- function (node=.self) {
  ids <- sapply (node$prenode[[1]]$postnode,
                 function (x) x$id)
  i <- which (ids != node$id)
  node$prenode[[1]]$postnode[[i]]
}
.Node_postDist <- function (node, res=0) {
  i <- 1
  for (n in node$postnode) {
    res <- c (res[-i], .Node_postDist (n, res[i]))
    i <- 1 + length (res)
    res[i] <- 0
  }
  max (res + node$span)
}
.Node_preDist <- function (node, res=0) {
  i <- 1
  for (n in node$prenode) {
    res[i] <- res[i] + node$span
    res[i] <- .Node_preDist (n, res[i])
    i <- 1 + length (res)
    res[i] <- 0
  }
  max (res)
}
.Node_pd <- function (node, res=0) {
  for (n in node$postnode) {
    res <- res + n$span
    res <- .Node_pd (n, res)
  }
  res
}
.Node_prenodes <- function (node, ids=vector()) {
  for (n in node$prenode) {
    ids <- append (ids, n$id)
    ids <- .Node_prenodes (n, ids)
  }
  ids
}
.Node_postnodes <- function (node, ids=vector()) {
  for (n in node$postnode) {
    ids <- append (ids, n$id)
    ids <- .Node_postnodes (n, ids)
  }
  ids
}
.Node_children <- function (node, ids=vector()) {
  for (n in node$postnode) {
    if (length (n$postnode) == 0) {
      ids <- append (ids, n$id)
    }
    ids <- .Node_children (n, ids)
  }
  ids
}

#' @name Node
#' @title Node Class
#' @docType class
#' @section No details yet
#' @exportClass Node
#' @examples
#' # No examples yet

Node <- setRefClass ('Node',
                     fields=list (
                       id='character',  # id of current node
                       span='numeric',  # span of preceding edge
                       prenode='vector',  # vector of immediate rootwards node(s)
                       postnode='vector'),  # vector of immediate tipwards node(s)
                     methods=list (
                       'sister'=function() { # return sister Node
                         .Node_sister (.self)
                         },
                       'children'=function() {  # return children IDs
                         .Node_children(.self)
                       },
                       'preDist'=function() {  # return max dist root-wards
                         .Node_preDist(.self)
                         },
                       'postDist'=function() {  # return max dist tip-wards
                         .Node_postDist(.self)
                       },
                       'pd'=function() {  # return pd
                         .Node_pd(.self)
                       },
                       'prenodes'=function() {  # return all nodes root-wards
                         .Node_prenodes(.self)
                       }
                       ,
                       'postnodes'=function() {  # return all nodes tip-wards
                         .Node_postnodes(.self)
                       }))

# Display methods
setGeneric ('as.character')
setMethod ('as.character', c('x'='Node'),
           function(x) {
             paste0 ('Node Object: $id = [', x@id,']')
           })
setGeneric ('print')
setMethod ('print', c('x'='Node'),
           function(x){
             .listnodes <- function (x, e) {
               nodes <- sapply (x[[e]], function (x) x$id)
               paste (nodes, collapse=', ')
             }
             .concat <- function (msg, n, e) {
               paste0 (msg, '  $', n, ' = ', e, '\n') 
             }
             msg <- paste0('Node Object with fields:\n')
             msg <- .concat (msg, 'id', x$id)
             if (length (x$span) > 0) {
               msg <- .concat (msg, 'span', x$span)
             }
             if (length (x$prenode) > 0) {
               nodes <- .listnodes (x, 'prenode')
               msg <- .concat (msg, 'prenode', nodes)
             }
             if (length (x$postnode) > 0) {
               nodes <- .listnodes (x, 'postnode')
               msg <- .concat (msg, 'postnode', nodes)
             }
             cat (msg)
           })
setMethod ('show', 'Node',
           function(object){
             print (object)
           })
setMethod ('str', c('object'='Node'),
           function (object, max.level=2L, ...) {
             # Prevent inf regression
             if (is.na (max.level)) {
               stop ('max.level must be numeric')
             }
             str@default (object, max.level=max.level, ...)
           })