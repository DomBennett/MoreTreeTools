#' @name NodeList
#' @title NodeList Class
#' @docType class
#' @section No details yet
#' @exportClass NodeList
#' @examples
#' # No examples yet
nodelist_representation <- representation (nodes='list')  # list of Node objects
setClass ('NodeList', representation=nodelist_representation)
setGeneric ('print')
setMethod ('print', c('x'='NodeList'),
           function(x){
             msg <- 'NodeList Object representing:\n'
             msg <- paste0 (msg, '  - [', length (x@nodes), '] nodes\n')
             #              msg <- paste0 (msg, '  - [', length (x@tip_indexes), '] tips\n')
             #              if (length (x@root_index) == 1) {
             #                msg <- paste0 (msg, '  - [', x@root_index, '] root node index\n')
             #              } else {
             #                msg <- paste0 (msg, '  - unrooted')
             #              }
             cat (msg)
           })
setMethod ('show', 'NodeList',
           function(object){
             print (object)
           })
setMethod ('str', c('object'='NodeList'),
           function (object, max.level=2L, ...) {
             str@default (object, max.level=max.level, ...)
           })
setMethod ('[[', c ('NodeList', 'character', 'missing'),
           function(x, i, j, ...) {
             x@nodes[[i]]
           })
setGeneric ('add', signature=c('x', 'node'),
            function(x, node) {
              genericFunction ('add')
            })
setMethod ('add', c ('NodeList', 'Node'),
           function (x, node) {
             x@nodes[[node$id]] <- node
             initialize (x)
           })