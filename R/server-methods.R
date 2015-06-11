#' @name taxaResolve
#' @title Resolve taxonmic names online
#' @description Resolve taxonomic names via the Global Names Resolver.
#' @details Returns dataframe containing GNR metadata for each name wames
#' that cannot be resolved are returned as NA. Various datasources are 
#' available, see \url{http://resolver.globalnames.biodinfo.org/data_sources} for a
#' list and IDs. Default is 4 for NCBI.
#' @param names vector of names
#' @param batch size of the batches to be queried
#' @param datasource ID number of the datasource
#' @param genus boolean, if true will search against GNR with just the genus
#'  name for names that failed to resolve using the full species name
#' @export
#' @examples
#' my.lovely.names <- c ('Gallus gallus', 'Pongo pingu', 'Homo sapiens',
#'                       'Arabidopsis thaliana', 'Macaca thibetana', 'Bacillus subtilis')
#' res <- taxaResolve (names=my.lovely.names)
#' length (colnames (res))  # 10 different metadata for returned names including original search name
#' # let's look at the lineages
#' lineages <- strsplit (as.vector (res$lineage), '\\|')
#' print (lineages[[6]])  # the bacteria has far fewer taxonomic levels

taxaResolve <- function (names, batch=100, datasource=4, genus=TRUE) {
  .replace <- function (i, slot.name) {
    # controlled extraction
    element <- data[[i]]$result[[1]][[slot.name]]
    if (!is.null (element)) {
      res <- element
    } else {
      res <- NA
    }
    res
  }
  batchResolve <- function (batch.names) {
    #create query from names
    url <- "http://resolver.globalnames.org/name_resolvers.json?"
    data_source_ids <- paste0 ("&data_source_ids=", datasource)
    names2 <- paste0 ("names=", paste0 (str_replace_all (
      batch.names, " ", "+"), collapse = "|"))
    query <- paste (compact (list (url, names2,
                                   data_source_ids)), collapse = "")
    #search via API
    data <- .safeFromJSON (query)$data
    return (data)
  }
  # avoid names -- names exist on database but mean nothing
  avoid <- c ('unidentified')
  # remove unwanted reserved characters from names
  names <- gsub ('&|\\?|\\+|\\=', ' ', names)
  # make sure names don't have '_'
  names <- gsub ('_', ' ', names)
  data <- list()
  # Split names into batch sized chunks
  #  http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  x <- seq_along (names)
  di <- split (names, ceiling (x/batch))
  for (d in di) {
    temp.data <- batchResolve (d)
    data <- c (data, temp.data)
  }
  #transform results into output
  search.name <- name.string <- canonical.form <-
    lineage <- lineage.ids <- rank <- taxid <-
    match.type <- prescore <- score <- rep (NA, length (names))
  for (i in 1:length (names)){
    if (!'results' %in% names (data[[i]])){
      search.name[i] <- data[[i]][[1]]
    } else if (data[[i]][[1]] %in% avoid) {
      search.name[i] <- data[[i]][[1]]
    } else {
      search.name[i] <- data[[i]][[1]]
      name.string[i] <- .replace(i, 'name_string')
      canonical.form[i] <- .replace(i, 'canonical_form')
      lineage[i] <- .replace(i, 'classification_path')
      lineage.ids[i] <- .replace(i, 'classification_path_ids')
      rank[i] <- .replace(i, 'classification_path_ranks')
      taxid[i] <- .replace(i, 'taxon_id')
      match.type[i] <- .replace(i, 'match_type')
      prescore[i] <- .replace(i, 'prescore')
      score[i] <- data[[i]]$results[[1]]$score
    }
  }
  res <- data.frame (search.name=search.name,
                     name.string=name.string,
                     canonical.form=canonical.form,
                     lineage=lineage, lineage.ids=lineage.ids,
                     rank=rank, taxid=taxid,
                     match.type=match.type, prescore=prescore,
                     score=score, stringsAsFactors=FALSE)
  failed <- which (is.na (res$name.string))
  if (genus & length (failed) > 0) {
    # if genus, search just genus names
    genus.names <- sub ('\\s+.*', '', res$search.name[failed])
    genus.res <- taxaResolve (genus.names, batch, datasource, genus=FALSE)
    # replace in original results, all slots except search.name
    res[failed,-1] <- genus.res[ ,-1]
  }
  return (res)
}

.safeFromJSON <- function (url, max.trys=12, power=2) {
  # Safe wrapper for fromJSON
  trys <- 0
  waittime <- 2
  while (trys < max.trys) {
    json.obj <- try (fromJSON (url), silent = TRUE)
    if (class (json.obj) == 'try-error') {
      cat ('---- Connection failed: trying again in [', waittime,
           's]----\n', sep='')
      trys <- trys + 1
      Sys.sleep (waittime)
      waittime <- waittime*power
    } else {
      return (json.obj)
    }
  }
  stop ("Failed to connect, server may be down.")
}