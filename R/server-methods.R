#' @name taxaResolve
#' @title Resolve taxonmic names online
#' @description Resolve taxonomic names via the Global Names Resolver.
#' @details Returns dataframe containing GNR metadata for each name wames
#' that cannot be resolved are returned as NA. Various datasources are 
#' available, see \url{http://resolver.globalnames.biodinfo.org/data_sources} for a
#' list and IDs. Default is 4 for NCBI.
#' @param nms vector of names
#' @param batch size of the batches to be queried
#' @param datasource ID number of the datasource
#' @param genus boolean, if true will search against GNR with just the genus
#'  name for names that failed to resolve using the full species name
#' @param cache T/F, create a local cache of downloaded names?
#' @param parent specify parent of all names to prevent false names
#' @export
#' @examples
#' my.lovely.names <- c ('Gallus gallus', 'Pongo pingu', 'Homo sapiens',
#'                       'Arabidopsis thaliana', 'Macaca thibetana', 'Bacillus subtilis')
#' res <- taxaResolve (nms=my.lovely.names)
#' length (colnames (res))  # 10 different metadata for returned names including original search name
#' # let's look at the lineages
#' lineages <- strsplit (as.vector (res$lineage), '\\|')
#' print (lineages[[6]])  # the bacteria has far fewer taxonomic levels

taxaResolve <- function (nms, batch=100, datasource=4, genus=TRUE,
                         cache=FALSE, parent=NULL) {
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
  batchResolve <- function (batch.nms) {
    #create query from nms
    url <- "http://resolver.globalnames.org/name_resolvers.json?"
    data_source_ids <- paste0 ("&data_source_ids=", datasource)
    nms2 <- paste0 ("names=", paste0 (stringr::str_replace_all (
      batch.nms, " ", "+"), collapse = "|"))
    query <- paste (plyr::compact (list (url, nms2,
                                   data_source_ids)), collapse = "")
    #search via API
    data <- .safeFromJSON (query)$data
    return (data)
  }
  # avoid names -- names exist on database but mean nothing
  avoid <- c ('unidentified')
  # make sure names don't have '_'
  trms <- gsub ('_', ' ', nms)
  # remove all non alphanumerics
  trms <- gsub ('\\W+', ' ', trms)
  # remove trailing whitespace
  trms <- gsub ('^\\s+|\\s+$', '', trms)
  # any missing trms, replace with stubs
  trms[trms==''] <- 'invalid'
  deja_vues <- rep(FALSE, length(trms))
  data <- vector("list", length=length(trms))
  names(data) <- nms
  if(cache) {
    if(!file.exists("gnr_cache")) {
      dir.create("gnr_cache")
    }
    for(i in 1:length(nms)) {
      fp <- file.path("gnr_cache", paste0(nms[i], ".RData"))
      if(file.exists(fp)) {
        load(fp)
        data[[nms[i]]] <- nd
        deja_vues[i] <- TRUE
      }
    }
  }
  if(sum(!deja_vues) > 0) {
    # Split nms into batch sized chunks
    #  http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
    x <- seq_along (trms)
    btrms <- split (trms[!deja_vues], ceiling (x/batch))
    bnms <- split (nms[!deja_vues], ceiling (x/batch))
    for (i in 1:length(btrms)) {
      temp.data <- batchResolve(btrms[[i]])
      data[bnms[[i]]] <- temp.data
    }
  }
  #transform results into output
  search.name <- name.string <- canonical.form <-
    lineage <- lineage.ids <- rank <- taxid <-
    match.type <- prescore <- score <- rep (NA, length (nms))
  for (i in 1:length (data)){
    parent_test <- TRUE
    nd <- data[[i]]
    if(cache & !deja_vues[i]) {
      fp <- file.path("gnr_cache", paste0(names(data)[i], ".RData"))
      save(nd, file=fp)
    }
    if (!'results' %in% names (nd)){
      search.name[i] <- nms[i]
    } else if (nd[[1]] %in% avoid) {
      search.name[i] <- nms[i]
    } else {
      search.name[i] <- nms[i]
      lng <- .replace(i, 'classification_path')
      if(!is.null(parent)) {
        parent_test <- grepl(parent, lng)
      }
      if(parent_test) {
        name.string[i] <- .replace(i, 'name_string')
        canonical.form[i] <- .replace(i, 'canonical_form')
        lineage[i] <- lng
        lineage.ids[i] <- .replace(i, 'classification_path_ids')
        rank[i] <- .replace(i, 'classification_path_ranks')
        taxid[i] <- .replace(i, 'taxon_id')
        match.type[i] <- .replace(i, 'match_type')
        prescore[i] <- .replace(i, 'prescore')
        score[i] <- nd$results[[1]]$score
      }
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
    genus.nms <- sub ('\\s+.*', '', res$search.name[failed])
    genus.res <- taxaResolve(genus.nms, batch, datasource, genus=FALSE,
                             parent=parent, cache=cache)
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
    json.obj <- try (RJSONIO::fromJSON(url), silent = TRUE)
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