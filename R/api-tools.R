.safeFromJSON <- function (url, max.trys = 10) {
  # Safe wrapper for fromJSON
  trys <- 0
  while (trys < max.trys) {
    json.obj <- try (fromJSON (url), silent = TRUE)
    if (class (json.obj) == 'try-error') {
      cat ('---- Connection failed: trying again ----\n')
      trys <- trys + 1
      Sys.sleep (10)
    } else {
      return (json.obj)
    }
  }
  stop ("Failed to connect, server may be down.")
}

.taxaResolve <- function (names, batch = 100, datasource = 4){
  # Resolve taxonomic names via the Global Names Resolver.
  #  Names that cannot be resolved are returned as NA.
  #
  # Args:
  #  names: vector of names
  #  batch: the max batch number of names to be searched
  #
  # Return:
  #  dataframe containing GNR metadata for each name
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
    #print(i)
    if (length (data[[i]]) == 1){
      search.name[i] <- data[[i]][[1]]
    } else {
      search.name[i] <- data[[i]][[1]]
      name.string[i] <-
        data[[i]]$results[[1]]$name_string
      canonical.form[i] <-
        data[[i]]$results[[1]]$canonical_form
      lineage[i] <-
        data[[i]]$results[[1]]$classification_path
      lineage.ids[i] <-
        data[[i]]$results[[1]]$classification_path_ids
      rank[i] <-
        data[[i]]$results[[1]]$classification_path_ranks
      taxid[i] <-
        data[[i]]$results[[1]]$taxon_id
      match.type[i] <-
        data[[i]]$results[[1]]$match_type
      prescore[i] <-
        data[[i]]$results[[1]]$prescore
      score[i] <- data[[i]]$results[[1]]$score
    }
  }
  res <- data.frame (search.name = search.name,
                     name.string = name.string,
                     canonical.form = canonical.form,
                     lineage = lineage, lineage.ids =
                       lineage.ids, rank = rank,
                     taxid = taxid, match.type =
                       match.type, prescore = prescore,
                     score = score)
  return (res)
}