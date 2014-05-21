.findClade <- function (lineages) {
  # for a list of lineages, find the clade shared by all
  subj <- lineages[[1]]
  for (i in 2:length (lineages)) {
    query <- lineages[[i]]
    subj <- subj[subj %in% query]
  }
  subj[length (subj)]
}