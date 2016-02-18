

library (MoreTreeTools)

## Function


# explore
tree <- rtree (20)
Rprof(tmp <- tempfile())
chromatophylo (tree)
Rprof()
summaryRprof(tmp)
