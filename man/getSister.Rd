\name{getSister}
\alias{getSister}
\title{Return sister node}
\usage{
getSister(tree, node = "all")
}
\arguments{
  \item{tree}{phylogenetic tree (ape class)}

  \item{node}{number indicating node}
}
\description{
Return sister node of node(s) given, default will return
sister nodes of all internal nodes.
}
\details{
A sister node is defined as the other node descending from
the parent node. This functions finds all sister nodes. It
does not handle polytomies.
}
\examples{
# example.var <- exampleFun (test.data)
}
\references{
No references yet!
\url{https://github.com/DomBennett/MoreTreeTools}.
}
\keyword{ecology}
\keyword{evolution}
\keyword{manipulation}
\keyword{palaeontology}
\keyword{phylogeny}

