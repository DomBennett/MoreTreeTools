\name{calcED}
\alias{calcED}
\title{Calculate Evolutionary Distinctiveness (Fair Proportion or Pendant Edge)}
\usage{
calcED(tree, tips = "all", type = "FP", progress = "none")
}
\arguments{
  \item{tree}{phylogenetic tree (ape class)}

  \item{tips}{vector of tips for which ED is calculated,
  else 'all'}

  \item{type}{method of ED calculation, either 'FP' or
  'PE'}

  \item{progress}{bar showing the progress of the
  calculation. Either time, text, tk or none. Uses
  \code{plyr}'s \code{create_progress_bar}.}
}
\description{
Calculate species evolutionary distinctiveness (ED) using
either of two methods: Fair Proportion (FP) or Pendant Edge
(PE).
}
\details{
Evolutionary distinctiveness is a measure of how much
independent evolution a species represents. Multiple
methods exist for its calculation all of which require an
ultrametric phylogenetic tree. The methods used here are
Pendant Edge (PE) the length of a species branch that
connects it to the tree (Altschul and Lipman, 1990) or Fair
Proportion (FP) the total proportion of the phylogenetic
tree that a species represents where each branch is equally
divided between all descdendants (Isaac et al. 2007).

N.B. \code{picante} already has a function for doing this.
\code{calcED}, however, uses \code{plyr} vectorisation and
should be faster.
}
\examples{
# example.var <- exampleFun (test.data)
}
\references{
No references yet!
\url{https://github.com/DomBennett/MoreTreeTools}.

Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and
Baillie, J.E.M. (2007).  Mammals on the EDGE: conservation
priorities based on threat and phylogeny. PLoS ONE, 2,
e296.

Altschul, S., and Lipman, D. (1990). Equal animals. Nature.
}
\keyword{ecology}
\keyword{evolution}
\keyword{manipulation}
\keyword{palaeontology}
\keyword{phylogeny}

