\name{catarrhines}
\alias{catarrhines}
\docType{data}
\title{Phylogenetic tree of Catarrhini}
\description{
  Phylogenetic tree of the Parvorder Catarrhini, sister group of the New World Monkeys consisting of the Old World Monkeys and Apes, trimmed from Bininda-Edmonds et al. (2007)'s Mammalian supertree. Here used for testing and demonstrating MoreTreeTools functions.
}
\usage{
   data(catarrhines)
}
\format{
  \code{catarrhines} is a \code{phylo} from the pacakge \code{ape}.
}
\source{
  Bininda-Emonds, O. R. P., Cardillo, M., Jones, K. E., MacPhee, R. D. E., Beck, R. M. D., Grenyer, R., … Purvis, A. (2007). The delayed rise of present-day mammals. Nature, 446(7135), 507–12. doi:10.1038/nature05634
}
\examples{
data (catarrhines) # load object
plot (catarrhines) # plot
}
\keyword{datasets}
