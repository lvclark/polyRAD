\name{CanDoGetWeightedMeanGeno}
\alias{CanDoGetWeightedMeanGeno}
\title{
Check Whether \code{GetWeightedMeanGenotypes} Can Be Run
}
\description{
This function is used internally by \code{\link{AddPCA}}, 
\code{\link{AddAlleleFreqByTaxa}}, and the internal function \code{.alleleFreq}
to test whether \code{\link{GetWeightedMeanGenotypes}} can be run on a
\code{"\link{RADdata}"} object.
}
\usage{
CanDoGetWeightedMeanGeno(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}

\value{
A single Boolean value.  To be \code{TRUE}, \code{object$posteriorProb} must be
non-null, and either there must be only one possible ploidy, or 
\code{object$ploidyChiSq} must be non-null.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypePosteriorProb}}, \code{\link{AddPloidyChiSq}}
}
\examples{
data(exampleRAD)

CanDoGetWeightedMeanGeno(exampleRAD)

exampleRAD <- AddAlleleFreqHWE(exampleRAD)
exampleRAD <- AddGenotypePriorProb_HWE(exampleRAD)
exampleRAD <- AddGenotypeLikelihood(exampleRAD)
exampleRAD <- AddPloidyChiSq(exampleRAD)
exampleRAD <- AddGenotypePosteriorProb(exampleRAD)

CanDoGetWeightedMeanGeno(exampleRAD)
}

\keyword{ methods }
\keyword{ utilities }
