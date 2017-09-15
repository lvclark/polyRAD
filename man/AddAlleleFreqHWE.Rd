\name{AddAlleleFreqHWE}
\alias{AddAlleleFreqHWE}
\alias{AddAlleleFreqHWE.RADdata}
\title{
Estimate Allele Frequencies in a RADdata Object Assuming Hardy-Weinberg Equilibrium
}
\description{
Allele frequencies are estimated based on the best parameters available.
\code{object$alleleFreqByTaxa} is used if available.  If \code{object$alleleFreqByTaxa}
is null, \code{\link{GetWeightedMeanGenotypes}} is used, and if that isn't possible
\code{object$depthRatio} is used.  From whichever of the three options is used, 
column means are taken, the output of which is stored as \code{object$alleleFreq}.
}
\usage{
AddAlleleFreqHWE(object, ...)
\method{AddAlleleFreqHWE}{RADdata}(object, excludeTaxa = GetBlankTaxa(object), \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{excludeTaxa}{
A character vector indicating taxa that should be excluded from the calculation.
}
  \item{\dots}{
Included to allow more arguments in the future, although none are currently used.
}
}
\value{
A \code{"RADdata"} object identical to the one passed to the function, but 
with allele frequencies added to \code{object$alleleFreq}, and \code{"HWE"}
as the value for the \code{"alleleFreqType"} attribute.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddAlleleFreqMapping}}, \code{\link{AddGenotypePriorProb_HWE}}
}
\examples{
# load in an example dataset
data(exampleRAD)
exampleRAD

# add allele frequencies
exampleRAD <- AddAlleleFreqHWE(exampleRAD)
exampleRAD$alleleFreq
}

\keyword{ arith }
\keyword{ methods }
