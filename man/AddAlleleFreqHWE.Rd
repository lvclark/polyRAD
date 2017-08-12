\name{AddAlleleFreqHWE}
\alias{AddAlleleFreqHWE}
\title{
Estimate Allele Frequencies in a RADdata Object Assuming Hardy-Weinberg Equilibrium
}
\description{
A simple estimator of allele frequencies in the dataset is performed based on
read depth. The column means of \code{object$depthRatio} are calculated and
added to \code{object} as the \code{alleleFreq} slot.
}
\usage{
AddAlleleFreqHWE(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
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
\code{\link{AddAlleleFreqMapping}}
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
