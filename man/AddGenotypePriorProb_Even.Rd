\name{AddGenotypePriorProb_Even}
\alias{AddGenotypePriorProb_Even}
\title{
Add Uniform Priors to a RADdata Object
}
\description{
To estimate genotype posterior probabilities based on read depth alone, without
taking any population parameters into account, this function can be used to set
a uniform prior probability on all possible genotypes.  This function is not
part of any pipeline but can be used for very rough and quick genotype
estimates, when followed by \code{\link{AddGenotypeLikelihood}},
\code{\link{AddGenotypePosteriorProb}}, \code{\link{AddPloidyChiSq}}, and
\code{\link{GetWeightedMeanGenotypes}} or \code{\link{GetProbableGenotypes}}.
}
\usage{
AddGenotypePriorProb_Even(object, ...)
}
\arguments{
  \item{object}{
A \code{\link{RADdata}} object.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\value{
A \dQuote{RADdata} object identical that passed to the function, but with
data stored in one new slot:

\item{priorProb }{A two-dimensional list of matrices, with rows corresponding
to \code{object$possiblePloidies} and columns corresponding to unique values in
\code{object$taxaPloidy}. Each item in the list is a matrix.  
For each matrix, allele copy number (from zero to the total
ploidy) is in rows, and alleles are in columns.  Each value is 
\eqn{1 / (ploidy + 1).}}
}

\note{
Values in \code{object$ploidyChiSq} may not be particularly meaningful
under uniform priors.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypePriorProb_HWE}}
}
\examples{
data(exampleRAD)

exampleRAD <- AddGenotypePriorProb_Even(exampleRAD)
exampleRAD$priorProb

# finish protocol to get genotypes
exampleRAD <- AddGenotypeLikelihood(exampleRAD)
exampleRAD <- AddPloidyChiSq(exampleRAD)
exampleRAD <- AddGenotypePosteriorProb(exampleRAD)

genmat <- GetWeightedMeanGenotypes(exampleRAD)
genmat
}

\keyword{ methods }
\keyword{ distribution }
