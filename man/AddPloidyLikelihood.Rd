\name{AddPloidyLikelihood}
\alias{AddPloidyLikelihood}
\alias{AddPloidyLikelihood.RADdata}
\title{
Likelihoods for Possible Ploidies Based on Genotype Distributions
}
\description{
Given prior genotype probabilities, and a set of high-confidence 
genotypes estimated with \code{\link{GetLikelyGen}}, this function
estimates the probability of observing that distribution of genotypes
and stores the probability in the \code{$ploidyLikelihood} slot of the
\code{"RADdata"} object.
}
\usage{
AddPloidyLikelihood(object, ...)
\method{AddPloidyLikelihood}{RADdata}(object, excludeTaxa = GetBlankTaxa(object), 
                            minLikelihoodRatio = 50, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  Prior genotype probabilities and
genotype likelihoods should have already been added using the appropriate
functions.
}
  \item{\dots}{
Additional arguments to be passed to the method for \code{"RADdata"}.
}
\item{excludeTaxa}{
A character vector indicating taxa that should be excluded from calculations.
}
\item{minLikelihoodRatio}{
A number, one or higher, to be passed to \code{GetLikelyGen}.
}
}
\details{
The purpose of this function is to estimate the correct inheritance mode for
each locus.  This function may be deleted in the future in favor of better
alternatives.
}
\value{
A \code{"RADdata"} object identical to that passed to the function, but with
results added to the \code{$ploidyLikelihood} slot.  This has one row for each
possible ploidy (each ploidy with data in \code{$priorProb}), and one column
for each allele.  Each element of the matrix is the multinomial probability 
of seeing that distribution of genotypes given the prior probabilities.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddPloidyChiSq}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    UseMethod("AddPloidyLikelihood", object)
  }
}
\keyword{ distribution }

