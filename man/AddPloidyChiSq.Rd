\name{AddPloidyChiSq}
\alias{AddPloidyChiSq}
\alias{AddPloidyChiSq.RADdata}
\title{
Chi-Square Test on Genotype Likelihood Distributions
}
\description{
This function is intended to help identify the correct inheritance mode for
each locus in a \code{"\link{RADdata}"} object.  Expected genotype frequencies
are taken from \code{object$priorProb}.  Observed genotype frequencies are
estimated from \code{object$genotypeLikelihood}, where each taxon has a
partial assignment to each genotype, proportional to genotype likelihoods.
A \eqn{\chi ^ 2}{chi-squared} statistic is then estimated.
}
\usage{
AddPloidyChiSq(object, ...)
\method{AddPloidyChiSq}{RADdata}(object, excludeTaxa = GetBlankTaxa(object),
               \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  Genotype prior probabilities and
likelihoods should have been added.
}
\item{excludeTaxa}{
A character vector indicating names of taxa to exclude from calculations.
}
  \item{\dots}{
Additional arguments to be passed to other methods (none currently in use).
}
}
\details{
Parents (in mapping populations) and blank taxa are automatically excluded
from calculations.

Genotypes with zero prior probability would result in an infinite 
A \eqn{\chi ^ 2}{chi-squared} statistic and therefore are excluded from
the calculation.  However, the total number of observations (total number
of taxa) remains the same, so that if there are many taxa with high
likelihood for a genotype with zero prior probability, 
\eqn{\chi ^ 2}{chi-squared} will be high.
}
\value{
A \code{"RADdata"} object identical to the one passed to the function,
but with a matrix added to the \code{$ploidyChiSq}
slot.  This matrix has inheritance rows (matching \code{object$priorProb}) in
rows and alleles in columns.  \code{object$ploidyChiSq} contains the 
\eqn{\chi ^ 2}{chi-squared} values.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypeLikelihood}}, \code{\link{AddPloidyLikelihood}}
}
\examples{
# load dataset and set some parameters
data(exampleRAD_mapping)
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")
exampleRAD_mapping <- AddAlleleFreqMapping(exampleRAD_mapping,
                                           expectedFreqs = c(0.25, 0.75),
                                           allowedDeviation = 0.08)
exampleRAD_mapping <- AddGenotypeLikelihood(exampleRAD_mapping)
exampleRAD_mapping <- AddGenotypePriorProb_Mapping2Parents(exampleRAD_mapping,
                                                      n.gen.backcrossing = 1)
                                                      
# get chi-squared values
exampleRAD_mapping <- AddPloidyChiSq(exampleRAD_mapping)
# view chi-squared and p-values (diploid only)
exampleRAD_mapping$ploidyChiSq
}

\keyword{ methods }
\keyword{ distribution }
