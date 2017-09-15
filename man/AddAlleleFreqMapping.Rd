\name{AddAlleleFreqMapping}
\alias{AddAlleleFreqMapping}
\alias{AddAlleleFreqMapping.RADdata}
\title{
Estimate Allele Frequencies in a Mapping Population
}
\description{
Estimate allele frequencies using data from a mapping population, assuming
a fixed set of allele frequencies are possible.
}
\usage{
AddAlleleFreqMapping(object, ...)
\method{AddAlleleFreqMapping}{RADdata}(object, expectedFreqs = seq(0, 1, 0.25),
                     allowedDeviation = 0.05,
                     excludeTaxa = c(GetDonorParent(object),
                                     GetRecurrentParent(object),
                                     GetBlankTaxa(object)), \dots)
}

\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  The donor and recurrent parent should have
been assigned with \code{\link{SetDonorParent}} and 
\code{\link{SetRecurrentParent}}, respectively.  If this is not a backcross
population, it does not matter which is the donor or recurrent parent.
}
\item{expectedFreqs}{
A numeric vector listing all expected allele frequencies in the mapping
population.
}
\item{allowedDeviation}{
A value indicating how far an observed allele frequency can deviate from an
expected allele frequency and still be categorized as that allele frequency.
Must be no more than half the smallest interval seen in \code{expectedFreqs}.
}
\item{excludeTaxa}{
A character vector indicating taxa that should be excluded from the allele 
frequency estimate.
}
  \item{\dots}{
Arguments to be passed to the method for \code{"RADdata"}.
}
}
\details{
Allele frequencies are first estimated as the column means of 
\code{object$depthRatio} (unless posterior genotype probabilities and ploidy
chi-squared values have already been calculated, in which case 
\code{\link{GetWeightedMeanGenotypes}} is run and the column means of its output are
taken), excluding any taxa listed in \code{excludeTaxa}.
These are then categorized based on which, if any, expected allele frequency
they match with, based on the intervals described by \code{expectedFreqs} and
\code{allowedDeviation}.  If an allele frequency does not fall within any of
these intervals it is classified as \code{NA}; otherwise it is converted to the
matching value in \code{expectedFreqs}.
}
\value{
A \code{"RADdata"} object identical to the one passed to the function, but with
allele frequencies added to \code{object$alleleFreq}, and \code{"mapping"}
as the \code{"alleleFreqType"} attribute.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddAlleleFreqHWE}}
}
\examples{
# load example dataset
data(exampleRAD_mapping)
exampleRAD_mapping

# specify parents
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")

# estimate allele frequencies in diploid BC1 population
exampleRAD_mapping <- AddAlleleFreqMapping(exampleRAD_mapping,
                        expectedFreqs = c(0.25, 0.75),
                        allowedDeviation = 0.08)
                        
exampleRAD_mapping$alleleFreq
}

\keyword{ methods }
\keyword{ arith }
