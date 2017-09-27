\name{PipelineMapping2Parents}
\alias{PipelineMapping2Parents}
\title{
Run polyRAD Pipeline on a Mapping Population
}
\description{
This function is a wrapper for \code{\link{AddAlleleFreqMapping}},
\code{\link{AddGenotypeLikelihood}}, 
\code{\link{AddGenotypePriorProb_Mapping2Parents}},
\code{\link{AddPloidyChiSq}}, and \code{\link{AddGenotypePosteriorProb}}.
It covers the full pipeline for estimating genotype posterior probabilities
from read depth in a \code{"\link{RADdata}"} object containing data from
a mapping population.
}
\usage{
PipelineMapping2Parents(object, donorParent = GetDonorParent(object), 
                        recurrentParent = GetRecurrentParent(object), 
                        n.gen.backcrossing = 0, n.gen.intermating = 0, 
                        n.gen.selfing = 0, 
                        donorParentPloidies = object$possiblePloidies, 
                        recurrentParentPloidies = object$possiblePloidies, 
                        minLikelihoodRatio = 10, freqAllowedDeviation = 0.05, 
                        freqExcludeTaxa = c(GetDonorParent(object), 
                                            GetRecurrentParent(object), 
                                            GetBlankTaxa(object)))
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
  \item{donorParent}{
A character string indicating which taxon is the donor parent.  If 
backcrossing was not performed, it does not matter which taxon is listed
as the donor parent or recurrent parent.
}
  \item{recurrentParent}{
A character string indicating which parent is the recurrent parent.
}
  \item{n.gen.backcrossing}{
An integer, zero or greater, indicating how many generations of 
backcrossing to the recurrent parent were performed.
}
  \item{n.gen.intermating}{
An integer, zero or greater, indicating how many generations of 
intermating within the population were performed.
}
  \item{n.gen.selfing}{
An integer, zero or greater, indicatin ghow many generations of
selfing were performed.
}
  \item{donorParentPloidies}{
A list, where each item in the list is an integer vector indicating a
potential inheritance mode that could be observed among loci in the 
donor parent.  \code{2} indicates diploid, \code{4} indicates 
autotetraploid, \code{c(2, 2)} indicates, allotetraploid, \emph{etc}.
}
  \item{recurrentParentPloidies}{
A list in the same format as \code{donorParentPloidies} indicating
inheritance modes that could be observed among loci in the recurrent
parent.
}
  \item{minLikelihoodRatio}{
The minimum likelihood ration for determining parental genotypes with
confidence, to be passed to \code{\link{GetLikelyGen}} for both
parental taxa.
}
  \item{freqAllowedDeviation}{
For \code{\link{AddAlleleFreqMapping}}, the amount by which an allele
frequency can deviate from an expected allele frequency in order to be
counted as that allele frequency.
}
  \item{freqExcludeTaxa}{
A character vector indicating taxa to exclude from allele frequency
estimates and ploidy \eqn{\chi ^ 2}{chi-squared} estimates.
}
}
\details{
Unlike \code{\link{IterateHWE}} and \code{\link{IteratePopStruct}},
\code{PipelineMapping2Parents} only runs through each function once,
rather than iteratively until convergence.
}
\value{
A \code{"RADdata"} object identical to that passed to the function, with
the following slots added: \code{$alleleFreq}, \code{$genotypeLikelihood},
\code{$priorProb}, \code{$priorProbPloidies}, \code{$ploidyChiSq},
\code{$ploidyChiSqP}, and \code{$posteriorProb}.  See the documentation
for the functions listed in the description for more details on the
data contained in these slots.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{GetWeightedMeanGenotypes}} for exporting genotypes from 
the resulting object.
}
\examples{
# load data for the example
data(exampleRAD_mapping)

# specify donor and recurrent parents
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")

# run the pipeline
exampleRAD_mapping <- PipelineMapping2Parents(exampleRAD_mapping,
                                 n.gen.backcrossing = 1)
                                 
# export results
wmgeno <- GetWeightedMeanGenotypes(exampleRAD_mapping)[-(1:2),]
wmgeno
}
\keyword{ misc }
