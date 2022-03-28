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
PipelineMapping2Parents(object, n.gen.backcrossing = 0,
                        n.gen.intermating = 0, n.gen.selfing = 0, 
                        minLikelihoodRatio = 10, freqAllowedDeviation = 0.05, 
                        freqExcludeTaxa = c(GetDonorParent(object), 
                                            GetRecurrentParent(object), 
                                            GetBlankTaxa(object)),
                        useLinkage = TRUE, linkageDist = 1e7,
                        minLinkageCorr = 0.5, overdispersion = 9)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
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
An integer, zero or greater, indicating how many generations of
selfing were performed.
}
  \item{minLikelihoodRatio}{
The minimum likelihood ratio for determining parental genotypes with
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
  \item{useLinkage}{
Boolean.  Should genotypes at nearby loci (according to genomic alignment
data) be used for updating genotype priors?
}
  \item{linkageDist}{
A number, in basepairs, indicating the maximum distance for linked loci.
Ignored if \code{useLinkage = FALSE}.
}
  \item{minLinkageCorr}{
A number ranging from zero to one.  Indicates the minimum correlation 
coeffienct between weighted mean genotypes at two alleles in order for linkage
data to be used for updating genotype priors. Ignored if 
\code{useLinkage = FALSE}.
}
\item{overdispersion}{
Overdispersion parameter; see \code{\link{AddGenotypeLikelihood}}.
}
}
\details{
Unlike \code{\link{IterateHWE}} and \code{\link{IteratePopStruct}},
\code{PipelineMapping2Parents} only runs through each function once,
rather than iteratively until convergence.
}
\value{
A \code{"RADdata"} object identical to that passed to the function, with
the following slots added: \code{$alleleFreq}, \code{depthSamplingPermutations},
\code{$genotypeLikelihood},
\code{likelyGeno_donor}, \code{likelyGeno_recurrent}, \code{$priorProb},
\code{$ploidyChiSq}, \code{$posteriorProb}, and if \code{useLinkage = TRUE},
\code{$alleleLinkages} and \code{$priorProbLD}.  See the documentation
for the functions listed in the description for more details on the
data contained in these slots.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{SetDonorParent}} and \code{SetRecurrentParent} to indicate which
individuals are the parents before running the function.

\code{\link{AddGenotypePriorProb_Mapping2Parents}} for how ploidy of parents and
progeny is interpreted.

\code{\link{GetWeightedMeanGenotypes}} or \code{\link{Export_polymapR}} for
exporting genotypes from the resulting object.

\code{\link{StripDown}} to remove memory-hogging slots that are no longer
needed after the pipeline has been run.
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
