\name{AddGenotypePriorProb_Mapping2Parents}
\alias{AddGenotypePriorProb_Mapping2Parents}
\alias{AddGenotypePriorProb_Mapping2Parents.RADdata}
\alias{EstimateParentalGenotypes}
\alias{EstimateParentalGenotypes.RADdata}
\title{
Expected Genotype Frequencies in Mapping Populations
}
\description{
\code{EstimateParentalGenotypes} estimates the most likely genotypes of two
parent taxa.  Using those parental genotypes,
\code{AddGenotypePriorProb_Mapping2Parents} estimates expected genotype
frequencies for a population of progeny, which are added to the
\code{"\link{RADdata}"} object in the \code{$priorProb} slot.
}
\usage{
AddGenotypePriorProb_Mapping2Parents(object, ...)
\method{AddGenotypePriorProb_Mapping2Parents}{RADdata}(object,
    donorParent = GetDonorParent(object), 
    recurrentParent = GetRecurrentParent(object), 
    n.gen.backcrossing = 0, n.gen.intermating = 0, n.gen.selfing = 0,
    minLikelihoodRatio = 10, \dots)
    
EstimateParentalGenotypes(object, ...)
\method{EstimateParentalGenotypes}{RADdata}(object,
    donorParent = GetDonorParent(object), 
    recurrentParent = GetRecurrentParent(object), 
    n.gen.backcrossing = 0, n.gen.intermating = 0, n.gen.selfing = 0,
    minLikelihoodRatio = 10, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  Ideally this should be set up as a mapping
population using \code{\link{SetDonorParent}}, \code{\link{SetRecurrentParent}},
and \code{\link{AddAlleleFreqMapping}}.
}
  \item{\dots}{
Additional arguments, listed below, to be passed to the method for 
\code{"RADdata"} objects.
}
\item{donorParent}{A character string indicating which taxon is the donor 
parent.  If backcrossing was not performed, it does not matter which was the
donor or recurrent parent.}
\item{recurrentParent}{A character string indicating which taxon is the
recurrent parent.}
\item{n.gen.backcrossing}{
An integer, zero or greater, indicating how many generations of backcrossing
to the recurrent parent were performed.
}
\item{n.gen.intermating}{
An integer, zero or greater, indicating how many generations of intermating
within the population were performed.  (Values above one should not have
an effect on the genotype priors that are output, \emph{i.e.} 
genotype probabilities after one generation of random mating are identical
to genotype probabilities after >1 generation of random mating, assuming no
genetic drift or selection).
}
\item{n.gen.selfing}{
An integer, zero or greater, indicating how many generations of selfing were
performed.
}
\item{minLikelihoodRatio}{
The minimum likelihood ratio for determining parental genotypes with 
confidence, to be passed to \code{\link{GetLikelyGen}} for both
parental taxa.
}
}
\details{
\code{AddGenotypePriorProb_Mapping2Parents} examines the parental and progeny
ploidies stored in \code{object$taxaPloidy} and throws an error if they do not
meet expectations.  In particular, all progeny must be the same ploidy, and that
must be the ploidy that would be expected if the parents produced normal gametes.
For example in an F1 cross, if one parent was diploid and the other tetraploid,
all progeny must be triploid.  If both parents are tetraploid, all progeny must
be tetraploid.

The most likely genotypes for the two parents are estimated by
\code{EstimateParentalGenotypes} using
\code{\link{GetLikelyGen}}.  If parental gentoypes don't match progeny allele
frequencies, the function attempts to correct the parental genotypes to the
most likely combination that matches the allele frequency.

For each ploidy being examined, F1 genotype probabilities are then calculated
by \code{AddGenotypePriorProb_Mapping2Parents}.
Genotype probabilities are updated for each backcrossing generation, then each
intermating generation, then each selfing generation.

The default, with \code{n.gen.backcrossing = 0}, \code{n.gen.intermating = 0}
and \code{n.gen.selfing = 0}, will simulate an F1 population.  A BC1F2
population, for example, would have \code{n.gen.backcrossing = 1}, 
\code{n.gen.intermating = 0} and \code{n.gen.selfing = 1}.  A typical F2 
population would have \code{n.gen.selfing = 1} and the other two parameters
set to zero.  However, in a self-incompatible species where many F1 are 
intermated to produce the F2, one would instead use 
\code{n.gen.intermating = 1} and set the other parameters to zero.
}
\value{
A \code{"RADdata"} object identical to that passed to the function, but with
data stored in three new slots:

\item{priorProb }{A two-dimensional list of matrices, with rows corresponding to
\code{object$possiblePloidies} and columns corresponding to unique values in
\code{object$taxaPloidy}. Each item in the list is a matrix.  
For each matrix, allele copy number (from zero to the total
ploidy) is in rows, and alleles are in columns.  Each value is the probability
of sampling an individual with that allele copy number from the population.
For any taxa ploidies other than the progeny ploidy, even priors are returned.}
\item{likelyGeno_donor }{A matrix of the donor parent genotypes that were
used for estimating genotype prior probabilities.  Formatted like the
output of \code{\link{GetLikelyGen}}.}
\item{likelyGeno_recurrent }{A matrix of the recurrent parent genotypes
that were use for estimating gentoype prior probabilities.}
}
\author{
Lindsay V. Clark
}
\note{
For the time being, in allopolyploids it is assumed that copies of an allele
are distributed among as few isoloci as possible.  For example, if an
autotetraploid genotype had two copies of allele A and two copies of allele
B, it is assumed to be AA BB rather than AB AB.  This may be remedied in the
future by examining distribution of genotype likelihoods.
}

\seealso{
\code{\link{AddGenotypeLikelihood}}, \code{\link{AddGenotypePriorProb_HWE}}
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

# examine the dataset
exampleRAD_mapping
exampleRAD_mapping$alleleFreq

# estimate genotype priors for a BC1 population
exampleRAD_mapping <- AddGenotypePriorProb_Mapping2Parents(exampleRAD_mapping,
                                                      n.gen.backcrossing = 1)
exampleRAD_mapping$priorProb
}
\keyword{ methods }
\keyword{ distribution }
