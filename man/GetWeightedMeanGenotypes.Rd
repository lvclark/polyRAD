\name{GetWeightedMeanGenotypes}
\alias{GetWeightedMeanGenotypes}
\alias{GetWeightedMeanGenotypes.RADdata}
\alias{GetProbableGenotypes}
\alias{GetProbableGenotypes.RADdata}
\title{
Export Numeric Genotype Values from Posterior Probabilities
}
\description{
These functions calculate numerical genotype values using posterior
probabilities in a \code{"\link{RADdata}"} object, and output
those values as a matrix of taxa by alleles.
\code{GetWeightedMeanGenotypes} returns continuous genotype values, 
weighted by posterior genotype probabilities.  
\code{GetProbableGenotypes} returns discrete genotype values indicating
the most probable genotype.  If the \code{"\link{RADdata}"} 
object includes more than one possible inheritance mode, the 
\code{$ploidyChiSq} slot is used for selecting or weighting 
inheritance modes for each allele.
}
\usage{
GetWeightedMeanGenotypes(object, ...)
\method{GetWeightedMeanGenotypes}{RADdata}(object, minval = 0, maxval = 1, 
                         omit1allelePerLocus = TRUE, 
                         omitCommonAllele = TRUE,
                         naIfZeroReads = FALSE, 
                         onePloidyPerAllele = FALSE, \dots)

GetProbableGenotypes(object, ...)
\method{GetProbableGenotypes}{RADdata}(object, omit1allelePerLocus = TRUE,
                     omitCommonAllele = TRUE,
                     naIfZeroReads = FALSE, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  Posterior genotype probabilities should
have been added with \code{\link{AddGenotypePosteriorProb}}, and if there is
more than one possible ploidy,
ploidy chi-squared values should have been added with 
\code{\link{AddPloidyChiSq}}.
}
  \item{\dots}{
Additional arguments, listed below, to be passed to the method for
\code{"RADdata"}.
}
\item{minval}{
The number that should be used for indicating that a taxon has zero
copies of an allele.
}
\item{maxval}{
The number that should be used for indicating that a taxon has the
maximum copies of an allele (equal to the ploidy of the locus).
}
\item{omit1allelePerLocus}{
A logical indicating whether one allele per locus should be omitted
from the output, in order to reduce the number of variables and
prevent singularities for 
genome-wide association and genomic prediction.  The value for one
allele can be predicted from the values from all other alleles at
its locus.
}
\item{omitCommonAllele}{
A logical, passed to the \code{commonAllele} argument of 
\code{\link{OneAllelePerMarker}}, indicating whether the most common allele
for each locus should be omitted (as opposed to simply the first allele
for each locus).  Ignored if \code{omit1allelePerLocus = FALSE}.
}
\item{naIfZeroReads}{
A logical indicating whether \code{NA} should be inserted into the 
output matrix for any taxa and loci where the total read depth for
the locus is zero.  If \code{FALSE}, the output for these genotypes is
essentially calculated using prior genotype probabilities, since
prior and posterior genotype probabilities are equal when there are no
reads.
}
\item{onePloidyPerAllele}{
Logical.  If \code{TRUE}, for each allele the inheritance mode with the
lowest \eqn{\chi ^ 2}{chi-squared} value is selected and is assumed to be
the true inheritance mode.  If \code{FALSE}, inheritance modes are weighted
by inverse \eqn{\chi ^ 2}{chi-squared} values for each allele, and mean
genotypes that have been weighted across inheritance modes are returned.
}
}

\details{
For each inheritance mode \eqn{m}, taxon \eqn{t}, allele \eqn{a}, allele copy number 
\eqn{i}, total ploidy \eqn{k}, and posterior genotype probability \eqn{p_{i,t,a,m}},
weighted mean genotype \eqn{g_{t,a,m}} is estimated by \code{GetWeightedMeanGenotypes} as:

\deqn{g_{t,a,m} = \sum_{i = 0}^k p_{i,t,a,m} * \frac{i}{k}}{g_{t,a,m} = \sum_{i = 0}^k p_{i,t,a,m} * i/k}

For \code{GetProbableGenotypes}, the genotype is the one with the maximum posterior
probability:

\deqn{g_{t,a,m} = i | \max_{i = 0}^k{p_{i,t,a,m}}}

When there are multiple inheritance modes and \code{onePloidyPerAllele = FALSE},
the weighted genotype is estimated by \code{GetWeightedMeanGenotypes} as:

\deqn{g_{t,a} = \sum_m [ g_{t,a,m} * \frac{1}{\chi^2_{m,a}} / \sum_m \frac{1}{\chi^2_{m,a}}]}{
\sum_m [ g_{t,a,m} * 1/\chi^2_{m,a} / \sum_m 1/\chi^2_{m,a}]
}

In \code{GetProbableGenotypes}, or \code{GetWeightedMeanGenotypes}
when there are multiple inheritance modes and \code{onePloidyPerAllele = TRUE},
the genotype is simply the one corresponding to the inheritance mode with the minimum
\eqn{\chi ^2}{chi-squared} value:

\deqn{g_{t,a} = g_{t,a,m} | \min_m{\chi^2_{m,a}}}
}

\value{
For \code{GetWeightedMeanGenotypes},
a named matrix, with taxa in rows and alleles in columns,
and values ranging from \code{minval} to \code{maxval}.
These values can be treated as continuous genotypes.

For \code{GetProbableGenotypes}, a list:

\item{genotypes}{
A named integer matrix, with taxa in rows and alleles in columns,
and values ranging from zero to the maximum ploidy for
each allele.  These values can be treated as discrete genotypes.
}
\item{ploidy_index}{
A vector with one value per allele.  It contains the index
of the most inheritance mode of that allele in 
\code{object$priorProbPloidies}.
}
}
\author{
Lindsay V. Clark
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
exampleRAD_mapping <- AddGenotypePosteriorProb(exampleRAD_mapping)
exampleRAD_mapping <- AddPloidyChiSq(exampleRAD_mapping)

# get weighted mean genotypes
wmg <- GetWeightedMeanGenotypes(exampleRAD_mapping)
# examine the results
wmg[1:10,]

# get most probable genotypes
pg <- GetProbableGenotypes(exampleRAD_mapping, naIfZeroReads = TRUE)
# examine the results
pg$genotypes[1:10,]
}

\keyword{ methods }
\keyword{ array }
