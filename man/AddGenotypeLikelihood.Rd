\name{AddGenotypeLikelihood}
\alias{AddGenotypeLikelihood}
\alias{AddGenotypeLikelihood.RADdata}
\alias{AddDepthSamplingPermutations}
\title{
Estimate Genotype Likelihoods in a RADdata object
}
\description{
For each possible allele copy number across each possible ploidy in each taxon,
\code{AddGenotypeLikelihood} estimates the probability of observing the 
distribution of read counts that are recorded for that taxon and locus.
\code{AddDepthSamplingPermutations} is called by \code{AddGenotypeLikelihood}
the first time it is run, so that part of the likelihood calcluation is
stored in the \code{RADdata} object and doesn't need to be re-run on each
iteration of the \link[=IterateHWE]{pipeline} functions.
}
\usage{
AddGenotypeLikelihood(object, ...)

\method{AddGenotypeLikelihood}{RADdata}(object, overdispersion = 9, \dots)

AddDepthSamplingPermutations(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{overdispersion}{
An overdispersion parameter.  Higher values will cause the expected read depth
distribution to more resemble the binomial distribution.  Lower values indicate
more overdispersion, \emph{i.e.} sample-to-sample variance in the probability
of observing reads from a given allele.
}
  \item{\dots}{
Other arguments; none are currently used.
}
}
\details{
If allele frequencies are not already recorded in \code{object}, they will 
be added using \code{\link{AddAlleleFreqHWE}}.  Allele frequencies are then
used for estimating the probability of sampling an allele from a genotype due
to sample contamination.  Given a known genotype with \eqn{x} copies of 
allele \eqn{i}, ploidy \eqn{k}, allele frequency \eqn{p_i} in the population used for
making sequencing libraries, and contamination rate \eqn{c}, the probabiity of 
sampling a read \eqn{r_i} from that locus corresponding to that allele is

\deqn{P(r_i | x) = \frac{x}{k} * (1 - c) + p_i * c}{P(r_i | x) = x/k * (1 - c) + p_i * c}

To estimate the genotype likelihood, where \eqn{nr_i} is the number of reads 
corresponding to allele \eqn{i} for a given taxon and locus and \eqn{nr_j} is the
number of reads corresponding to all other alleles for that taxon and locus:

\deqn{P(nr_i, nr_j | x) = {{nr_i + nr_j}\choose{nr_i}} * \frac{B[P(r_i | x) * d + nr_i, [1 - P(r_i | x)] * d + nr_j]]}{B[P(r_i | x) * d, [1 - P(r_i | x)] * d]}}{
P(nr_i, nr_j | x) = [(nr_i + nr_j) choose nr_i] * B[P(r_i | x) * d + nr_i, [1 - P(r_i | x)] * d + nr_j]] / B[P(r_i | x) * d, [1 - P(r_i | x)] * d]
}

where

\deqn{{{nr_i + nr_j}\choose{nr_i}} = \frac{(nr_i + nr_j)!}{nr_i! * nr_j!}}{
[(nr_i + nr_j) choose nr_i] = (nr_i + nr_j)!/(nr_i! * nr_j!)
}

B is the beta function, and \eqn{d} is the overdispersion parameter set by 
\code{overdispersion}.  \eqn{{{nr_i + nr_j}\choose{nr_i}}}{[(nr_i + nr_j) choose nr_i]}
is calculated by \code{AddDepthSamplingPermutations}.
}
\value{
A \code{"RADdata"} object identical to that passed to the function, but with
genotype likelihoods stored in \code{object$genotypeLikelihood}.  This item is a 
two dimensional list, with one row for each ploidy listed
in \code{object$possiblePloidies}, ignoring differences between 
autopolyploids and allopolyploids, and one column for each ploidy listed in
\code{object$taxaPloidy}.  Each item in the list is a three-dimensional 
array with number of allele copies in the first dimension, taxa in the second dimension,
and alleles in the third dimension.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddAlleleFreqMapping}}
}
\examples{
# load example dataset and add allele frequency
data(exampleRAD)
exampleRAD <- AddAlleleFreqHWE(exampleRAD)

# estimate genotype likelihoods
exampleRAD <- AddGenotypeLikelihood(exampleRAD)

# inspect the results
# the first ten individuals and first two alleles, assuming diploidy
exampleRAD$alleleDepth[1:10,1:2]
exampleRAD$genotypeLikelihood[[1]][,1:10,1:2]
# assuming tetraploidy
exampleRAD$genotypeLikelihood[[2]][,1:10,1:2]
}

\keyword{ methods }

