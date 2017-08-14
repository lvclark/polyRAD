\name{AddGenotypeLikelihood}
\alias{AddGenotypeLikelihood}
\title{
Estimate Genotype Likelihoods in a RADdata object
}
\description{
For each possible allele copy number across each possible ploidy in each taxon,
\code{AddGenotypeLikelihood} estimates the probability of observing the 
distribution of read counts that are recorded for that taxon and locus.
}
\usage{
AddGenotypeLikelihood(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
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

\deqn{P(nr_i, nr_j | x) = {{nr_i + nr_j}\choose{nr_i}}P(r_i | x)^{nr_i} * [1 - P(r_i | x)]^{nr_j}}{
P(nr_i, nr_j | x) = [(nr_i + nr_j) choose nr_i] * P(r_i | x)^{nr_i} * [1 - P(r_i | x)]^{nr_j}
}
}
\value{
A \code{"RADdata"} object identical to that passed to the function, but with
genotype likelihoods stored in \code{object$genotypeLikelihood}.  This item is a 
list, with one item for each possible ploidy, ignoring differences between 
autopolyploids and allopolyploids.  For each ploidy there is a three-dimensional 
array with number of allele copies in the first dimension, taxa in the second dimension,
and alleles in the third dimension.
}

\author{
Lindsay V. Clark
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    UseMethod("AddGenotypeLikelihood", object)
  }
}

\keyword{ methods }

