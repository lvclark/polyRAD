\name{HindHe}
\alias{HindHe}
\alias{HindHe.RADdata}
\alias{HindHeMapping}
\alias{HindHeMapping.RADdata}
\title{
Identify Non-Mendelian Loci and Taxa that Deviate from Ploidy Expectations
}
\description{
\code{HindHe} and \code{HindHeMapping} both generate a matrix of values, with
taxa in rows and loci in columns.  The mean value of the matrix is expected to
be a certain value depending on the ploidy and, in the case of natural
populations and diversity panels, the inbreeding coefficient.  \code{colMeans}
of the matrix can be used to filter non-Mendelian loci from the dataset.
\code{rowMeans} of the matrix can be used to identify taxa that are not the
expected ploidy, are interspecific hybrids, or are a mix of multiple samples.
}
\usage{
HindHe(object, ...)

\method{HindHe}{RADdata}(object, omitTaxa = GetBlankTaxa(object), \dots)

HindHeMapping(object, ...)

\method{HindHeMapping}{RADdata}(object, n.gen.backcrossing = 0, n.gen.intermating = 0,
              n.gen.selfing = 0, ploidy = object$possiblePloidies[[1]],
              minLikelihoodRatio = 10,
              omitTaxa = c(GetDonorParent(object), GetRecurrentParent(object), 
                           GetBlankTaxa(object)), \dots)
}
\arguments{
  \item{object}{
A \code{\link{RADdata}} object.  Genotype calling does not need to have been
performed yet.  If the population is a mapping population,
\code{\link{SetDonorParent}} and \code{\link{SetRecurrentParent}} should have
been run already.
}
  \item{omitTaxa}{
A character vector indicating names of taxa not to be included in the output.
For \code{HindHe}, these taxa will also be omitted from allele frequency
estimations.
}
  \item{n.gen.backcrossing}{
The number of generations of backcrossing performed in a mapping population.  
}
  \item{n.gen.intermating}{
The number of generations of intermating performed in a mapping population.
Included for consistency with \code{\link{PipelineMapping2Parents}}, but
currently will give an error if set to any value other than zero.  If the most
recent generation in your mapping population was random mating among all
progeny, use \code{HindHe} instead of \code{HindHeMapping}.
}
  \item{n.gen.selfing}{
The number of generations of self-fertilization performed in a mapping population.
}
  \item{ploidy}{
A single value indicating the assumed ploidy to test.  Currently, only
autopolyploid and diploid inheritance modes are supported.
}
  \item{minLikelihoodRatio}{
Used internally by \code{\link{EstimateParentalGenotypes}} as a threshold for
certainty of parental genotypes.  Decrease this value if too many markers are
being discarded from the calculation.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
These functions are especially useful for highly duplicated genomes, in which
RAD tag alignments may have been incorrect, resulting in groups of alleles that
do not represent true Mendelian loci.  The statistic that is calculated is based
on the principle that observed heterozygosity will be higher than expected
heterozygosity if a "locus" actually represents two or more collapsed paralogs.
However, the statistic uses read depth in place of genotypes, eliminating the
need to perform genotype calling before filtering.

For a given taxon * locus, \eqn{H_{ind}}{Hind} is the probability that two
sequencing reads, sampled without replacement, are different alleles (RAD tags).

In \code{HindHe}, \eqn{H_E}{He} is the expected heterozygosity, estimated from
allele frequencies by taking the column means of \code{object$depthRatios}.
This is also the estimated probability that if two alleles were sampled at
random from the population at a given locus, they would be different alleles.

In \code{HindHeMapping}, \eqn{H_E}{He} is the average probability that in
a random progeny, two alleles sampled without replacement would be different.
The number of generations of backcrossing and self-fertilization, along with the
ploidy and estimated parental genotypes, are needed to make this calculation.
The function essentially simulates the mapping population based on parental
genotypes to determine \eqn{H_E}{He}.

The expectation is that

\deqn{H_{ind}/H_E = \frac{ploidy - 1}{ploidy} * (1 - F)}{Hind/He = (ploidy - 1)/ploidy * (1 - F)}

in a diversity panel, where \eqn{F} is the inbreeding coefficient, and

\deqn{H_{ind}/H_E = \frac{ploidy - 1}{ploidy}}{Hind/He = (ploidy - 1)/ploidy}

in a mapping population.  Loci that have much higher average values likely
represent collapsed paralogs that should be removed from the dataset.  Taxa
with much higher average values may be higher ploidy than expected,
interspecific hybrids, or multiple samples mixed together.
}
\value{
A named matrix, with taxa in rows and loci in columns.  For \code{HindHeMapping},
loci are omitted if consistent parental genotypes could not be determined across
alleles.
}
\references{
A preprint is in the works.  In the meantime, see a seminar describing
\eqn{H_{ind}/H_E}{Hind/He} at \url{https://youtu.be/Z2xwLQYc8OA?t=1678}.
}
\author{
Lindsay V. Clark
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(exampleRAD)

hhmat <- HindHe(exampleRAD)
colMeans(hhmat, na.rm = TRUE) # near 0.5 for diploid loci, 0.75 for tetraploid loci

data(exampleRAD_mapping)
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")

hhmat2 <- HindHeMapping(exampleRAD_mapping, n.gen.backcrossing = 1)
colMeans(hhmat2, na.rm = TRUE) # near 0.5; all loci diploid
}

\keyword{ distribution }
