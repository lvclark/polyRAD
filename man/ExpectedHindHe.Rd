\name{ExpectedHindHe}
\alias{ExpectedHindHe}
\alias{SimGenotypes}
\alias{SimAlleleDepth}
\alias{SimGenotypesMapping}
\alias{ExpectedHindHeMapping}
\title{
Simulate Data to Get Expected Distribution of Hind/He
}
\description{
These functions were created to help users determine an appropriate cutoff for
filtering loci based on \eqn{H_{ind}/H_E}{Hind/He} after running
\code{\link{HindHe}} and \code{\link{InbreedingFromHindHe}}.
\code{ExpectedHindHe} takes allele frequencies, sample size, and read depths from
a \code{\link{RADdata}} object, simulates genotypes and allelic read depths from
these assuming Mendelian inheritance, and then estimates
\eqn{H_{ind}/H_E}{Hind/He} for each simulated locus.
\code{ExpectedHindHeMapping} performs similar simulation and estimation, but
in mapping populations based on parental genotypes and expected distribution
of progeny genotypes.
\code{SimGenotypes}, \code{SimGenotypesMapping}, and
\code{SimAlleleDepth} are internal functions used by \code{ExpectedHindHe}
and \code{ExpectedHindHeMapping}
but are provided at the user level since they may be more broadly useful.
}
\usage{
ExpectedHindHe(object, ploidy = object$possiblePloidies[[1]], inbreeding = 0,
               overdispersion = 20, contamRate = 0, errorRate = 0.001,
               reps = ceiling(5000/nLoci(object)),
               quiet = FALSE, plot = TRUE)

ExpectedHindHeMapping(object, ploidy = object$possiblePloidies[[1]],
                      n.gen.backcrossing = 0, n.gen.selfing = 0,
                      overdispersion = 20, contamRate = 0, errorRate = 0.001,
                      freqAllowedDeviation = 0.05, 
                      minLikelihoodRatio = 10, reps = ceiling(5000/nLoci(object)),
                      quiet = FALSE, plot = TRUE)

SimGenotypes(alleleFreq, alleles2loc, nsam, inbreeding, ploidy)

SimGenotypesMapping(donorGen, recurGen, alleles2loc, nsam,
                    ploidy.don, ploidy.rec,
                    n.gen.backcrossing, n.gen.selfing)

SimAlleleDepth(locDepth, genotypes, alleles2loc, overdispersion = 20,
               contamRate = 0, errorRate = 0.001)
}
\arguments{
  \item{object}{
A \code{RADdata} object.
}
  \item{ploidy}{
A single integer indicating the ploidy to use for genotype simulation.
For \code{ExpectedHindHe} and \code{ExpectedHindHeMapping},
this number will be multiplied by the values in \code{GetTaxaPloidy(object)}
then divided by two to determine the ploidy of each individual for simulation.
}
  \item{inbreeding}{
A number ranging from 0 to 1 indicating the amount of inbreeding (\eqn{F}).
This represents inbreeding from all sources (population structure,
self-fertilization, etc.) and can be estimated with \code{InbreedingFromHindHe}.
}
  \item{overdispersion}{
Overdispersion parameter as described in \code{\link{AddGenotypeLikelihood}}.
Lower values will cause allelic read depth distributions to deviate further from
expectations based on allele copy number.
}
  \item{contamRate}{
Sample cross-contamination rate to simulate.  Although 0 is the default, 0.001
is also reasonable.
}
  \item{errorRate}{
Sequencing error rate to simulate. For Illumina reads, 0.001 is a reasonable
value. An error is assumed to have an equal chance of converting an allele to
any other allele at the locus, although this is somewhat of an
oversimplification.
}
  \item{reps}{
The number of times to simulate the data and estimate \eqn{H_{ind}/H_E}{Hind/He}.
This can generally be left at the default, but
set it higher than 1 if you want to see within-locus variance in the estimate.
}
  \item{quiet}{
Boolean indicating whether to suppress messages and results printed to console.
}
  \item{plot}{
Boolean indicating whether to plot a histogram of \eqn{H_{ind}/H_E}{Hind/He}
values.
}
\item{n.gen.backcrossing}{
An integer indicating the number of generations of backcrossing.
}
\item{n.gen.selfing}{
An integer indicating the number of generations of self-fertilization.
}
\item{freqAllowedDeviation}{
The amount by which allele frequencies are allowed to deviate from expected
allele frequencies.  See \code{\link{AddAlleleFreqMapping}}.
}
\item{minLikelihoodRatio}{
Minimum likelihood ratio for determining the most likely parental genotypes.
See \code{\link{GetLikelyGen}}.
}
\item{alleleFreq}{
A vector of allele frequencies, as can be found in the \code{$alleleFreq} slot
of a \code{RADdata} object after running \code{\link{AddAlleleFreqHWE}}.
}
\item{alleles2loc}{
An integer vector assigning alleles to loci, as can be found in the
\code{$alleles2loc} slot of a \code{RADdata} object.
}
\item{nsam}{
An integer indicating the number of samples (number of taxa) to simulate.
}
\item{donorGen}{
A vector indicating genotypes of the donor parent (which can be either parent
if backcrossing was not performed), with one value for each allele in the
dataset, and numbers indicating the copy number of each allele.
}
\item{recurGen}{
A vector indicating genotypes of the recurrent parent, as with \code{donorGen}.
}
\item{ploidy.don}{
A single integer indicating the ploidy of the donor parent.
}
\item{ploidy.rec}{
A single integer indicating the ploidy of the recurrent parent.
}
\item{locDepth}{
An integer matrix indicating read depth at each taxon and locus.  Formatted as
the \code{$locDepth} slot of a \code{RADdata} object, notably with columns
named by locus number rather than locus name.
}
\item{genotypes}{
A numeric matrix, formatted as the output of \code{\link{GetProbableGenotypes}}
or \code{SimGenotypes}, indicating genotypes as allele copy number.
}
}
\details{
To prevent highly inflated values in the output, \code{ExpectedHindHe} filters
loci with minor allele frequencies below five times the sequencing error rate.
}
\value{
\code{ExpectedHindHe} and \code{ExpectedHindHeMapping} invisibly return a
matrix, with loci in rows and reps in
columns, containing \eqn{H_{ind}/H_E}{Hind/He} from the simulated loci.

\code{SimGenotypes} and \code{SimGenotypesMapping} return a numeric matrix of
allele copy number, with samples
in rows and alleles in columns, similar to that produced by
\code{\link{GetProbableGenotypes}}.

\code{SimAlleleDepth} returns an integer matrix of allelic read depth, with
samples in rows and alleles in columns, similar to the \code{$alleleDepth}
slot of a \code{RADdata} object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Lindsay V. Clark
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# Load dataset for the example
data(exampleRAD)
exampleRAD <- AddAlleleFreqHWE(exampleRAD)

# Simulate genotypes
simgeno <- SimGenotypes(exampleRAD$alleleFreq, exampleRAD$alleles2loc, 10, 0.2, 2)

# Simulate reads
simreads <- SimAlleleDepth(exampleRAD$locDepth[1:10,], simgeno, exampleRAD$alleles2loc)

# Get expected Hind/He distribution if all loci in exampleRAD were well-behaved
ExpectedHindHe(exampleRAD, reps = 10)

# Mapping population example
data(exampleRAD_mapping)
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")
exampleRAD_mapping <- AddAlleleFreqMapping(exampleRAD_mapping,
                                           expectedFreqs = c(0.25, 0.75),
                                           allowedDeviation = 0.08)
exampleRAD_mapping <- AddGenotypeLikelihood(exampleRAD_mapping)
exampleRAD_mapping <- EstimateParentalGenotypes(exampleRAD_mapping,
                                                n.gen.backcrossing = 1)

simgenomap <- SimGenotypesMapping(exampleRAD_mapping$likelyGeno_donor[1,],
                                  exampleRAD_mapping$likelyGeno_recurrent[1,],
                                  exampleRAD_mapping$alleles2loc,
                                  nsam = 10, ploidy.don = 2, ploidy.rec = 2,
                                  n.gen.backcrossing = 1,
                                  n.gen.selfing = 0)
                                  
ExpectedHindHeMapping(exampleRAD_mapping, n.gen.backcrossing = 1, reps = 10)
}

\keyword{ datagen }
