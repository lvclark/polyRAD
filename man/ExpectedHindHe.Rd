\name{ExpectedHindHe}
\alias{ExpectedHindHe}
\alias{SimGenotypes}
\alias{SimAlleleDepth}
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
\eqn{H_{ind}/H_E}{Hind/He} for each simulated locus. \code{SimGenotypes} and
\code{SimAlleleDepth} are internal functions used by \code{ExpectedHindHe}
but are provided at the user level since they may be more broadly useful.
}
\usage{
ExpectedHindHe(object, ploidy = object$possiblePloidies[[1]], inbreeding = 0,
               overdispersion = 20, reps = ceiling(5000/nLoci(object)),
               quiet = FALSE, plot = TRUE)

SimGenotypes(alleleFreq, alleles2loc, nsam, inbreeding, ploidy)

SimAlleleDepth(locDepth, genotypes, alleles2loc, overdispersion = 20)
}
\arguments{
  \item{object}{
A \code{RADdata} object.
}
  \item{ploidy}{
A single integer indicating the ploidy to use for genotype simulation.
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
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{ExpectedHindHe} invisibly returns a matrix, with loci in rows and reps in
columns, containing \eqn{H_{ind}/H_E}{Hind/He} from the simulated loci.

\code{SimGenotypes} returns a numeric matrix of allele copy number, with samples
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
}

\keyword{ datagen }
