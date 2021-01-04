\name{AddGenotypePriorProb_HWE}
\alias{AddGenotypePriorProb_HWE}
\alias{AddGenotypePriorProb_HWE.RADdata}
\title{
Estimate Genotype Prior Probabilities In the Absence of Population Structure
}
\description{
Assuming Hardy-Weinberg Equilibrium, this function uses allele frequencies
and possible ploidies stored in a \dQuote{\link{RADdata}} object to estimate
genotype frequencies in the population, then stores these genotype 
frequencies in the \code{$priorProb} slot.  Inbreeding can also be simulated
using the \code{selfing.rate} argument.
}
\usage{
AddGenotypePriorProb_HWE(object, ...)
\method{AddGenotypePriorProb_HWE}{RADdata}(object, selfing.rate = 0, \dots)
}
\arguments{
  \item{object}{
A \dQuote{RADdata} object that has had allele frequencies added with
\code{\link{AddAlleleFreqHWE}}.
}
\item{selfing.rate}{
A number ranging from zero to one indicating the frequency of
self-fertilization in the species.
}
  \item{\dots}{
Additional arguments (none currently implemented).
}
}

\details{
For an autopolyploid, or within one subgenome of an allopolyploid,
genotype prior probabilities are estimated as:

\deqn{P(G_i) = {k \choose i} p^i * (1 - p)^{k - i}}{P(G_i) = (k choose i) * p^i * (1 - p)^(k - i)}

where \eqn{k} is the ploidy, \eqn{i} is the copy number of a given allele,
and \eqn{p} is the allele frequency in the population.

If the selfing rate is above zero, genotype prior probabilities are asjusted
according to Equation 6 of de Silva et al. (2005):

\deqn{P(G_{self}) = (1 - s)(I - sA)^{-1}P(G)}{P(G_{self}) = (1 - s) * (I - sA)^{-1}P(G)}

where \eqn{s} is the selfing rate.  \eqn{A} is a \eqn{k + 1 \times k + 1}{k + 1 x k + 1} matrix,
with each column representing the allele copy number from 0 to \eqn{k} of a
parental genotype, and each row representing the allele copy number from 0 to
\eqn{k} of a progeny genotype, and matrix elements representing the frequencies
of progeny after self-fertilization (each column summing to one).
}

\value{
A \dQuote{RADdata} object identical that passed to the function, but with
data stored in two new slots:

\item{priorProb }{A list of matrices, with one matrix per possible ploidy.  
For each matrix, allele copy number (from zero to the total
ploidy) is in rows, and alleles are in columns.  Each value is the probability
of sampling an individual with that allele copy number from the population.}
\item{priorProbPloidies }{A list identical to \code{object$possiblePloidies}.
It is in the same order as \code{$priorProb}, with each item indicating the
inheritance mode for the corresponding prior probability matrix.}
}

\author{
Lindsay V. Clark
}

\references{
De Silva, H. N., Hall, A. J., Rikkerink, E., and Fraser, L. G. (2005)
Estimation of allele frequencies in polyploids under
certain patterns of inheritance.  \emph{Heredity} \bold{95}, 327--334.
\doi{10.1038/sj.hdy.6800728}
}

\seealso{
\code{\link{AddGenotypePriorProb_Mapping2Parents}}, \code{\link{AddGenotypeLikelihood}},
\code{\link{AddGenotypePriorProb_ByTaxa}}
}
\examples{
# load in an example dataset
data(exampleRAD)
# add allele frequencies
exampleRAD <- AddAlleleFreqHWE(exampleRAD)
# add inheritance modes
exampleRAD$possiblePloidies <- list(2L, 4L, c(2L, 2L))

# estimate genotype prior probabilities
exampleRAD <- AddGenotypePriorProb_HWE(exampleRAD)

# examine results
exampleRAD$alleleFreq
exampleRAD$priorProb

# try it with inbreeding
exampleRAD2 <- AddGenotypePriorProb_HWE(exampleRAD, selfing.rate = 0.5)
exampleRAD2$priorProb
}

\keyword{ methods }
\keyword{ distribution }
