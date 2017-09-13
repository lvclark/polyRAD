\name{AddGenotypePriorProb_HWE}
\alias{AddGenotypePriorProb_HWE}
\title{
Estimate Genotype Prior Probabilities Under Hardy-Weinberg Equilibrium
}
\description{
Assuming Hardy-Weinberg Equilibrium, this function uses allele frequencies
and possible ploidies stored in a \dQuote{\link{RADdata}} object to estimate
genotype frequencies in the population, then stores these genotype 
frequencies in the \code{$priorProb} slot.
}
\usage{
AddGenotypePriorProb_HWE(object, ...)
}
\arguments{
  \item{object}{
A \dQuote{RADdata} object that has had allele frequencies added with
\code{\link{AddAlleleFreqHWE}}.
}
  \item{\dots}{
Additional arguments (none currently implemented).
}
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


\seealso{
\code{\link{AddGenotypePriorProb_Mapping2Parents}}, \code{\link{AddGenotypeLikelihood}}
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
}

\keyword{ methods }
\keyword{ distribution }
