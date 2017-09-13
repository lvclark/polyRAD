\name{IterateHWE}
\alias{IterateHWE}
\title{
Iteratively Estimate Population Parameters and Genotypes Under Hardy-Weinberg Equilibrium
}
\description{
This is a wrapper function that iteratively runs \code{\link{AddAlleleFreqHWE}}, 
\code{\link{AddGenotypePriorProb_HWE}}, \code{\link{AddGenotypeLikelihood}},
\code{\link{AddPloidyChiSq}}, and \code{\link{AddGenotypePosteriorProb}} until
allele frequencies stabilize to within a user-defined threshold.  Genotype
posterior probabilities can then be exported for downstream analysis.
}
\usage{
IterateHWE(object, tol = 1e-08, excludeTaxa = GetBlankTaxa(object))
}

\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
  \item{tol}{
A number indicating when the iteration should end.  It indicates the maximum mean 
difference in allele frequencies between iterations that is tolerated.  Larger
numbers will lead to fewer iterations.
}
  \item{excludeTaxa}{
A character vector indicating names of taxa that should be excluded from allele 
frequency estimates and chi-squared estimates.
}
}

\value{
A \code{"RADdata"} object identical to that passed to the function, but with
\code{$alleleFreq}, \code{$priorProb}, \code{$genotypeLikelihood},
\code{$ploidyChiSq}, and \code{$posteriorProb} slots added.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{GetWeightedMeanGenotypes}} for outputting genotypes in a 
useful format after iteration is completed.
}
\examples{
# load dataset
data(exampleRAD)

# iteratively estimate parameters
exampleRAD <- IterateHWE(exampleRAD)

# export results
GetWeightedMeanGenotypes(exampleRAD)
}

\keyword{ iteration }

