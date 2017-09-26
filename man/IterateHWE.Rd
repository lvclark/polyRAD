\name{IterateHWE}
\alias{IterateHWE}
\alias{IteratePopStruct}
\title{
Iteratively Estimate Population Parameters and Genotypes In a Diversity Panel
}
\description{
These are wrapper function that iteratively run other \pkg{polyRAD} functions 
until allele frequencies stabilize to within a user-defined threshold.  Genotype
posterior probabilities can then be exported for downstream analysis.
}
\usage{
IterateHWE(object, tol = 1e-08, excludeTaxa = GetBlankTaxa(object))
IteratePopStruct(object, tol = 1e-03, excludeTaxa = GetBlankTaxa(object),
                 nPcsInit = 50)
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

\details{
For \code{IterateHWE}, the following functions are run iteratively,
assuming Hardy-Weinberg Equilibrium:
\code{\link{AddAlleleFreqHWE}}, 
\code{\link{AddGenotypePriorProb_HWE}}, \code{\link{AddGenotypeLikelihood}},
\code{\link{AddPloidyChiSq}}, and \code{\link{AddGenotypePosteriorProb}}.

For \code{IteratePopStruct}, the following functions are run iteratively,
modeling population structure:
\code{\link{AddPCA}}, \code{\link{AddAlleleFreqByTaxa}}, 
\code{\link{AddAlleleFreqHWE}}, \code{\link{AddGenotypePriorProb_ByTaxa}},
\code{\link{AddGenotypeLikelihood}}, and \code{\link{AddGenotypePosteriorProb}}.
}

\value{
A \code{"RADdata"} object identical to that passed to the function, but with
\code{$alleleFreq}, \code{$priorProb}, \code{$genotypeLikelihood},
and \code{$posteriorProb} slots added.  For \code{IterateHWE}, 
\code{$ploidyChiSq}, is also added.  For \code{IteratePopStruct}, 
\code{$alleleFreqByTaxa} is also added.
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

# re-load to run pipeline assuming population structure
data(exampleRAD)

# for now, IteratePopStruct only works with one ploidy
exampleRAD$possiblePloidies <- list(2L)

# run pipeline
exampleRAD <- IteratePopStruct(exampleRAD, nPcsInit = 3)

# export results
GetWeightedMeanGenotypes(exampleRAD)
}

\keyword{ iteration }

