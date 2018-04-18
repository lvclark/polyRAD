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
                 nPcsInit = 10, minfreq = 0.0001)
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
\item{nPcsInit}{
An integer indicating the number of principal component axes to initially 
estimate from \code{object$depthRatio}.  Passed to \code{\link{AddPCA}}.
}
\item{minfreq}{
A number indicating the minimum allele frequency allowed.  Passed to
\code{\link{AddAlleleFreqByTaxa}}.
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
\code{\link{AddGenotypeLikelihood}}, \code{\link{AddPloidyChiSq}}, and
\code{\link{AddGenotypePosteriorProb}}.
After the first PCA analysis, the number of principal component axes is not
allowed to decrease, and can only increase by one from one round to the next,
in order to help the algorithm converge.
}

\value{
A \code{"RADdata"} object identical to that passed to the function, but with
\code{$alleleFreq}, \code{$priorProb}, \code{$genotypeLikelihood},
and \code{$posteriorProb} slots added.  For \code{IterateHWE}, 
\code{$ploidyChiSq}, is also added.  For \code{IteratePopStruct}, 
\code{$alleleFreqByTaxa} is also added.
}

\note{
If you see the error 

\code{Error in if (rel_ch < threshold & count > 5) \{ : 
  missing value where TRUE/FALSE needed}
  
try lowering \code{nPcsInit}.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{GetWeightedMeanGenotypes}} for outputting genotypes in a 
useful format after iteration is completed.

\code{\link{StripDown}} to remove memory-hogging slots that are no longer
needed after the pipeline has been run.

\code{\link{PipelineMapping2Parents}} for mapping populations.
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

# run pipeline
exampleRAD <- IteratePopStruct(exampleRAD, nPcsInit = 3)

# export results
GetWeightedMeanGenotypes(exampleRAD)
}

\keyword{ iteration }

