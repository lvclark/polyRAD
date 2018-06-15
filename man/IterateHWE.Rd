\name{IterateHWE}
\alias{IterateHWE}
\alias{IterateHWE_LD}
\alias{IteratePopStruct}
\alias{IteratePopStructLD}
\title{
Iteratively Estimate Population Parameters and Genotypes In a Diversity Panel
}
\description{
These are wrapper function that iteratively run other \pkg{polyRAD} functions 
until allele frequencies stabilize to within a user-defined threshold.  Genotype
posterior probabilities can then be exported for downstream analysis.
}
\usage{
IterateHWE(object, selfing.rate = 0, tol = 1e-05,
           excludeTaxa = GetBlankTaxa(object))
IterateHWE_LD(object, selfing.rate = 0, tol = 1e-05, 
              excludeTaxa = GetBlankTaxa(object),
              LDdist = 1e4, minLDcorr = 0.2)
IteratePopStruct(object, selfing.rate = 0, tol = 1e-03,
                 excludeTaxa = GetBlankTaxa(object),
                 nPcsInit = 10, minfreq = 0.0001)
IteratePopStructLD(object, selfing.rate = 0, tol = 1e-03,
                   excludeTaxa = GetBlankTaxa(object),
                   nPcsInit = 10, minfreq = 0.0001, LDdist = 1e4, 
                   minLDcorr = 0.2)
}

\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{selfing.rate}{
A number ranging from zero to one indicating the frequency of
self-fertilization in the species.
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
\item{LDdist}{
The distance, in basepairs, within which to search for alleles that may be in
linkage disequilibrium with a given allele.
}
\item{minLDcorr}{
The minimum correlation coefficient between two alleles 
for linkage disequilibrium between those alleles to be used by the pipeline
for genotype estimation; see \code{\link{AddAlleleLinkages}}.
}
}

\details{
For \code{IterateHWE}, the following functions are run iteratively,
assuming no population structure:
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

\code{IteratePopStructLD} runs each of the functions listed for 
\code{IteratePopStruct} once, then runs \code{\link{AddAlleleLinkages}}.
It then runs  
\code{\link{AddAlleleFreqHWE}}, \code{\link{AddGenotypePriorProb_ByTaxa}},
\code{\link{AddGenotypePriorProb_LD}},
\code{\link{AddGenotypeLikelihood}}, \code{\link{AddPloidyChiSq}}, 
\code{\link{AddGenotypePosteriorProb}}, \code{\link{AddPCA}}, 
and \code{\link{AddAlleleFreqByTaxa}} iteratively until convergence of
allele frequencies.
}

\value{
A \code{"RADdata"} object identical to that passed to the function, but with
\code{$alleleFreq}, \code{$priorProb}, \code{$genotypeLikelihood}, 
\code{$ploidyChiSq}, and \code{$posteriorProb} slots added.  
For \code{IteratePopStruct} and \code{IteratePopStructLD}, 
\code{$alleleFreqByTaxa} and \code{$PCA} are also added.  For 
\code{IteratePopStructLD} and \code{IterateHWE_LD}, \code{$alleleLinkages} 
and \code{$priorProbLD} are also added.
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

# dataset for LD pipeline
require(Biostrings)
data(Msi01genes)

# run HWE + LD pipeline
mydata1 <- IterateHWE_LD(Msi01genes)

# run pop. struct + LD pipeline
mydata2 <- IteratePopStructLD(Msi01genes)
}

\keyword{ iteration }

