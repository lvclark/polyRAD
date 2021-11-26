\name{GetLikelyGen}
\alias{GetLikelyGen}
\title{
Output the Most Likely Genotype
}
\description{
For a single taxon in a \code{"\link{RADdata}"} object, \code{GetLikelyGen}
returns the most likely genotype (expressed in allele copy number) for each
allele and each possible ploidy.  The likelihoods used for determining 
genotypes are those stored in \code{object$genotypeLikelihood}.
}
\usage{
GetLikelyGen(object, taxon, minLikelihoodRatio = 10)
}
\arguments{
  \item{object}{
A \code{"RADdata"} object.
}
  \item{taxon}{
A character string indicating the taxon for which genotypes should be
returned.
}
\item{minLikelihoodRatio}{
A number indicating the minimum ratio of the likelihood of the most likely
genotype to the likelihood of the second-most likely genotype for any
genotype to be output for a given allele.  If this number is one or less,
all of the most likely genotypes will be output regardless of likelihood 
ratio.  Where filtering is required so that only high confidence genotypes
are retained, this number should be increased.

}
}
\value{
A matrix with ploidies in rows (named with ploidies converted to character 
format) and alleles in columns.  Each value indicates the most likely number
of copies of that allele that the taxon has, assuming that ploidy.
}
\author{
Lindsay V. Clark
}


\seealso{
\code{\link{AddGenotypeLikelihood}}
}
\examples{
# load dataset for this example
data(exampleRAD)
# add allele frequencies and genotype likelihoods
exampleRAD <- AddAlleleFreqHWE(exampleRAD)
exampleRAD <- AddGenotypeLikelihood(exampleRAD)

# get most likely genotypes
GetLikelyGen(exampleRAD, "sample001")
GetLikelyGen(exampleRAD, "sample082")

# try different filtering
GetLikelyGen(exampleRAD, "sample001", minLikelihoodRatio = 1)
GetLikelyGen(exampleRAD, "sample001", minLikelihoodRatio = 100)
}

\keyword{ methods }

