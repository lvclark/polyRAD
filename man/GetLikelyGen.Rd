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
GetLikelyGen(object, taxon)
}
\arguments{
  \item{object}{
A \code{"RADdata"} object.
}
  \item{taxon}{
A character string indicating the taxon for which genotypes should be
returned.
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
GetLikelyGen(exampleRAD, "sample002")
}

\keyword{ methods }

