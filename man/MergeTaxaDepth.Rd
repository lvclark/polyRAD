\name{MergeTaxaDepth}
\alias{MergeTaxaDepth}
\alias{MergeTaxaDepth.RADdata}
\title{
Combine Read Depths from Multiple Taxa into One Taxon
}
\description{
This function should be used in situations where data that were imported as
separate taxa should be merged into a single taxon.  The function should be
used before any of the pipeline functions for genotype calling.  Read depths
are summed across duplicate taxa and output as a single taxon.
}
\usage{
MergeTaxaDepth(object, ...)

\method{MergeTaxaDepth}{RADdata}(object, taxa, \dots)
}
\arguments{
  \item{object}{
A \code{\link{RADdata}} object.
}
\item{taxa}{
A character vector indicating taxa to be merged.  The first taxon in the vector
will be used to name the combined taxon in the output.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
Examples of reasons to use this function:

\itemize{
\item Duplicate samples across different libraries were given different
names so that preliminary analysis could confirm that they were truly the same
(\emph{i.e.} no mix-ups) before combining them.
\item Typos in the key file for the SNP mining software (TASSEL, Stacks, etc.)
caused duplicate samples to have different names when they really should have
had the same name.
}

To merge multiple sets of taxa into multiple combined taxa, this function can
be run multiple times or in a loop.
}
\value{
A \code{RADdata} object derived from \code{object}.  The \code{alleleDepth},
\code{antiAlleleDepth}, \code{locDepth}, \code{depthRatio}, and
\code{depthSamplingPermutation} slots, and \code{"taxa"} and \code{"nTaxa"}
attributes, have been changed accordingly to reflect the merge.
}
\author{
Lindsay V. Clark
}
\seealso{
\code{\link{SubsetByTaxon}}
}
\examples{
# dataset for this example
data(exampleRAD)

# merge the first three taxa into one
exampleRADm <- MergeTaxaDepth(exampleRAD, c("sample001", "sample002", "sample003"))

# inspect read depth
exampleRAD$alleleDepth[1:3,]
exampleRADm$alleleDepth[1:3,]
}
\keyword{ methods }
\keyword{ manip }
