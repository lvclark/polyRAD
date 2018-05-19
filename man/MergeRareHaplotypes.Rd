\name{MergeRareHaplotypes}
\alias{MergeRareHaplotypes}
\title{
Consolidate Reads from Rare Alleles
}
\description{
\code{MergeRareHaplotypes} searches for rare alleles in a
\code{"\link{RADdata}"} object, and merges them into the most similar allele
at the same locus based on nucleotide sequence (or the most common if multiple
are equally similar).  Read
depth is summed across merged alleles, and the \code{alleleNucleotides} slot
of the \code{"RADdata"} object contains IUPAC ambiguity codes to indicate 
nucleotide differences across merged alleles.  This function is designed to be
used immediately after data import.
}
\usage{
MergeRareHaplotypes(object, ...)
\method{MergeRareHaplotypes}{RADdata}(object, min.ind.with.haplotype = 10, \dots)

}
\arguments{
  \item{object}{
A \code{"RADdata"} object.
}
\item{min.ind.with.haplotype}{
The minimum number of taxa having reads from a given allele for that allele to
not be merged.
}
  \item{\dots}{
Additional arguments; none implemented.
}
}

\value{
A \code{"RADdata"} object identical to \code{object}, but with its 
\code{$alleleDepth}, \code{$antiAlleleDepth}, \code{$depthRatio},
\code{$depthSamplingPermutations}, \code{$alleleNucleotides}, and
\code{$alleles2loc} arguments adjusted after merging alleles.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{SubsetByLocus}}, \code{\link{VCF2RADdata}}, \code{\link{readStacks1}}
}
\examples{
data(exampleRAD)
exampleRAD2 <- MergeRareHaplotypes.RADdata(exampleRAD, 
                                           min.ind.with.haplotype = 12)
exampleRAD$alleleDepth[21:30,3:5]
exampleRAD2$alleleDepth[21:30,3:4]
exampleRAD$alleleNucleotides
exampleRAD2$alleleNucleotides
}

\keyword{ manip }
\keyword{ methods }
