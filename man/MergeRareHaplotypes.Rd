\name{MergeRareHaplotypes}
\alias{MergeRareHaplotypes}
\alias{MergeRareHaplotypes.RADdata}
\title{
Consolidate Reads from Rare Alleles
}
\description{
\code{MergeRareHaplotypes} searches for rare alleles in a
\code{"\link{RADdata}"} object, and merges them into the most similar allele
at the same locus based on nucleotide sequence (or the most common allele if
multiple are equally similar).  Read
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

\details{
Alleles with zero reads across the entire dataset are removed by
\code{MergeRareHaplotypes} without merging nucleotide sequences.  After
merging, at least one allele is left, even if it has fewer than 
\code{min.ind.with.haplotype} taxa with reads, as long as it has more than zero
taxa with reads.
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
\code{\link{SubsetByLocus}}, \code{\link{VCF2RADdata}}, \code{\link{readStacks}}
}
\examples{
data(exampleRAD)
exampleRAD2 <- MergeRareHaplotypes(exampleRAD, 
                                   min.ind.with.haplotype = 20)
exampleRAD$alleleDepth[21:30,6:7]
exampleRAD2$alleleDepth[21:30,6,drop=FALSE]
exampleRAD$alleleNucleotides
exampleRAD2$alleleNucleotides
}

\keyword{ manip }
\keyword{ methods }
