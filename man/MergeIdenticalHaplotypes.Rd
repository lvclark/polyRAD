\name{MergeIdenticalHaplotypes}
\alias{MergeIdenticalHaplotypes}

\title{
Merge Alleles with Identical DNA Sequences
}
\description{
If any alleles within a locus have identical \code{alleleNucleotides} values,
this function merges those alleles, summing their read depths.  This function is
primarily intended to be used internally in cases where tags vary in length
within a locus, resulting in truncated \code{alleleNucleotides}.
}
\usage{
MergeIdenticalHaplotypes(object, ...)
}

\arguments{
  \item{object}{
A \code{\link{RADdata}} object.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
%\details{
%}
\value{
A \code{RADdata} object identical to \code{object}, but with alleles merged.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{MergeRareHaplotypes}}, \code{\link{readProcessIsoloci}}
}
\examples{
data(exampleRAD)
# change a haplotype for this example
exampleRAD$alleleNucleotides[5] <- "GC"

nAlleles(exampleRAD)
exampleRAD <- MergeIdenticalHaplotypes(exampleRAD)
nAlleles(exampleRAD)
}

\keyword{ manip }
