\name{RADdata}
\alias{RADdata}
\title{
RADdata object constructor
}
\description{
\code{RADdata} is used internally to generate objects of the S3 class 
\dQuote{RADdata} by \pkg{polyRAD} functions for importing read depth data.
It is also available at the user level for cases where the data for import are
not already in a format supported by \pkg{polyRAD}.
}
\usage{
RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies, contamRate,
        alleleNucleotides)
}
\arguments{
  \item{alleleDepth}{
An integer matrix, with taxa in rows and alleles in columns.  Taxa names should
be included as row names.  Each value indicates the number of reads for a given
allele in a given taxon.  There should be no \code{NA} values; use zero to 
indicate no reads.
}
  \item{alleles2loc}{
An integer vector with one value for each column of \code{alleleDepth}.  The 
number indicates the identity of the locus to which the allele belongs.  A 
locus can have any number of alleles assigned to it (including zero).
}
  \item{locTable}{
A data frame, where locus names are row names.  There must be at least as 
many rows as the highest value of \code{alleles2loc}; each number in 
\code{alleles2loc} corresponds to a row index in \code{locTable}.  No columns
are required, although if provided a column named \dQuote{Chr} will be used for
indicating chromosome identities and a column named \dQuote{Pos} will be used
for indicating physical position.
}
  \item{possiblePloidies}{
A list, where each item in the list is an integer vector (or a numeric vector
that can be converted to integer).  Each vector indicates an inheritance 
pattern that SNPs in the dataset might obey.  \code{2} indicates diploid,
\code{4} indicates autotetraploid, \code{c(2, 2)} indicates allotetraploid,
\emph{etc.}
}
  \item{contamRate}{
A number ranging from zero to one (although in practice probably less than 
0.01) indicating the expected sample cross-contamination rate.
}
  \item{alleleNucleotides}{
A character vector or \code{DNAStringSet} with one value for each column of 
\code{alleleDepth}, indicating the DNA sequence for that allele.  Typically 
only the sequence at variable sites is provided, although intervening 
non-variable sequence can also be provided.
}
}
\details{
If \pkg{Biostrings} is installed, the \code{alleleNucleotides} slot will be
converted to a \code{DNAStringSet}.  Otherwise it will be left as a character
vector.
}
\value{
An object of the S3 class \dQuote{RADdata}.  The following slots are available
using the \code{$} operator:

\item{alleleDepth }{Identical to the argument provided to the function.}
\item{alleles2loc }{Identical to the argument provided to the function.}
\item{locTable }{Identical to the argument provided to the function.}
\item{possiblePloidies }{The \code{possiblePloidies} argument, converted to
integer.}
\item{locDepth }{A matrix with taxa in rows and loci in columns, with read 
depth summed across all alleles for each locus.  Column names are locus 
numbers rather than locus names.  See \code{\link{GetLocDepth}} for retrieving
the same matrix but with locus names as column names.}
\item{depthSamplingPermutations }{An integer matrix with taxa in rows and 
alleles in columns.  It is calculated as \eqn{locDepth choose alleleDepth}.
This is used as a coefficient for likelihood estimations done by other
\pkg{polyRAD} functions (\emph{i.e.} \code{\link{AddGenotypeLikelihood}}).}
\item{depthRatio }{A numeric matrix with taxa in rows and alleles in columns.
Calculated as \eqn{alleleDepth / locDepth}.  Used by other \pkg{polyRAD}
functions for rough estimation of genotypes and allele frequency.}
\item{antiAlleleDepth }{An integer matrix with taxa in rows and alleles in
columns.  For each allele, the number of reads from the locus that do NOT
belong to that allele.  Calculated as \eqn{locDepth - alleleDepth}. Used for
likelihood estimations by other \pkg{polyRAD} functions.}
\item{alleleNucleotides }{Identical to the argument provided to the 
function.}

The object additionally has several attributes (see \code{\link{attr}}):

\item{taxa }{A character vector listing all taxa names, in the same order as
the rows of \code{alleleDepth}.}
\item{nTaxa }{An integer indicating the number of taxa.}
\item{nLoc }{An integer indicating the number of loci in \code{locTable}.}
\item{contamRate }{Identical to the argument provided to the function.}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{readHMC}}, \code{\link{readTagDigger}}
}
\examples{
# create the dataset
mydepth <- matrix(sample(100, 16), nrow = 4, ncol = 4,
                  dimnames = list(paste("taxon", 1:4, sep = ""),
                  paste("loc", c(1,1,2,2), "_", c(0,1,0,1), sep = "")))
mydata <- RADdata(mydepth, c(1L,1L,2L,2L), 
                  data.frame(row.names = c("loc1", "loc2"), Chr = c(1,1),
                             Pos = c(2000456, 5479880)),
                  list(2, c(2,2)), 0.001, c("A", "G", "G", "T"))

# inspect the dataset
mydata
mydata$alleleDepth
mydata$locDepth
mydata$depthRatio

# the S3 class structure is flexible; other data can be added
mydata$GPS <- data.frame(row.names = attr(mydata, "taxa"),
                         Lat = c(43.12, 43.40, 43.05, 43.27),
                         Long = -c(70.85, 70.77, 70.91, 70.95))
mydata$GPS
}

\keyword{ methods }

