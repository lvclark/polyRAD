\name{SubsetByLocus}
\alias{SubsetByLocus}
\alias{SubsetByLocus.RADdata}
\alias{SplitByChromosome}
\alias{SplitByChromosome.RADdata}
\alias{RemoveMonomorphicLoci}
\alias{RemoveMonomorphicLoci.RADdata}
\alias{RemoveHighDepthLoci}
\alias{RemoveHighDepthLoci.RADdata}
\alias{RemoveUngenotypedLoci}
\alias{RemoveUngenotypedLoci.RADdata}
\title{
Create RADdata Objects with a Subset of Loci
}
\description{
These functions take a \code{\link{RADdata}} object as input and generate smaller \code{RADdata}
objects containing only the specified loci.  \code{SubsetByLocus} allows the
user to specify which loci are kept, whereas \code{SplitByChromosome} creates
multiple \code{RADdata} objects representing chromosomes or sets of chromosomes.
\code{RemoveMonomorphicLoci} eliminates any loci with fewer than two alleles.
\code{RemoveHighDepthLoci} eliminates loci that have especially high read
depth in order to eliminate false loci originating from repetitive sequence.
\code{RemoveUngenotypedLoci} is intended for datasets that have been run
through \code{\link{PipelineMapping2Parents}} and may have some genotypes that
are missing or non-variable due to how priors were determined.
}
\usage{
SubsetByLocus(object, ...)
\method{SubsetByLocus}{RADdata}(object, loci, \dots)

SplitByChromosome(object, ...)
\method{SplitByChromosome}{RADdata}(object, chromlist = NULL, chromlist.use.regex = FALSE, 
                  fileprefix = "splitRADdata", \dots)
                  
RemoveMonomorphicLoci(object, ...)
\method{RemoveMonomorphicLoci}{RADdata}(object, verbose = TRUE, \dots)

RemoveHighDepthLoci(object, ...)
\method{RemoveHighDepthLoci}{RADdata}(object, max.SD.above.mean = 2, verbose = TRUE, \dots)

RemoveUngenotypedLoci(object, ...)
\method{RemoveUngenotypedLoci}{RADdata}(object, removeNonvariant = TRUE, \dots)
}
\arguments{
  \item{object}{
A \code{RADdata} object.
}
  \item{loci}{
A character or numeric vector indicating which loci to include in the output
\code{RADdata} object.  If numeric, it refers to row numbers in 
\code{object$locTable}.  If character, it refers to row names in 
\code{object$locTable}.
}
  \item{chromlist}{
An optional list indicating how chromosomes should be split into separate
\code{RADdata} objects.  Each item in the list is a vector of the same 
class as \code{object$locTable$Chr} (character or numeric) containing the
names of chromosomes that should go into one group.  If not provided, each
chromosome will be sent to a separate \code{RADdata} object.
}
  \item{chromlist.use.regex}{
If \code{TRUE}, the character strings in \code{chromlist} will be treated
as regular expressions for searching chromosome names.  For example, if
one wanted all chromosomes beginning with the string "scaffold" to go 
into one \code{RADdata} object, one could include the string
\code{"^scaffold"} as an item in \code{chromlist} and set 
\code{chromlist.use.regex = TRUE}.  If \code{FALSE}, exact matches to
chromosome names will be used.
}
  \item{fileprefix}{
A character string indicating the prefix of .RData files to export.  
}
\item{max.SD.above.mean}{
The maximum number of standard deviations above the mean read depth that a
locus can be in order to be retained.
}
\item{verbose}{
If \code{TRUE}, print out information about the original number of loci and
the number of loci that were retained.  For \code{RemoveHighDepthLoci}, a
histogram is also plotted showing mean depth per locus, and the cutoff for
removing loci.
}
\item{removeNonvariant}{
If \code{TRUE}, in addition to removing loci where posterior probabilities are
missing, loci will be removed where posterior probabilities are uniform across
the population.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
\code{SubsetByLocus} may be useful if the user has used their own filtering
criteria to determine a set of loci to retain, and wants to create a new
dataset with only those loci.  It can be used at any point in the analysis
process.

\code{SplitByChromosome} is intended to make large datasets more manageable
by breaking them into smaller datasets that can be processed independently,
either in parallel computing jobs on a cluster, or one after another on a
computer with limited RAM.  Generally it should be used immediately after
data import.  Rather than returning new \code{RADdata} objects, it saves
them individually to separate workspace image files, which can than be 
loaded one at a time to run analysis pipelines such as \code{\link{IteratePopStruct}}.
\code{\link{GetWeightedMeanGenotypes}} or one of the export functions can be 
run on each resulting \code{RADdata} object, and the resulting matrices
concatenated with \code{cbind}.

\code{SplitByChromosome}, \code{RemoveMonomorphicLoci}, and 
\code{RemoveHighDepthLoci} use \code{SubsetByLocus} internally.
}
\value{
\code{SubsetByLocus}, \code{RemoveMonomorphicLoci}, 
\code{RemoveHighDepthLoci}, and \code{RemoveUngenotypedLoci} return a \code{RADdata}
object with all the slots and attributes of \code{object}, but only 
containing the loci listed in \code{loci}, only loci with two or more
alleles, only loci without abnormally high depth, or only loci where posterior
probabilities are non-missing and variable, respectively.

\code{SplitByChromosome} returns a character vector containing file names
where .RData files have been saved.  Each .RData file contains one 
\code{RADdata} object named \code{splitRADdata}.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{VCF2RADdata}}, \code{\link{SubsetByTaxon}}
}

\examples{
# load a dataset for this example
data(exampleRAD)
exampleRAD

# just keep the first and fourth locus
subsetRAD <- SubsetByLocus(exampleRAD, c(1, 4))
subsetRAD

# split by groups of chromosomes
exampleRAD$locTable
tf <- tempfile()
splitfiles <- SplitByChromosome(exampleRAD, list(c(1, 4), c(6, 9)),
                                fileprefix = tf)
load(splitfiles[1])
splitRADdata

# filter out monomorphic loci (none removed in example)
filterRAD <- RemoveMonomorphicLoci(exampleRAD)

# filter out high depth loci (none removed in this example)
filterRAD2 <- RemoveHighDepthLoci(filterRAD)

# filter out loci with missing or non-variable genotypes 
# (none removed in this example)
filterRAD3 <- IterateHWE(filterRAD2)
filterRAD3 <- RemoveUngenotypedLoci(filterRAD3)
}

\keyword{ manip }
