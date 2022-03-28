\name{SubsetByPloidy}
\alias{SubsetByPloidy}
\alias{SubsetByPloidy.RADdata}
\title{
Create a RADdata object with a Subset of Possible Ploidies
}
\description{
This function is used for removing some of the ploidies (i.e. inheritance modes
possible across loci)
stored in a \code{\link{RADdata}} object.  If genotype calling has already
been performed, all of the relevant slots will be subsetted to only keep the
ploidies that the user indicates.
}
\usage{
SubsetByPloidy(object, ...)
\method{SubsetByPloidy}{RADdata}(object, ploidies, \dots)
}
\arguments{
  \item{object}{
A \code{RADdata} object.
}
  \item{ploidies}{
A list, formatted like \code{object$possiblePloidies}, indicating ploidies to
retain.  Each item in the list is a vector, where \code{2} indicates diploid,
\code{c(2, 2)} allotetraploid, \code{4} autotetraploid, etc.
}
  \item{\dots}{
Other arguments (none implemented).
}
}
\details{
Note that slots of \code{object} are subsetted but not recalculated.  For
example, \code{\link{GetWeightedMeanGenotypes}} takes a weighted mean across
ploidies, which is in turn used for estimating allele frequencies and 
performing PCA.  If the values in \code{object$ploidyChiSq} are considerably
higher for the ploidies being removed than for the ploidies being retained,
this difference is likely to be small and not substantially impact genotype
calling.  Otherwise, it may be advisable to 
\link[=IterateHWE]{re-run genotype calling} after running \code{SubsetByPloidy}.
}
\value{
A \code{RADdata} object identical to \code{object}, but only containing data
relevant to the inheritance modes listed in \code{ploidies}.
}
\note{
If you only wish to retain taxa of a certain ploidy, instead do

\code{object <- SubsetByTaxon(object, GetTaxaByPloidy(object, 4))}

to, for example, only retain tetraploid taxa.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{SubsetByTaxon}}, \code{\link{SubsetByLocus}}
}
\examples{
# Example dataset assuming diploidy or autotetraploidy
data(exampleRAD)
exampleRAD <- IterateHWE(exampleRAD)
# Subset to only keep tetraploid results
exampleRAD <- SubsetByPloidy(exampleRAD, ploidies = list(4))
}
\keyword{ manip }
