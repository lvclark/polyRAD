\name{OneAllelePerMarker}
\alias{OneAllelePerMarker}
\alias{OneAllelePerMarker.RADdata}
\title{
Return the Index of One Allele for Each Locus
}
\description{
This function exists primarily to be called by functions such as 
\code{\link{AddPCA}} and \code{\link{GetWeightedMeanGenotypes}}
that may need to exclude one allele per locus to avoid mathematical
singularities.  For a \code{"\link{RADdata}"} object, it returns
the indices of one allele per locus.
}
\usage{
OneAllelePerMarker(object, ...)
\method{OneAllelePerMarker}{RADdata}(object, commonAllele = FALSE, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
  \item{commonAllele}{
If \code{TRUE}, the index of the most common allele for each locus is 
returned, according to \code{object$alleleFreq}.  If \code{FALSE}, 
the index of the first allele for each locus is returned.
  }
  \item{\dots}{
Additional arguments (none implemented).
}
}
\value{
An integer vector indicating the index of one allele for each locus
in \code{object}.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{GetTaxa}} for a list of accessors.
}
\examples{
data(exampleRAD)

OneAllelePerMarker(exampleRAD)

OneAllelePerMarker(exampleRAD, commonAllele = TRUE)
}

\keyword{ utilities }
\keyword{ methods }
