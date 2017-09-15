\name{Accessors}
\alias{GetTaxa}
\alias{GetLoci}
\alias{GetLocDepth}
\alias{GetContamRate}
\alias{SetContamRate}
\alias{nTaxa}
\alias{nAlleles}
\alias{GetAlleleNames}
\title{
Accessor Functions for RADdata Objects
}
\description{
These functions can be used for accessing and replacing data within
a \code{"\link{RADdata}"} object.  Data slots that do not yet have
accessors can be accessed and replaced using the \code{$} operator or the
\code{attr} function.
}
\usage{
GetTaxa(object, ...)
GetLoci(object, ...)
GetLocDepth(object, ...)
GetContamRate(object, ...)
SetContamRate(object, value, ...)
nTaxa(object, ...)
nAlleles(object, ...)
GetAlleleNames(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{value}{
A value to assign.  For \code{SetContamRate}, a number generally ranging 
from zero to 0.01 indicating the expected rate of sample cross-contamination.
}
  \item{\dots}{
Additional arguments (none currently supported).
}
}
\value{
For \code{GetTaxa} and \code{GetLoci}, a character vector listing taxa names
or loci names, respectively.  For \code{GetLocDepth}, a named matrix with
taxa in rows and loci in columns, giving the total read depth for each taxon
and locus.  For \code{GetContamRate}, a number indicating the expected 
contamination rate that is stored in the object.  For \code{SetContamRate}, a
\code{"RADdata"} object with an updated contamination rate.
For \code{nTaxa}, the number of taxa in the object.  For \code{nAlleles}, the
number of alleles across all loci in the object.  For \code{GetAlleleNames},
the names of all alleles.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{SetBlankTaxa}} for functions that assign taxa to particular roles.
}
\examples{
data(exampleRAD)
GetTaxa(exampleRAD)
GetLoci(exampleRAD)
GetLocDepth(exampleRAD)
GetContamRate(exampleRAD)
exampleRAD <- SetContamRate(exampleRAD, 0.0000001)
GetContamRate(exampleRAD)
nTaxa(exampleRAD)
nAlleles(exampleRAD)
GetAlleleNames(exampleRAD)
}

\keyword{ methods }
\keyword{ utilities }
