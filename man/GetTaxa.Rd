\name{Accessors}
\alias{GetTaxa}
\alias{GetLoci}
\alias{GetLocDepth}
\alias{GetContamRate}
\alias{SetContamRate}
\alias{nTaxa}
\alias{nLoci}
\alias{nAlleles}
\alias{GetAlleleNames}
\alias{GetTaxaPloidy}
\alias{SetTaxaPloidy}
\alias{GetTaxaByPloidy}
\alias{GetTaxaByPloidy.RADdata}
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
nLoci(object, ...)
nAlleles(object, ...)
GetAlleleNames(object, ...)
GetTaxaPloidy(object, ...)
SetTaxaPloidy(object, value, ...)
GetTaxaByPloidy(object, ...)

\method{GetTaxaByPloidy}{RADdata}(object, ploidy, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{value}{
A value to assign.  For \code{SetContamRate}, a number generally ranging 
from zero to 0.01 indicating the expected rate of sample cross-contamination.
For \code{SetTaxaPloidy}, a vector of integers indicating ploidy, with one
value for each taxon.  If the vector for \code{SetTaxaPloidy} is named, the
names should correspond to taxa names in the object.
}
\item{ploidy}{
An integer indicating a single ploidy for which to return taxa.
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
For \code{nTaxa}, the number of taxa in the object.  For \code{nLoci}, the 
number of loci in the object.  For \code{nAlleles}, the
number of alleles across all loci in the object.  For \code{GetAlleleNames},
the names of all alleles.  For \code{GetTaxaPloidy}, a named integer vector
indicating the ploidy of each taxon.  For \code{SetTaxaPloidy}, a
\code{"RADdata"} object with the taxa ploidies updated.  For
\code{GetTaxaByPloidy}, a character vector listing taxa.
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
GetTaxaPloidy(exampleRAD)
exampleRAD <- SetTaxaPloidy(exampleRAD, rep(c(2, 5), time = c(75, 25)))
GetTaxaByPloidy(exampleRAD, 2)
}

\keyword{ methods }
\keyword{ utilities }
