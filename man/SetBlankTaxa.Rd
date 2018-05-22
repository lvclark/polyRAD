\name{SetBlankTaxa}
\alias{SetBlankTaxa}
\alias{GetBlankTaxa}
\alias{SetDonorParent}
\alias{GetDonorParent}
\alias{SetRecurrentParent}
\alias{GetRecurrentParent}
\title{
Functions to Assign Taxa to Specific Roles
}
\description{
These functions are used for assigning and retrieving taxa from a 
\code{"\link{RADdata}"} object that serve particular roles in the dataset.
Blank taxa can be used for estimating the contamination rate (see 
\code{\link{EstimateContaminationRate}}), and the donor and recurrent parents
are used for determining
expected genotype distributions in mapping populations.  Many functions
in \pkg{polyRAD} will automatically exclude taxa from analysis if they
have been assigned to one of these roles.
}
\usage{
SetBlankTaxa(object, value)
GetBlankTaxa(object, ...)
SetDonorParent(object, value)
GetDonorParent(object, ...)
SetRecurrentParent(object, value)
GetRecurrentParent(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
  \item{value}{
A character string (or a character vector for \code{SetBlankTaxa}) indicating
the taxon or taxa to be assigned to the role.
}
\item{\dots}{
Other arguments (none currently supported).
}
}
\value{
For the \dQuote{Get} functions, a character vector indicating the taxon or taxa
that have been assigned to that role.  For the \dQuote{Set} functions, a
\code{"\link{RADdata}"} object identical to the one passed to the function, but
with new taxa assigned to that role.
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypePriorProb_Mapping2Parents}}
}
\examples{
# assign parents in a mapping population
data(exampleRAD_mapping)
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")
GetDonorParent(exampleRAD_mapping)
GetRecurrentParent(exampleRAD_mapping)

# assign blanks
exampleRAD_mapping <- SetBlankTaxa(exampleRAD_mapping, 
                                   c("progeny019", "progeny035"))
GetBlankTaxa(exampleRAD_mapping)
}

\keyword{ methods }
\keyword{ utilities }
