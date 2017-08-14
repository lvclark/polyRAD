\name{exampleRAD}
\alias{exampleRAD}
\alias{exampleRAD_mapping}
\docType{data}
\title{
Miniature Datasets for Testing polyRAD Functions
}
\description{
\code{exampleRAD} and \code{exampleRAD_mapping} are two very small
simulated \code{"\link{RADdata}"} datasets for testing polyRAD
functions.  Each has four loci.  \code{exampleRAD} is a 
natural population of 100 individuals with a mix of diploid and tetraploid 
loci.
\code{exampleRAD_mapping} is a diploid BC1 mapping population with two parents
and 100 progeny.
}
\usage{
data(exampleRAD)
data(exampleRAD_mapping)
}
\format{
  See the format described in \code{"\link{RADdata}"}.
}

\source{
Randomly generated using a script available in
polyRAD/extdata/simulate_rad_data.R.
}

\examples{
data(exampleRAD)
exampleRAD
data(exampleRAD_mapping)
exampleRAD_mapping
}
\keyword{datasets}
