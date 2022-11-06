\name{exampleRAD}
\alias{exampleRAD}
\alias{exampleRAD_mapping}
\alias{Msi01genes}
\docType{data}
\title{
Miniature Datasets for Testing polyRAD Functions
}
\description{
\code{exampleRAD} and \code{exampleRAD_mapping} are two very small
simulated \code{"\link{RADdata}"} datasets for testing polyRAD
functions.  Each has four loci.  \code{exampleRAD} is a 
natural population of 100 individuals with a mix of diploid and tetraploid 
loci, with 80 individuals diploid and 20 individuals triploid.
\code{exampleRAD_mapping} is a diploid BC1 mapping population with two parents
and 100 progeny.
\code{Msi01genes} is a \code{"RADdata"} object with 585 taxa and 24 loci,
containing real data from \emph{Miscanthus sinensis}, obtained by using
\code{\link{VCF2RADdata}} on the file Msi01genes.vcf.  Most individuals
in \code{Msi01genes} are diploid, with three haploids and one triploid.
}
\usage{
data(exampleRAD)
data(exampleRAD_mapping)
data(Msi01genes)
}
\format{
  See the format described in \code{"\link{RADdata}"}.
}

\source{
Randomly generated using a script available in
polyRAD/extdata/simulate_rad_data.R.

\emph{M. sinensis} sequencing data available at
\url{https://www.ncbi.nlm.nih.gov//bioproject/PRJNA207721}, with full genotype
calls at \doi{10.13012/B2IDB-1402948_V1}.
}

\examples{
data(exampleRAD)
exampleRAD
data(exampleRAD_mapping)
exampleRAD_mapping
data(Msi01genes)
Msi01genes
}
\keyword{datasets}
