\name{SubsetByTaxon}
\alias{SubsetByTaxon}
\alias{SubsetByTaxon.RADdata}
\title{
Create RADdata Object with a Subset of Taxa
}
\description{
This function is used for removing some of the taxa from a dataset stored in a
\code{\link{RADdata}} object.
}
\usage{
SubsetByTaxon(object, ...)
\method{SubsetByTaxon}{RADdata}(object, taxa, \dots)
}

\arguments{
  \item{object}{
A \code{RADdata} object.
}
  \item{taxa}{
A character or numeric vector indicating which taxa to retain in the output.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
This function may be used for subsetting a \code{RADdata} object either
immediately after data import, or after additional analysis has been
performed.  Note however that estimation of allele frequencies, genotype
prior probabilities, PCA, \emph{etc.} are very dependent on what samples
are included in the dataset.  If those calculations have already been 
performed, the results will be transferred to the new object but not
recalculated.
}
\value{
A \code{RADdata} object containing only the taxa listed in \code{taxa}.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{SubsetByLocus}}
}
\examples{
# load data for this example
data(exampleRAD)
exampleRAD

# just keep the first fifty taxa
subsetRAD <- SubsetByTaxon(exampleRAD, 1:50)
subsetRAD

}

\keyword{ manip }
