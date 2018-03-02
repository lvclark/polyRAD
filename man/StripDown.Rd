\name{StripDown}
\alias{StripDown}
\alias{StripDown.RADdata}
\title{
Remove Unneeded Slots to Conserve Memory
}
\description{
This function is designed to be used after a \code{\link{RADdata}} object
has been processed by one of the \link[dest=IterateHWE]{pipeline} functions.  
Slots that are no longer needed are removed in order to conserve memory.
}
\usage{
StripDown(object, ...)
\method{StripDown}{RADdata}(object, 
          remove.slots = c("depthSamplingPermutations",
                           "depthRatio", "antiAlleleDepth",
                           "genotypeLikelihood", "priorProb"),
          \dots)
}
\arguments{
  \item{object}{
A \code{RADdata} object.
}
\item{remove.slots}{
A character vector listing slots that will be removed.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
The default slots that are removed take up a lot of memory but are not used
by the export functions.  Other slots to consider removing are 
\code{alleleFreq}, \code{alleleFreqByTaxa}, \code{PCA}, \code{locDepth}, and 
\code{alleleDepth}.  Of course, if you
have custom uses for some of the slots that are removed by default, you can 
change the \code{remove.slots} vector to not include them.

The function will throw an error if the user attempts to remove key slots that
are needed for export and downstream analysis, including:

\itemize{
\item alleles2loc
\item alleleNucleotides
\item locTable
\item priorProbPloidies
\item possiblePloidies
\item ploidyChiSq
\item posteriorProb
}

}
\value{
A \code{RADdata} object 
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{SubsetByTaxon}}, \code{\link{SubsetByLocus}}
}
\examples{
# load a dataset for this example
data(exampleRAD)

# run a pipeline
exampleRAD <- IterateHWE(exampleRAD)

# check the size of the resulting object
object.size(exampleRAD)

# remove unneeded slots
exampleRAD <- StripDown(exampleRAD)

# check object size again
object.size(exampleRAD)
}

\keyword{ manip }
