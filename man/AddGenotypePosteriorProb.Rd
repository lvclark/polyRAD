\name{AddGenotypePosteriorProb}
\alias{AddGenotypePosteriorProb}
\title{
Estimate Posterior Probabilities of Genotypes
}
\description{
Given a \code{"\link{RADdata}"} object containing genotype prior probabilities
and genotype likelihoods, this function estimates genotype posterior
probabilities and adds them to the \code{$posteriorProb} slot of the object.
}
\usage{
AddGenotypePosteriorProb(object, ...)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  Prior genotype probabilities and
genotype likelihood should have already been added.
}
  \item{\dots}{
Potential future arguments (none currently in use).
}
}
%\details{
%If \code{\link{AddPriorTimesLikelihood}} has not already been run on the
%object, it will be run by \code{AddGenotypePosteriorProb} in order to
%perform the necessary calculations.
%}
\value{
A \code{"RADdata"} object identical to that passed to the function, but 
with a two-dimensional list added to the \code{$posteriorProb} slot.  Rows of
the list correspont to \code{object$possiblePloidies}, and columns to unique
values in \code{object$taxaPloidy}, similarly to \code{object$priorProb}. Each
item of the list is a three dimensional array, with allele copy number in the
first dimension, taxa in the second dimension, and alleles in the third
dimension. For each allele and taxa, posterior probabilities will sum to one
across all potential allele copy numbers.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypeLikelihood}}, 
\code{\link{AddGenotypePriorProb_Mapping2Parents}}
}
\examples{
# load dataset and set some parameters
data(exampleRAD_mapping)
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")
exampleRAD_mapping <- AddAlleleFreqMapping(exampleRAD_mapping,
                                           expectedFreqs = c(0.25, 0.75),
                                           allowedDeviation = 0.08)
exampleRAD_mapping <- AddGenotypeLikelihood(exampleRAD_mapping)
exampleRAD_mapping <- AddGenotypePriorProb_Mapping2Parents(exampleRAD_mapping,
                                                      n.gen.backcrossing = 1)
                                                      
# estimate posterior probabilities
exampleRAD_mapping <- AddGenotypePosteriorProb(exampleRAD_mapping)
# examine the results
exampleRAD_mapping$posteriorProb[[1,1]][,3,]
}

\keyword{ methods }
\keyword{ array }
