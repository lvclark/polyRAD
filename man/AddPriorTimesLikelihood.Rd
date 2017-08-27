\name{AddPriorTimesLikelihood}
\alias{AddPriorTimesLikelihood}
\title{
Multiply Genotype Prior Probabilities by Genotype Likelihoods
}
\description{
In preparation for calculations to be done by other functions, the genotype
likelihoods for each taxon, allele, and ploidy are multiplied by the
corresponding genotype prior probabilities, then stored back in the 
\code{"RADdata"} object in the \code{$priorTimesLikelihood} slot.
}
\usage{
AddPriorTimesLikelihood(object, ...)
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
\value{
A \code{"RADdata"} object identical to that passed to the function, but with
the \code{$priorTimesLikelihood} slot added.  This slot contains a list.  Each
element in the list corresponds to an element in \code{object$priorProbs}.
Each element is an array, with allele copy number in the first dimension, 
taxa in the second dimension, and alleles in the third dimension.
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
                                                      
# perform the multiplication
exampleRAD_mapping <- AddPriorTimesLikelihood(exampleRAD_mapping)
# examine the results
exampleRAD_mapping$priorTimesLikelihood[[1]][,50,] # for one progeny
exampleRAD_mapping$priorTimesLikelihood[[1]][,1,]
 # --> for the donor parent; not a good idea to use since the priors 
 #    aren't appropriate
}

\keyword{ arith }
\keyword{ array }
\keyword{ methods }
