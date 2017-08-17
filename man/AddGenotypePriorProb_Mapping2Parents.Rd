\name{AddGenotypePriorProb_Mapping2Parents}
\alias{AddGenotypePriorProb_Mapping2Parents}
\alias{AddGenotypePriorProb_Mapping2Parents.RADdata}
\title{
Expected Genotype Frequencies in Mapping Populations
}
\description{
Given the most likely genotypes of two parent taxa, expected genotype
frequencies are estimated for a population of progeny and are added to the
\code{"\link{RADdata}"} object in the \code{$priorProb} slot.
}
\usage{
AddGenotypePriorProb_Mapping2Parents(object, ...)
AddGenotypePriorProb_Mapping2Parents.RADdata(object,
    donorParent = GetDonorParent(object), 
    recurrentParent = GetRecurrentParent(object), n.gen.backcrossing = 0,
    n.gen.selfing = 0, donorParentPloidies = object$possiblePloidies,
    recurrentParentPloidies = object$possiblePloidies)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  Ideally this should be set up as a mapping
population using \code{\link{SetDonorParent}}, \code{\link{SetRecurrentParent}},
and \code{\link{AddAlleleFreqMapping}}.
}
  \item{\dots}{
Additional arguments, listed below, to be passed to the method for 
\code{"RADdata"} objects.
}
\item{donorParent}{A character string indicating which taxon is the donor 
parent.  If backcrossing was not performed, it does not matter which was the
donor or recurrent parent.}
\item{recurrentParent}{A character string indicating which taxon is the
recurrent parent.}
\item{n.gen.backcrossing}{
An integer, zero or greater, indicating how many generations of backcrossing
to the recurrent parent were performed.
}
\item{n.gen.selfing}{
An integer, zero or greater, indicating how many generations of selfing were
performed.
}
\item{donorParentPloidies}{
A list, where each item in the list is an integer vector indicating a 
potential inheritance mode that could be observed among loci in the donor
parent.  \code{2} indicates diploid, \code{4} indicates autotetraploid, 
\code{c(2, 2)} indicates allotetraploid, \emph{etc.}
}
\item{recurrentParentPloidies}{
A list in the same format as \code{donorParentPloidies} indicating 
inheritance modes that could be observed among loci in the recurrent
parent.
}
}
\details{
The function first determines which combinations of inheritance modes from
the two parents should be examined in the progeny.  The expected progeny
ploidy must be in \code{object$possiblePloidies} for a given combination
to be examined.

The most likely genotypes for the two parents are estimated using
\code{\link{GetLikelyGen}}.  For each ploidy being examined, F1 genotype
probabilities are then calculated.

If both selfing and backcrossing were performed, selfing is assumed to have 
been done after backcrossing.  Genotype probabilities are updated for each
backcrossing generation, then each selfing generation.
}
\value{
A \code{"RADdata"} object identical to that passed to the function, but with
data stored in two new slots:

\item{priorProb }{A list of matrices, with one matrix per possible ploidy
of offspring.  For each matrix, allele copy number (from zero to the total
ploidy) is in rows, and alleles are in columns.  Each value is the probability
of sampling an individual with that allele copy number from the population.}
\item{priorProbPloidies }{A list in the same format as 
\code{object$possiblePloidies}, and the same length as 
\code{object$priorProb}.  Each item in the list is a vector indicating the
inheritance mode for the corresponding matrix in \code{object$priorProb}.}
}
\author{
Lindsay V. Clark
}
\note{
For the time being, in allopolyploids it is assumed that copies of an allele
are distributed among as few isoloci as possible.  For example, if an
autotetraploid genotype had two copies of allele A and two copies of allele
B, it is assumed to be AA BB rather than AB AB.  This may be remedied in the
future by examining distribution of genotype likelihoods.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    UseMethod("AddGenotypePriorProb_Mapping2Parents", object)
  }
}
\keyword{ methods }
\keyword{ distribution }
