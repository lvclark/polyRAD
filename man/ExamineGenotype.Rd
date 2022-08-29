\name{ExamineGenotype}
\alias{ExamineGenotype}
\alias{ExamineGenotype.RADdata}
\title{
Plots to Examine Genotype Calling at a Single Taxon and Allele
}
\description{
For a given taxon and allele, this function generates barplots showing read
depth ratio, posterior mean genotype, genotype prior probabilities, genotype
likelihoods, and genotype posterior probabilities. It is intended as a sanity
check on genotype calling, as well as a means to visually demonstrate the
concept of Bayesian genotype calling.
}
\usage{
ExamineGenotype(object, ...)

\method{ExamineGenotype}{RADdata}(object, taxon, allele, pldindex = 1, \dots)
}
\arguments{
  \item{object}{
A \code{RADdata} object for which genotype calling has already been performed.
}
\item{taxon}{
A single character string indicating the taxon to show.
}
\item{allele}{
A single character string indicating the allele to show.
}
\item{pldindex}{
An index of which inheritance mode to use within \code{object$possiblePloidies}.
}
\item{...}{
Other arguments (none implemented).
}
}
\value{
A barplot is generated. Invisibly, a list is returned:
\item{alleleDepth }{Sequence read depth for the selected allele.}
\item{antiAlleleDepth }{Sequence read depth for all other alleles at the locus.}
\item{depthRatio }{Proportion of reads at this taxon and locus belonging to this allele.}
\item{priorProb }{A vector of genotype prior probabilities.}
\item{genotypeLikelhood }{A vector of genotype likelihoods.}
\item{posteriorProb }{A vector of genotype posterior probabilities.}
\item{postMean }{The posterior mean genotype on a scale of 0 to 1.}
}
\author{
Lindsay V. Clark
}
\examples{
data(exampleRAD)

exampleRAD <- IterateHWE(exampleRAD)

eg <- ExamineGenotype(exampleRAD, "sample088", "loc1_T")
}
