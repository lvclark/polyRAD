\name{AddAlleleLinkages}
\alias{AddAlleleLinkages}
\alias{AddAlleleLinkages.RADdata}
\alias{AddGenotypePriorProb_LD}
\title{
Identify and Utilize Linked Alleles for Estimating Genotype Priors
}
\description{
\code{AddAlleleLinkages} finds alleles, if any, in linkage disequilibrium 
with each allele in a \code{\link{RADdata}} object, and computes a correlation
coefficient representing the strength of the linkage.  
\code{AddGenotypePriorProb_LD} adds a second set of prior genotype 
probabilities to a \code{RADdata} object based on the genotype posterior
probabilities at linked alleles.
}
\usage{
AddAlleleLinkages(object, ...)
\method{AddAlleleLinkages}{RADdata}(object, type, linkageDist, minCorr,
                  excludeTaxa = character(0), \dots)
                  
AddGenotypePriorProb_LD(object, ...)
\method{AddGenotypePriorProb_LD}{RADdata}(object, type, \dots)
}
\arguments{
  \item{object}{
A \code{RADdata} object with genomic alignment data stored in 
\code{object$locTable$Chr} and \code{object$locTable$pos}.
}
\item{type}{
A character string, either \dQuote{mapping}, \dQuote{hwe}, or 
\dQuote{popstruct}, to indicate the type of population being analyzed.
}
\item{linkageDist}{
A number, indicating the distance in basepairs from a locus within which to
search for linked alleles.
}
\item{minCorr}{
A number ranging from zero to one indicating the minimum correlation needed
for an allele to be used for genotype prediction at another allele.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
These functions are primarily designed to be used internally by the
\link[IterateHWE]{pipeline} functions.

\code{AddAlleleLinkages} obtains genotypic values using 
\code{GetWeightedMeanGenotypes}, then regresses those values for a given 
allele against those values for nearby alleles to obtain correlation coefficients.
For the population structure model, the genotypic values for an allele are
first regressed on the PC axes from \code{object$PCA}, then the residuals are
regressed on the genotypic values at nearby alleles to obtain correlation
coefficients.

\code{AddGenotypePriorProb_LD} makes a second set of priors in addition to
\code{object$priorProb}.  This second set of priors has one value per
inheritance mode per taxon per allele per possible allele copy number.
Where \eqn{K} is the ploidy, with allele copy number \eqn{c} ranging from 0 to
\eqn{K}, \eqn{i} is an allele, \eqn{j} is a linked allele at a different locus
out of \eqn{J} total alleles linked to \eqn{i},
\eqn{r_{ij}} is the correlation coefficient between those alleles, \eqn{t} is a
taxon, \eqn{post_{cjt}} is the posterior probability of a given allele copy
number for a given allele in a given taxon, and \eqn{prior_{cit}} is the
prior probability for a given allele copy number for a given allele in a given
taxon based on linkage alone:

\deqn{prior_{cit} = \frac{\prod_{j = 1}^J{post_{cjt} * r_{ij} + (1 - r_{ij})/(K + 1)}}{\sum_{c = 0}0^K{\prod_{j = 1}^J{post_{cjt} * r_{ij} + (1 - r_{ij})/(K + 1)}}}}{
prior_{cit} = (\prod_{j = 1}^J{post_{cjt} * r_{ij} + (1 - r_{ij})/(K + 1)})/(\sum_{c = 0}0^K{\prod_{j = 1}^J{post_{cjt} * r_{ij} + (1 - r_{ij})/(K + 1)}})}

}

\value{
A \code{RADdata} object is returned.  For \code{AddAlleleLinkages}, it has a new slot
called \code{$alleleLinkages} that is a list, with one item in the list for each
allele in the dataset.  Each item is a data frame, with indices for linked alleles
in the first column, and correlation coefficients in the second column.

For \code{AddGenotypePriorProb_LD}, the object has a new slot called 
\code{$priorProbLD}.  This is a list much like \code{$posteriorProb}, with one list
item per inheritance mode, and each item being an array with allele copy number in
the first dimension, taxa in the second dimension, and alleles in the third dimension.
Values indicate genotype prior probabilities based on linked alleles alone.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypePriorProb}}
}
\examples{
# load example dataset
mydata <- VCF2RADdata(system.file("extdata", "Msi01genes.vcf",
                                  package = "polyRAD"))
# Run non-LD pop structure pipeline
mydata <- IteratePopStruct(mydata, tol = 0.01, nPcsInit = 10)

# Add linkages
mydata <- AddAlleleLinkages(mydata, "popstruct", 1e4, 0.05)
# Get new posterior probabilities based on those linkages
mydata <- AddGenotypePriorProb_LD(mydata, "popstruct")

# Preview results
mydata$priorProbLD[[1]][,1:10,1:10]
}
\keyword{ array }
\keyword{ regression }
