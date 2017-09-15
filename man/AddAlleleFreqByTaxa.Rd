\name{AddAlleleFreqByTaxa}
\alias{AddAlleleFreqByTaxa}
\alias{AddAlleleFreqByTaxa.RADdata}
\title{
Estimate Local Allele Frequencies for Each Taxon Based on Population Structure
}
\description{
This function estimates allele frequencies per taxon, rather than for the whole
population.  The best estimated genotypes (either \code{object$depthRatio} or
\code{GetWeightedMeanGenotypes(object)}) are regressed against principal 
coordinate axes.  The regression coefficients are then in turn used to 
predict allele frequencies from PC axes.  Allele frequencies outside of a 
user-defined range are then adjusted so that they fall within that range.
}
\usage{
AddAlleleFreqByTaxa(object, ...)
\method{AddAlleleFreqByTaxa}{RADdata}(object, minfreq = 0.0001, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  \code{\link{AddPCA}} should have already
been run.
}
\item{minfreq}{
The minimum allowable allele frequency to be output.  The maximum allowable
allele frequency will be calculated as \code{1 - minfreq}.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
For every allele, all PC axes stored in \code{object$PCA} are used for 
generating regression coefficients and making predictions, regardless of whether
they are significantly associated with the allele.

\code{object$depthRatio} has missing data for loci with no reads; these missing
data are omitted on a per-allele basis when calculating regression coefficients.
However, allele frequencies are output for all taxa at all alleles, because
there are no missing data in the PC axes.  The output of 
\code{\link{GetWeightedMeanGenotypes}} has no missing data, so missing data are
not an issue when calculating regression coefficients using that method.

After predicting allele frequencies from the regression coefficients, the
function loops through all loci and taxa to adjust allele frequencies if necessary.
This is needed because otherwise some allele frequencies will be below zero or
above one (typically in subpopulations where alleles are near fixation), 
which interferes with prior genotype probability estimation.  For a 
given taxon and locus, any allele frequencies below \code{minfreq} are adjusted
to be equal to \code{minfreq}, and any allele frequencies above \code{1 - minfreq}
are adjusted to be \code{1 - minfreq}.  Remaining allele frequencies are adjusted
so that all allele frequencies for the taxon and locus sum to one.
}
\value{
A \code{"RADdata"} object identical to the one passed to the function, but with
a matrix of allele frequencies added to the \code{$alleleFreqByTaxa} slot.  Taxa
are in rows and alleles in columns.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypePriorProb_ByTaxa}}
}
\examples{
# load data
data(exampleRAD)
# do PCA
exampleRAD <- AddPCA(exampleRAD, nPcsInit = 3)

# get allele frequencies
exampleRAD <- AddAlleleFreqByTaxa(exampleRAD)

exampleRAD$alleleFreqByTaxa[1:10,]
}

\keyword{ regression }
\keyword{ methods }
