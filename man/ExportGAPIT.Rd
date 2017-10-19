\name{ExportGAPIT}
\alias{ExportGAPIT}
\alias{Export_rrBLUP_Amat}
\alias{Export_rrBLUP_GWAS}
\title{
Export RADdata Object for Use by Other R Packages
}
\description{
After a \code{"\link{RADdata}"} object has been run through a pipeline such as
\code{\link{IteratePopStruct}}, these functions can be used to export
the genotypes (weighted across all possible genotypes by
\code{\link{GetWeightedMeanGenotypes}}) to other R packages that can
perform genome-wide association and genomic prediction using genotypes
as a continuous numeric variable.
}
\usage{
ExportGAPIT(object, onePloidyPerAllele = FALSE)
Export_rrBLUP_Amat(object, naIfZeroReads = FALSE, 
                   onePloidyPerAllele = FALSE)
Export_rrBLUP_GWAS(object, naIfZeroReads = FALSE, 
                   onePloidyPerAllele = FALSE)
}
\arguments{
  \item{object}{
A \code{"RADdata"} object with posterior genotype probabilities
already estimated.
}
  \item{onePloidyPerAllele}{
Logical.  If \code{TRUE}, for each allele the inheritance mode with the
lowest \eqn{\chi ^ 2}{chi-squared} value is selected and is assumed to be
the true inheritance mode.  If \code{FALSE}, inheritance modes are weighted
by inverse \eqn{\chi ^ 2}{chi-squared} values for each allele, and mean
genotypes that have been weighted across inheritance modes are returned.
}
\item{naIfZeroReads}{
A logical indicating whether \code{NA} should be inserted into the 
output matrix for any taxa and loci where the total read depth for
the locus is zero.  If \code{FALSE}, the output for these genotypes is
essentially the weighted mean across prior genotype probabilities, since
prior and posterior genotype probabilities are equal when there are no
reads.
}
}
\details{
Both \pkg{GAPIT} and \pkg{rrBLUP} allow genotypes to be a continuous numeric variable.
\pkg{GAPIT} does not allow missing data, hence there is no \code{naIfZeroReads}
argument for \code{ExportGAPIT}.  Genotypes are exported on a scale of -1
to 1 for \pkg{rrBLUP}, and on a scale of 0 to 2 for \pkg{GAPIT}.

For all functions, one allele per marker is dropped.  For \code{ExportGAPIT}
and \code{Export_rrBLUP_GWAS}, chromosome and position are filled with dummy
data if they do not exist in \code{object$locTable}.

If the chromosomes in \code{object$locTable} are in character format,
\code{ExportGAPIT} will attempt to extract chromosome numbers.
}
\value{
For \code{ExportGAPIT}, a list:

\item{GD }{A data frame with taxa in the first column and alleles (markers)
in subsequent columns, containing the genotypes.  To be passed to the \code{GD}
argument for \code{GAPIT}.}
\item{GM }{A data frame with the name, chromosome number, and position of
every allele (marker).  To be passed to the \code{GM} argument for \code{GAPIT}.}

For \code{Export_rrBLUP_Amat}, a matrix with taxa in rows and alleles (markers)
in columns, containing genotype data.  This can be passed to \code{A.mat} in
\pkg{rrBLUP}.

For \code{Export_rrBLUP_Amat}, a data frame with alleles (markers) in rows.
The first three columns contain the marker names, chromosomes, and positions,
and the remaining columns each represent one taxon and contain the genotype
data.  This can be passed to the \code{GWAS} function in \pkg{rrBLUP}
}
\references{
\url{http://zzlab.net/GAPIT/}

Lipka, A. E., Tian, F., Wang, Q., Peiffer, J., Li, M., Bradbury, P. J., 
Gore, M. A., Buckler, E. S. and Zhang, Z. (2012) GAPIT: genome association 
and prediction integrated tool.  \emph{Bioinformatics} \bold{28}, 2397--2399.

Endelman, J.B. (2011) Ridge Regression and Other Kernels for Genomic 
Selection with R Package rrBLUP.  \emph{The Plant Genome} \bold{4}, 250--255.
}
\author{
Lindsay V. Clark
}
\note{
\pkg{rrBLUP} is available through CRAN, and \pkg{GAPIT} must be downloaded from
the Zhang lab website.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
