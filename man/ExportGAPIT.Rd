\name{ExportGAPIT}
\alias{ExportGAPIT}
\alias{Export_rrBLUP_Amat}
\alias{Export_rrBLUP_GWAS}
\alias{Export_TASSEL_Numeric}
\alias{Export_polymapR}
\alias{Export_polymapR_probs}
\alias{Export_MAPpoly}
\alias{Export_GWASpoly}
\alias{Export_Structure}
\alias{Export_adegenet_genind}
\title{
Export RADdata Object for Use by Other R Packages
}
\description{
After a \code{"\link{RADdata}"} object has been run through a pipeline such as
\code{\link{IteratePopStruct}}, these functions can be used to export
the genotypes to R packages and other software that can
perform genome-wide association and genomic prediction.  \code{ExportGAPIT},
\code{Export_rrBLUP_Amat}, \code{Export_rrBLUP_GWAS}, \code{Export_GWASpoly},
and \code{Export_TASSEL_Numeric} all export continuous numerical genotypes 
generated by \code{\link{GetWeightedMeanGenotypes}}.  \code{Export_polymapR},
\code{Export_Structure}, and \code{Export_adegenet_genind} use 
\code{\link{GetProbableGenotypes}} to export discrete 
genotypes.  \code{Export_MAPpoly} and \code{Export_polymapR_probs} export
genotype posterior probabilities.
}
\usage{
ExportGAPIT(object, onePloidyPerAllele = FALSE)

Export_rrBLUP_Amat(object, naIfZeroReads = FALSE, 
                   onePloidyPerAllele = FALSE)

Export_rrBLUP_GWAS(object, naIfZeroReads = FALSE, 
                   onePloidyPerAllele = FALSE)

Export_TASSEL_Numeric(object, file, naIfZeroReads = FALSE,
                      onePloidyPerAllele = FALSE)

Export_polymapR(object, naIfZeroReads = TRUE,
                progeny = GetTaxa(object)[!GetTaxa(object) \%in\% 
                  c(GetDonorParent(object), GetRecurrentParent(object),
                    GetBlankTaxa(object))])

Export_polymapR_probs(object, maxPcutoff = 0.9,
                      correctParentalGenos = TRUE,
                      multiallelic = "correct")

Export_MAPpoly(object, file, pheno = NULL, ploidyIndex = 1,
               progeny = GetTaxa(object)[!GetTaxa(object) \%in\% 
                 c(GetDonorParent(object), GetRecurrentParent(object),
                   GetBlankTaxa(object))],
               digits = 3)

Export_GWASpoly(object, file, naIfZeroReads = TRUE, postmean = TRUE, digits = 3,
                splitByPloidy = TRUE)

Export_Structure(object, file, includeDistances = FALSE, extraCols = NULL,
                 missingIfZeroReads = TRUE)

Export_adegenet_genind(object, ploidyIndex = 1)
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
essentially the mode (for \code{Export_polymapR} and \code{Export_GWASpoly})
or mean (for others) across prior genotype probabilities, since
prior and posterior genotype probabilities are equal when there are no
reads.
}
\item{file}{
A character string indicating a file path to which to write.
}
\item{pheno}{
A data frame or matrix of phenotypic values, with progeny in rows and traits
in columns.  Columns should be named.
}
\item{ploidyIndex}{
Index, within \code{object$possiblePloidies}, of the ploidy to be exported.
}
\item{progeny}{
A character vector indicating which individuals to export as progeny of the
cross.
}
\item{maxPcutoff}{
A cutoff for posterior probabilities, below which genotypes will be reported
as `NA` in the `geno` column.
}
\item{correctParentalGenos}{
Passed to \code{\link{GetProbableGenotypes}}. If \code{TRUE}, parental genotypes
are corrected based on progeny allele frequencies.
}
\item{multiallelic}{
Passed to \code{\link{GetProbableGenotypes}}. Under the default, genotypes at
multiallelic loci will be corrected to sum to the ploidy.
}
\item{digits}{
Number of decimal places to which to round genotype probabilities or posterior
mean genotypes in the output file.
}
\item{postmean}{
Logical. If \code{TRUE}, posterior mean genotypes will be output.  If
\code{FALSE}, discrete genotypes will be output.
}
\item{splitByPloidy}{
Logical. If \code{TRUE} and there are multiple \code{taxaPloidy} values in the
dataset, multiple files are written, one per ploidy.
}
\item{includeDistances}{
Logical. If \code{TRUE}, the second row of the Structure file will contain
distances between markers, which can be used by the linkage model in Structure.
}
\item{extraCols}{
An optional data frame, with one row per taxon, containing columns of data
to output to the left of the genotypes in the Structure file.
}
\item{missingIfZeroReads}{
See \code{naIfZeroReads}.  If \code{TRUE}, a value of \code{-9} will be output
for any genotypes with zero reads, indicating that those genotypes are missing.
}
}
\details{
\pkg{GAPIT}, \pkg{FarmCPU}, \pkg{rrBLUP}, \pkg{TASSEL}, and \pkg{GWASpoly} allow
genotypes to be a continuous numeric variable.  \pkg{MAPpoly} and \pkg{polymapR}
allow for import of genotype probabilities.
\pkg{GAPIT} does not allow missing data, hence there is no \code{naIfZeroReads}
argument for \code{ExportGAPIT}.  Genotypes are exported on a scale of -1
to 1 for \pkg{rrBLUP}, on a scale of 0 to 2 for \pkg{GAPIT} and \pkg{FarmCPU},
and on a scale of 0 to 1 for \pkg{TASSEL}.

For all functions except \code{Export_Structure} and \code{Export_adegenet_genind},
one allele per marker is dropped.  \code{Export_MAPpoly}
also drops alleles where one or both parental genotypes could not be determined,
and where both parents are homozygotes.

For \code{ExportGAPIT}
and \code{Export_rrBLUP_GWAS}, chromosome and position are filled with dummy
data if they do not exist in \code{object$locTable}.  For \code{Export_TASSEL_Numeric},
allele names are exported, but no chromosome or position information per se.

If the chromosomes in \code{object$locTable} are in character format,
\code{ExportGAPIT}, \code{Export_MAPpoly}, and \code{Export_GWASpoly} will
attempt to extract chromosome numbers.

For \pkg{polymapR} there must only be one possible inheritance mode across loci
(one value in \code{object$possiblePloidies}) in the \code{RADdata} object, although
triploid F1 populations derived from diploid and tetraploid parents are allowed.
See \code{\link{SubsetByPloidy}} for help reducing a \code{RADdata} object to a
single inheritance mode.

\pkg{MAPpoly} only
allows one ploidy, but \code{Export_MAPpoly} allows the user to select which
inheritance mode from the \code{RADdata} object to use.  (This is due to how internal
\pkg{polyRAD} functions are coded.)
}
\value{
For \code{ExportGAPIT}, a list:

\item{GD }{A data frame with taxa in the first column and alleles (markers)
in subsequent columns, containing the genotypes.  To be passed to the \code{GD}
argument for \code{GAPIT} or \code{FarmCPU}.}
\item{GM }{A data frame with the name, chromosome number, and position of
every allele (marker).  To be passed to the \code{GM} argument for \code{GAPIT}
or \code{FarmCPU}.}

For \code{Export_rrBLUP_Amat}, a matrix with taxa in rows and alleles (markers)
in columns, containing genotype data.  This can be passed to \code{A.mat} in
\pkg{rrBLUP}.

For \code{Export_rrBLUP_GWAS}, a data frame with alleles (markers) in rows.
The first three columns contain the marker names, chromosomes, and positions,
and the remaining columns each represent one taxon and contain the genotype
data.  This can be passed to the \code{GWAS} function in \pkg{rrBLUP}.

\code{Export_TASSEL_Numeric} and \code{Export_MAPpoly} write a file but does
not return an object.

For \code{Export_polymapR}, a matrix of integers indicating the most probable 
allele copy number, with markers in rows and individuals in columns.  The 
parents are listed first, followed by all progeny.

For \code{Export_polymapR_probs}, a data frame suitable to pass to the
\code{probgeno_df} argument of \link[polymapR]{checkF1}.  Note that under
default parameters, in some cases the \code{maxP}, \code{maxgeno}, and
\code{geno} columns may not actually reflect the maximum posterior probability
if genotype correction was performed.

For \code{Export_adegenet_genind}, a \code{\link[adegenet:genind]{"genind"}} object.

\code{Export_MAPpoly}, \code{Export_GWASpoly}, and \code{Export_Structure}
write files but do not return
an object.  Files output by \code{Export_GWASpoly} are comma delimited and
in numeric format.  Sample and locus names are included in the file output
by \code{Export_Structure}, and the number of rows for each sample is
equal to the highest ploidy as determined by the \code{taxaPloidy} slot and the
output of \code{GetProbableGenotypes}.
}
\references{
\bold{GAPIT and FarmCPU:}

\url{https://zzlab.net/GAPIT/}

Lipka, A. E., Tian, F., Wang, Q., Peiffer, J., Li, M., Bradbury, P. J., 
Gore, M. A., Buckler, E. S. and Zhang, Z. (2012) GAPIT: genome association 
and prediction integrated tool.  \emph{Bioinformatics} \bold{28}, 2397--2399.

\url{https://zzlab.net/FarmCPU/}

Liu, X., Huang, M., Fan, B., Buckler, E. S., Zhang, Z. (2016) Iterative usage
of fixed and random effects models for powerful and efficient genome-wide
association studies.  \emph{PLoS Genetics} \bold{12}, e1005767.

\bold{rrBLUP:}

Endelman, J.B. (2011) Ridge Regression and Other Kernels for Genomic 
Selection with R Package rrBLUP.  \emph{The Plant Genome} \bold{4}, 250--255.

\bold{TASSEL:}

\url{https://www.maizegenetics.net/tassel}

Bradbury, P. J., Zhang, Z., Kroon, D. E., Casstevens, T. M., Ramdoss, Y. 
and Buckler, E. S. (2007) TASSEL: Software for association mapping of complex 
traits in diverse samples. \emph{Bioinformatics} \bold{23}, 2633--2635.

\bold{polymapR:}

Bourke, P., van Geest, G., Voorrips, R. E., Jansen, J., Kranenberg, T.,
Shahin, A., Visser, R. G. F., Arens, P., Smulders, M. J. M. and Maliepaard, C.
(2018) polymapR: linkage analysis and genetic map construction from F1 
populations of outcrossing polyploids. \emph{Bioinformatics} \bold{34}, 
3496--3502.

\bold{MAPpoly:}

\url{https://github.com/mmollina/MAPpoly}

Mollinari, M. and Garcia, A. A. F. (2018) Linkage analysis and haplotype
phasing in experimental autopolyploid populations with high ploidy level using
hidden Markov models.  \emph{bioRxiv} doi: https://doi.org/10.1101/415232.

\bold{GWASpoly:}

\url{https://github.com/jendelman/GWASpoly}

Rosyara, U. R., De Jong, W. S., Douches, D. S., and Endelman, J. B.  (2016)
Software for Genome-Wide Association Studies in Autopolyploids and Its 
Application to Potato.  \emph{Plant Genome} \bold{9}.

\bold{Structure:}

\url{https://web.stanford.edu/group/pritchardlab/structure.html}

Hubisz, M. J., Falush, D., Stephens, M. and Pritchard, J. K. (2009)
Inferring weak population structure with the assistance of sample
group information.  \emph{Molecular Ecology Resources} \bold{9},
1322--1332.

Falush, D., Stephens, M. and Pritchard, J. K.  (2007)
Inferences of population structure using multilocus genotype data:
dominant markers and null alleles.  \emph{Molecular Ecology Notes}
\bold{7}, 574--578

Falush, D., Stephens, M. and Pritchard, J. K.  (2003)
Inferences of population structure using multilocus genotype data:
linked loci and correlated allele frequencies.
\emph{Genetics} \bold{164}, 1567--1587.

Pritchard, J. K., Stephens, M. and Donnelly, P.  (2000)
Inference of population structure using multilocus genotype data.
\emph{Genetics} \bold{155}, 945--959.
}
\author{
Lindsay V. Clark
}
\note{
\pkg{rrBLUP} and \pkg{polymapR} are available through CRAN, and \pkg{GAPIT}
and \pkg{FarmCPU} must be downloaded from the Zhang lab website.  \pkg{MAPpoly}
is available on GitHub but not yet on CRAN.  \pkg{GWASpoly} is available from
GitHub.

In my experience with \pkg{TASSEL 5}, numerical genotype files that are too large do 
not load/display properly.  If you run into this problem I recommend using
\code{\link{SplitByChromosome}} to split your \code{RADdata} object into 
multiple smaller objects, which can then be exported to separate files using
\code{Export_TASSEL_Numeric}.  If performing GWAS, you may also need to compute 
a kinship matrix using separate software such as \pkg{rrBLUP}.
}

\seealso{
\code{\link{GetWeightedMeanGenotypes}}, \code{\link{RADdata2VCF}}
}
\examples{
# load example dataset
data(exampleRAD)
# get genotype posterior probabilities
exampleRAD <- IterateHWE(exampleRAD)

# export to GAPIT
exampleGAPIT <- ExportGAPIT(exampleRAD)

# export to rrBLUP
example_rrBLUP_A <- Export_rrBLUP_Amat(exampleRAD)
example_rrBLUP_GWAS <- Export_rrBLUP_GWAS(exampleRAD)

# export to TASSEL
outfile <- tempfile() # temporary file for example
Export_TASSEL_Numeric(exampleRAD, outfile)

# for mapping populations
data(exampleRAD_mapping)

# specify donor and recurrent parents
exampleRAD_mapping <- SetDonorParent(exampleRAD_mapping, "parent1")
exampleRAD_mapping <- SetRecurrentParent(exampleRAD_mapping, "parent2")

# run the pipeline
exampleRAD_mapping <- PipelineMapping2Parents(exampleRAD_mapping)

# convert to polymapR format
examplePMR <- Export_polymapR(exampleRAD_mapping)

examplePMR2 <- Export_polymapR_probs(exampleRAD_mapping)

# export to MAPpoly
outfile2 <- tempfile() # temporary file for example
# generate a dummy phenotype matrix containing random numbers
mypheno <- matrix(rnorm(200), nrow = 100, ncol = 2,
                  dimnames = list(GetTaxa(exampleRAD_mapping)[-(1:2)],
                                  c("Height", "Yield")))
Export_MAPpoly(exampleRAD_mapping, file = outfile2, pheno = mypheno)

# load data into MAPpoly
# require(mappoly)
# mydata <- read_geno_prob(outfile2)

# export to GWASpoly
outfile3 <- tempfile() # temporary file for example
Export_GWASpoly(SubsetByPloidy(exampleRAD, list(2)), outfile3)

# export to Structure
outfile4 <- tempfile() # temporary file for example
Export_Structure(exampleRAD, outfile4)

# export to adegenet
if(requireNamespace("adegenet", quietly = TRUE)){
  mygenind <- Export_adegenet_genind(exampleRAD)
}
}

\keyword{ manip }
\keyword{ file }
