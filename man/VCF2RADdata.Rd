\name{VCF2RADdata}
\alias{VCF2RADdata}
\title{
Create a RADdata Object from a VCF File
}
\description{
This function reads a Variant Call Format (VCF) file containing allelic read depth
and SNP alignment positions, such as can be produced by TASSEL or GATK, and 
generates a \code{\link{RADdata}} dataset to be used for genotype calling in
\pkg{polyRAD}.
}
\usage{
VCF2RADdata(file, phaseSNPs = TRUE, tagsize = 80, refgenome = NULL, 
            tol = 0.01, al.depth.field = "AD", min.ind.with.reads = 200, 
            min.ind.with.minor.allele = 10, possiblePloidies = list(2), 
            taxaPloidy = 2L, contamRate = 0.001, 
            samples = VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(file)),
            svparam = VariantAnnotation::ScanVcfParam(fixed = "ALT", info = NA,
                                                      geno = al.depth.field,
                                                      samples = samples), 
            yieldSize = 5000, expectedAlleles = 5e+05, expectedLoci = 1e+05,
            maxLoci = NA)
}
\arguments{
  \item{file}{
The path to a VCF file to be read.  This can be uncompressed, bgzipped using
Samtools or Bioconductor, or a \code{TabixFile} object from Bioconductor.
}
  \item{phaseSNPs}{
If \code{TRUE}, markers that appear to have come from the same set of reads
will be phased and grouped into haplotypes.  Otherwise, each row of the file
will be kept as a distinct marker.
}
  \item{tagsize}{
The read length, minus any barcode sequence, that was used for genotyping.  In TASSEL,
this is the same as the kmerLength option.  This argument is used for grouping
SNPs into haplotypes and is ignored if \code{phaseSNPs = FALSE}.
}
  \item{refgenome}{
Optional.  The name of a FASTA file, or an \code{FaFile} object, containing
the reference genome.  When grouping SNPs into haplotypes, if provided this 
reference genome is used to insert non-variable nucleotides between the variable
nucleotides in the \code{alleleNucleotides} slot of the \code{RADdata} output.
Ignored if \code{phaseSNPs = FALSE}.  Useful if exact SNP positions need to be
retained for downstream analysis after genotype calling in \pkg{polyRAD}.
In particular this argument is necessary if you plan to export genotype calls
back to VCF.
}
  \item{tol}{
The proportion by which two SNPs can differ in read depth and still be merged
into one group for phasing.  Ignored if \code{phaseSNPs = FALSE}.
}
  \item{al.depth.field}{
The name of the genotype field in the VCF file that contains read depth at 
each allele.  This should be "AD" unless your format is very unusual.
}
  \item{min.ind.with.reads}{
Integer used for filtering SNPs.  To be retained, a SNP must have at least
this many samples with reads.
}
  \item{min.ind.with.minor.allele}{
Integer used for filtering SNPs.  To be retained, a SNP must have at least
this many samples with the minor allele.  When there are more than two
alleles, at least two alleles must have at least this many samples with
reads for the SNP to be retained.
}
  \item{possiblePloidies}{
A list indicating inheritance modes that might be encountered in the 
dataset.  See \code{\link{RADdata}}.
}
  \item{taxaPloidy}{
A single integer, or an integer vector with one value per taxon, indicating
ploidy.  See \code{\link{RADdata}}.
}
  \item{contamRate}{
A number indicating the expected sample cross-contamination rate.  See
\code{\link{RADdata}}.
}
  \item{samples}{
A character vector containing the names of samples from the file to 
export to the \code{RADdata} object.  The default is all samples.
If a subset is provided, filtering with \code{min.ind.with.reads} and
\code{min.ind.with.minor.allele} is performed within that subset.  Ignored
if a different \code{samples} argument is provided within \code{svparam}.
}
  \item{svparam}{
A \code{\link[VariantAnnotation:ScanVcfParam-class]{ScanVcfParam}} object to be
used with \code{\link[VariantAnnotation:readVcf-methods]{readVcf}}.  The primary
reasons to change this from the default would be 1) if you want additional
FIXED or INFO fields from the file to be exported to the \code{locTable}
slot of the \code{RADdata} object, and/or 2) if you only want to import
particular regions of the genome, as specified with the \code{which}
argument of \code{ScanVcfParam}.
}
  \item{yieldSize}{
An integer indicating the number of lines of the file to read at once.
Increasing this number will make the function faster but consume more RAM.
}
  \item{expectedAlleles}{
An integer indicating the approximate number of alleles that are expected
to be imported after filtering and phasing.  If this number is too low,
the function may slow down considerably.  Increasing this number 
increases the amount of RAM used by the function.
}
  \item{expectedLoci}{
An integer indicating the approximate number of loci that are expected
to be imported after filtering and phasing.  If this number is too low,
the function may slow down considerably.  Increasing this number 
increases the amount of RAM used by the function.
}
\item{maxLoci}{
An integer indicating the approximate maximum number of loci to return.  If
provided, the function will stop reading the file once it has found at least
this many loci that pass filtering and phasing.  This argument is intended to
be used for generating small \code{RADdata} objects for testing purposes, and
should be left \code{NA} under normal circumstances.
}
}
\details{
This function requires the BioConductor package \pkg{VariantAnnotation}.
See \url{https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}
for installation instructions.

If you anticipate running \code{VCF2RADdata} on the same file more than once,
it is recommended to run \code{bgzip} and \code{indexTabix} from the package
\pkg{Rsamtools} once before running \code{VCF2RADdata}.  See examples.

\code{min.ind.with.minor.allele} is used for filtering SNPs as the VCF file is
read.  Additionally, because phasing SNPs into haplotypes can cause some 
haplotypes to fail to pass this threshold, \code{VCF2RADdata} internally runs
\code{\link{MergeRareHaplotypes}} with 
\code{min.ind.with.haplotype = min.ind.with.minor.allele}, then
\code{\link{RemoveMonomorphicLoci}}, before returning the
final \code{RADdata} object.
}
\value{
A \code{\link{RADdata}} object.
}
\note{
In the \code{python} directory of the \pkg{polyRAD} installation, there is a
script called \code{tassel_vcf_tags.py} that can identify the full tag
sequence(s) for every allele imported by \code{VCF2RADdata}.
}
\references{
Variant Call Format specification: \url{http://samtools.github.io/hts-specs/}

TASSEL GBSv2 pipeline: \url{https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline}

GATK: \url{https://gatk.broadinstitute.org/hc/en-us}

Tassel4-Poly: \url{https://github.com/guilherme-pereira/tassel4-poly}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{MakeTasselVcfFilter}} for filtering to a smaller VCF file before
reading with \code{VCF2RADdata}.

To export to VCF: \code{\link{RADdata2VCF}}

Other data import functions:  \code{\link{readStacks}}, \code{\link{readHMC}},
\code{\link{readTagDigger}}, \code{\link{readTASSELGBSv2}},
\code{\link{readProcessIsoloci}}, \code{\link{readDArTag}}
}
\examples{
# get the example VCF installed with polyRAD
exampleVCF <- system.file("extdata", "Msi01genes.vcf", package = "polyRAD")

\donttest{
# loading VariantAnnotation namespace takes >10s,
# so is excluded from CRAN checks

require(VariantAnnotation)

# Compress and index the VCF before reading, if not already done
if(!file.exists(paste(exampleVCF, "bgz", sep = "."))){
  vcfBG <- bgzip(exampleVCF)
  indexTabix(vcfBG, "vcf")
}

# Read into RADdata object
myRAD <- VCF2RADdata(exampleVCF, expectedLoci = 100, expectedAlleles = 500)

# Example of subsetting by genomic region (first 200 kb on Chr01)
mysv <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD",
                     samples = samples(scanVcfHeader(exampleVCF)),
                     which = GRanges("01", IRanges(1, 200000)))
myRAD2 <- VCF2RADdata(exampleVCF, expectedLoci = 100, expectedAlleles = 500,
                      svparam = mysv, yieldSize = NA_integer_)
}
}

\keyword{ file }
