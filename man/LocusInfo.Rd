\name{LocusInfo}
\alias{LocusInfo}
\alias{LocusInfo.RADdata}
\title{
Get Information about a Single Locus
}
\description{
This function returns, and optionally prints, information about a single locus
with a \code{\link{RADdata}} object, including alignment position, allele
sequences, and genes overlapping the site.
}
\usage{
LocusInfo(object, ...)
\method{LocusInfo}{RADdata}(object, locus, genome = NULL,
          annotation = NULL, verbose = TRUE, ...)
}
\arguments{
  \item{object}{
A \code{RADdata} object.
}
\item{locus}{
A character string indicating the name of the locus to display.  Alternatively,
a character string indicating the name of an allele, for which the 
corresponding locus will be identified.
}
\item{genome}{
An optional \code{\link[Rsamtools:FaFile-class]{FaFile}} or 
\code{\link[BSgenome:BSgenome-class]{BSgenome}}
object containing the reference genome sequence.
}
\item{annotation}{
An optional \code{\link[GenomicFeatures:TxDb-class]{TxDb}} object containing the genome annotation.
}
\item{verbose}{
If \code{TRUE}, results will be printed to the console.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
The locus name, allele names, and allele sequences are always returned (although
allele names are not printed with \code{verbose}).  If the chromosome and 
position are known, those are also returned and printed.  If \code{annotation}
is provided, the function will return and print genes that overlap the locus.
If \code{annotation} and \code{genome} are provided, the function will attempt
to identify any amino acid changes caused by the alleles, using
\code{\link[VariantAnnotation:predictCoding-methods]{predictCoding}} 
internally.  Identification of
amino acid changes will work if the \code{RADdata} object was created with
\code{VCF2RADdata} using the \code{refgenome} argument to fill in non-variable
sites, and/or if the alleles are only one nucleotide long.
}
\value{
A list containing:

\item{Locus }{The name of the locus.}
\item{Chromosome }{The chromosome name, if present.}
\item{Position }{The position in base pairs on the chromosome, if present.}
\item{Alleles }{Allele names for the locus.}
\item{Haplotypes }{Allele sequences for the locus, in the same order.}
\item{Frequencies }{Allele frequencies, if present, in the same order.}
\item{Transcripts }{Transcripts overlapping the locus, if an annotation was
provided but it wasn't possible to predict amino acid changes.}
\item{PredictCoding }{The output of \code{predictCoding}, if it was run.}
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link[GenomicFeatures]{makeTxDbFromGFF}}, \code{\link{GetLoci}}
}
\examples{
data(exampleRAD)
exampleRAD <- AddAlleleFreqHWE(exampleRAD)
loc2info <- LocusInfo(exampleRAD, "loc2")
}

\keyword{ utilities }
