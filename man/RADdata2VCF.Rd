\name{RADdata2VCF}
\alias{RADdata2VCF}
\title{
Export RADdata Genotypes to VCF
}
\description{
Converts genotype calls from \pkg{polyRAD} into VCF format. The user may send
the results directly to a file, or to a
\code{\link[VariantAnnotation:VCF-class]{CollapsedVCF}} for further manipulation.
}
\usage{
RADdata2VCF(object, file = NULL, asSNPs = TRUE, hindhe = TRUE,
            sampleinfo = data.frame(row.names = GetTaxa(object)),
            contigs = data.frame(row.names = unique(object$locTable$Chr)))
}
\arguments{
  \item{object}{
A \code{\link{RADdata}} object in which genotype calling has been performed.
It is also important for the data to have been imported in a way that preserves
variant positions (i.e. \code{\link{readProcessIsoloci}},
\code{\link{readTASSELGBSv2}},
\code{\link{VCF2RADdata}} using the \code{refgenome} argument).
}
  \item{file}{
An optional character string or connection indicating where to write the file.
Append mode may be used with connections if multiple RADdata objects need to be
written to one VCF.
}
  \item{asSNPs}{
Boolean indicating whether to convert haplotypes to individual SNPs and indels.
}
  \item{hindhe}{
Boolean indicating whether to export a mean value of Hind/He
(see \code{\link{HindHe}}) for every sample and locus.
}
\item{sampleinfo}{
A data frame with optional columns indicating any sample metadata to export
to "SAMPLE" header lines.
}
\item{contigs}{
A data frame with optional columns providing information about contigs to
export to "contig" header lines.
}
}
\details{
Currently, the FORMAT fields exported are GT (genotype), AD (allelic read depth),
and DP (read depth).  Genotype posterior probabilities are not exported due to
the mathematical intractability of converting pseudo-biallelic probabilities to
multiallelic probabilities.

Genotypes exported to the GT field are obtained internally using
\code{\link{GetProbableGenotypes}}.

INFO fields exported include the standard fields NS (number of samples with more
than zero reads) and DP (total depth across samples) as well as the custom fields
LU (index of the marker in the original RADdata object) and HH (Hind/He statistic
for the marker).
}
\value{
A \code{\link[VariantAnnotation:VCF-class]{CollapsedVCF}} object.
}
\references{
\url{https://samtools.github.io/hts-specs/VCFv4.3.pdf}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{VCF2RADdata}}
}
\examples{
# Set up example dataset for export.
# You DO NOT need to adjust attr or locTable in your own dataset.
data(exampleRAD)
attr(exampleRAD$alleleNucleotides, "Variable_sites_only") <- FALSE
exampleRAD$locTable$Ref <- 
  exampleRAD$alleleNucleotides[match(1:nLoci(exampleRAD), exampleRAD$alleles2loc)]
exampleRAD <- IterateHWE(exampleRAD)

# An optional table of sample data
sampleinfo <- data.frame(row.names = GetTaxa(exampleRAD),
                         Population = rep(c("North", "South"), each = 50))

# Add contig information (fill in with actual data rather than random)
mycontigs <- data.frame(row.names = c("1", "4", "6", "9"), length = sample(1e8, 4),
                        URL = rep("ftp://mygenome.com/mygenome.fa", 4))

# Set up a file destination for this example
# (It is not necessary to use tempfile with your own data)
outfile <- tempfile(fileext = ".vcf")

# Export VCF
testvcf <- RADdata2VCF(exampleRAD, file = outfile, sampleinfo = sampleinfo,
                       contigs = mycontigs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ file }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
