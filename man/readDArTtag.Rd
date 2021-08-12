\name{readDArTtag}
\alias{readDArTtag}
\title{
Import Data from DArT Sequencing
}
\description{
\href{https://www.diversityarrays.com/}{Diversity Array Technologies (DArT)}
provides a tag-based genotyping-by-sequencing service. Together with
\href{https://breedinginsight.org/}{Breeding Insight}, a format was
developed indicting haplotype sequence and read depth, and that format is
imported by this function to make a \code{\link{RADdata}} object. The target
SNP and all off-target SNPs within the amplicon are imported as haplotypes.
Because the file format does not indicate strandedness of the tag, BLAST
results are used so that sequence and position are accurately stored in the
\code{RADdata} object.  See the \dQuote{extdata} folder of the \pkg{polyRAD}
installation for example files.
}
\usage{
readDArTtag(file, botloci = NULL, blastfile = NULL, excludeHaps = NULL,
            includeHaps = NULL, n.header.rows = 7, sample.name.row = 7,
            trim.sample.names = "_[^_]+_[ABCDEFGH][[:digit:]][012]?$",
            sep.counts = ",", sep.blast = "\t", possiblePloidies = list(2),
            contamRate = 0.001)
}
\arguments{
  \item{file}{
The file name of a spreadsheet from DArT indicating haplotype sequence and read
depth.
}
  \item{botloci}{
A character vector indicating the names of loci for which the sequence is on the
bottom strand with respect to the reference genome.  All other loci are assumed
to be on the top strand.  Only one of \code{blastfile} and \code{botloci} should
be provided.
}
  \item{blastfile}{
File name for BLAST results for haplotypes. The file should be in tabular format
with \code{qseqid}, \code{sseqid}, \code{sstart}, \code{send}, and \code{pident}
columns, indicated with column headers.  Only one of \code{blastfile} and
\code{botloci} should be provided.
}
  \item{excludeHaps}{
Optional.  Character vector with names of haplotypes (from the \dQuote{AlleleID}
column) that should not be imported.  Should not be used if \code{includeHaps}
is provided.
}
  \item{includeHaps}{
Optional.  Character vector with names of haplotypes (from the \dQuote{AlleleID}
column) that should be imported.  Should not be used if \code{excludeHaps} is
provided.
}
  \item{n.header.rows}{
Integer.  The number of header rows in \code{file}, not including the full row
of column headers.
}
  \item{sample.name.row}{
Integer.  The row within \code{file} from which sample names should be derived.
}
  \item{trim.sample.names}{
A regular expression indicating text to trim off of sample names.  Use \code{""}
if no trimming should be performed.
}
  \item{sep.counts}{
The field separator character for \code{file}.  The default assumes CSV.
}
  \item{sep.blast}{
The field separator character for the BLAST results.  The default assumes
tab-delimited.
}
  \item{possiblePloidies}{
A list indicating possible inheritance modes.  See \code{\link{RADdata}}.
}
  \item{contamRate}{
Expected sample cross-contamination rate. See \code{\link{RADdata}}.
}
}
\details{
The \dQuote{CloneID} column is used for locus names, and is assumed to contain
the chromosome (or scaffold) name and position, separated by an underscore.
The position is assumed to refer to the target SNP, which is identified by
comparing the \dQuote{Ref_001} and \dQuote{Alt_002} sequences.  The position
is then converted to refer to the beginning of the tag (which may have been
reverse complemented depending on BLAST results), since additional SNPs may
be present.  This facilitates accurate export to VCF using
\code{\link{RADdata2VCF}}.

Column names for the BLAST file can be \dQuote{Query}, \dQuote{Subject},
\dQuote{S_start}, \dQuote{S_end}, and \dQuote{\%Identity}, for compatibility
with Breeding Insight formats.
}
\value{
A \code{\link{RADdata}} object ready for QC and genotype calling. Assuming
the \dQuote{Ref_001} and \dQuote{Alt_002} alleles were not excluded, the
\code{locTable} slot will include columns for chromosome, position, strand, and
reference sequence.
}
\references{
\url{https://www.diversityarrays.com/}

\url{https://breedinginsight.org/}
}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{reverseComplement}}

\code{\link{readTagDigger}}, \code{\link{VCF2RADdata}}, 
\code{\link{readStacks}}, \code{\link{readTASSELGBSv2}},
\code{\link{readHMC}}

\code{\link{RADdata2VCF}}
}
\examples{
# Example files installed with polyRAD
dartfile <- system.file("extdata", "DArTtag_example.csv", package = "polyRAD")
blastfile <- system.file("extdata", "DArTtag_BLAST_example.txt",
                         package = "polyRAD")

# One haplotype doesn't seem to have correct alignment (see BLAST results)
exclude_hap <- c("Chr1_30668472|RefMatch_004")

# Import data
mydata <- readDArTtag(dartfile, blastfile = blastfile,
                      excludeHaps = exclude_hap,
                      possiblePloidies = list(4))
}

\keyword{ file }
