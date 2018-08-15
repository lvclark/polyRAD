\name{readStacks}
\alias{readStacks}
\title{
Import Read Depth from Stacks
}
\description{
Using the catalog files output by cstacks and matches file output by sstacks,
this function imports read depth into a \code{\link{RADdata}} object.  If
genomic alignments were used, alignment data can optionally be imported.
}
\usage{
readStacks(allelesFile, matchesFolder, version = 2,
           min.ind.with.reads = 200,
           min.ind.with.minor.allele = 10, readAlignmentData = FALSE, 
           sumstatsFile = "populations.sumstats.tsv",
           possiblePloidies = list(2), contamRate = 0.001)
}
\arguments{
  \item{allelesFile}{
Path to the "alleles" file from the Stacks catalog.
}
  \item{matchesFolder}{
Path to the folder containing "matches" files to import.
}
\item{version}{
Either the number 1 or 2, indicating the version of Stacks.
}

  \item{min.ind.with.reads}{
For filtering loci.  A locus must have at least this many samples with
reads in order to be retained.
}
  \item{min.ind.with.minor.allele}{
For filtering loci.  A locus must have at least this many samples with
reads for the minor allele in order to be retained.  For loci with more
than two alleles, at least two alleles must be present in at least this
many individuals.  This argument is also passed internally to the 
\code{min.ind.with.haplotype} argument of \code{\link{MergeRareHaplotypes}}
to consolidate reads from rare alleles.
}
  \item{readAlignmentData}{
If \code{TRUE} and \code{version = 1}, the "tags" file from the Stacks catalog will be read,
and chromosome, position, and strand will be imported to the \code{locTable}
slot of the output.  It is assumed that the "tags" file is in the same
directory as the "alleles" file.  If \code{TRUE} and \code{version = 2},
\code{sumstatsFile} will be used for import of chromosome and position data.
}
\item{sumstatsFile}{
The name of the file containing summary statistics for loci.  Ignored
unless \code{version = 2} and \code{readAlignmentData = TRUE}.
}
  \item{possiblePloidies}{
A list indicating possible inheritance modes in the dataset.
See \code{\link{RADdata}}.
}
  \item{contamRate}{
A number from 0 to 1 (generally very small) indicating the expected rate of
cross contamination between samples.
}
}

\value{
A \code{\link{RADdata}} object.
}
\references{
Stacks website: \url{http://catchenlab.life.illinois.edu/stacks/}

Rochette, N. and Catchen, J. (2017) Deriving genotypes from RAD-seq short-read data 
using Stacks. \emph{Nature Protocols} \bold{12}, 2640--2659.

Catchen, J., Hohenlohe, P. A., Bassham, S., Amores, A., and Cresko., W. A. (2013) 
Stacks: an analysis tool set for population genomics. \emph{Molecular Ecology}
\bold{22}, 3124--3140.

Catchen, J. M., Amores, A., Hohenlohe, P., Cresko, W., and Postlethwait, J. H. (2011)
Stacks: building and genotyping loci de novo from short-read sequences. 
\emph{G3: Genes, Genomes, Genetics} \bold{1}, 171--182.
}
\author{
Lindsay V. Clark
}
\note{
This function has been tested with output from Stacks 1.47.
}

\seealso{
\code{\link{VCF2RADdata}}, \code{\link{readTagDigger}},
\code{\link{readHMC}}, \code{\link{readTASSELGBSv2}}
}
\examples{
\dontrun{

# Assuming the working directory contains the catalog and all matches files:

myStacks <- readStacks("batch_1.catalog.alleles.tsv", ".", 
                       version = 1,
                       readAlignmentData = TRUE)
}
}

\keyword{ file }
