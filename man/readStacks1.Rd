\name{readStacks1}
\alias{readStacks1}
\title{
Import Read Depth from Stacks Version 1
}
\description{
Using the catalog files output by cstacks and matches file output by sstacks,
this function imports read depth into a \code{\link{RADdata}} object.  If
genomic alignments were used, alignment data can optionally be imported.
}
\usage{
readStacks1(allelesFile, matchesFolder, min.ind.with.reads = 200, 
            min.ind.with.minor.allele = 10, readAlignmentData = FALSE, 
            possiblePloidies = list(2), contamRate = 0.001)
}
\arguments{
  \item{allelesFile}{
Path to the "alleles" file from the Stacks catalog.
}
  \item{matchesFolder}{
Path to the folder containing "matches" files to import.
}
  \item{min.ind.with.reads}{
For filtering loci.  A locus must have at least this many samples with
reads in order to be retained.
}
  \item{min.ind.with.minor.allele}{
For filtering loci.  A locus must have at least this many samples with
reads for the minor allele in order to be retained.  For loci with more
than two alleles, at least two alleles must be present in at least this
many individuals.
}
  \item{readAlignmentData}{
If \code{TRUE}, the "tags" file from the Stacks catalog will be read,
and chromosome, position, and strand will be imported to the \code{locTable}
slot of the output.  It is assumed that the "tags" file is in the same
directory as the "alleles" file.
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
\code{\link{readHMC}}
}
\examples{
\dontrun{

# Assuming the working directory contains the catalog and all matches files:

myStacks <- readStacks1("batch_1.catalog.alleles.tsv", ".", 
                        readAlignmentData = TRUE)
}
}

\keyword{ file }
