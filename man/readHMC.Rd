\name{readHMC}
\alias{readHMC}
\title{
Import read depth from UNEAK
}
\description{
This function reads the \dQuote{HapMap.hmc.txt} and \dQuote{HapMap.fas.txt}
files output by the UNEAK pipeline and uses the data to generate a 
\dQuote{RADdata} object.
}
\usage{
readHMC(file, includeLoci = NULL, shortIndNames = TRUE, 
        possiblePloidies = list(2), taxaPloidy = 2L, contamRate = 0.001,
        fastafile = sub("hmc.txt", "fas.txt", file, fixed = TRUE))
}
\arguments{
  \item{file}{
Name of the file containing read depth (typically \dQuote{HapMap.hmc.txt}).
}
  \item{includeLoci}{
An optional character vector of loci to be included in the output.
}
  \item{shortIndNames}{
Boolean.  If TRUE, taxa names will be shortened with respect to those in the
file, eliminating all text after and including the first underscore.
}
  \item{possiblePloidies}{
A list of numeric vectors indicating potential inheritance modes of SNPs in the
dataset.  See \code{\link{RADdata}}.
}
  \item{taxaPloidy}{
A single integer, or an integer vector with one value per taxon, indicating
ploidy.  See \code{\link{RADdata}}.
}
  \item{contamRate}{
A number ranging from zero to one (typically small) indicating the expected 
rate of sample cross-contamination.
}
  \item{fastafile}{
Name of the file containing tag sequences (typically \dQuote{HapMap.fas.txt}).
}
}

\value{
A \code{\link{RADdata}} object containing read depth, taxa and locus names, and
nucleotides at variable sites.
}
\references{
Lu, F., Lipka, A. E., Glaubitz, J., Elshire, R., Cherney, J. H., Casler, M. D.,
Buckler, E. S. and Costich, D. E. (2013) Switchgrass genomic diversity, ploidy, 
and evolution: novel insights from a network-based SNP discovery protocol.  
\emph{PLoS Genetics} \bold{9}, e1003215.

\url{https://www.maizegenetics.net/tassel}

\url{https://tassel.bitbucket.io/TasselArchived.html}
}
\author{
Lindsay V. Clark
}
\note{
UNEAK is not able to report read depths greater than 127, which may be 
problematic for high depth data on polyploid organisms.  The UNEAK
pipeline is no longer being updated and is currently only available with
archived versions of TASSEL.
}

\seealso{
\code{\link{readTagDigger}}, \code{\link{VCF2RADdata}}, 
\code{\link{readStacks}}, \code{\link{readTASSELGBSv2}},
\code{\link{readDArTag}}
}
\examples{
# for this example we'll create dummy files rather than using real ones
hmc <- tempfile()
write.table(data.frame(rs = c("TP1", "TP2", "TP3"),
                       ind1_merged_X3 = c("15|0", "4|6", "13|0"),
                       ind2_merged_X3 = c("0|0", "0|1", "0|5"),
                       HetCount_allele1 = c(0, 1, 0),
                       HetCount_allele2 = c(0, 1, 0),
                       Count_allele1 = c(15, 4, 13),
                       Count_allele2 = c(0, 7, 5),
                       Frequency = c(0, 0.75, 0.5)), row.names = FALSE,
            quote = FALSE, col.names = TRUE, sep = "\t", file = hmc)
fas <- tempfile()
writeLines(c(">TP1_query_64",
             "TGCAGAAAAAAAACGCTCGATGCCCCCTAATCCGTTTTCCCCATTCCGCTCGCCCCATCGGAGT",
             ">TP1_hit_64",
             "TGCAGAAAAAAAACGCTCGATGCCCCCTAATCCGTTTTCCCCATTCCGCTCGCCCCATTGGAGT",
             ">TP2_query_64",
             "TGCAGAAAAACAACACCCTAGGTAACAACCATATCTTATATTGCCGAATAAAAAACAACACCCC",
             ">TP2_hit_64",
             "TGCAGAAAAACAACACCCTAGGTAACAACCATATCTTATATTGCCGAATAAAAAATAACACCCC",
             ">TP3_query_64",
             "TGCAGAAAACATGGAGAGGGAGATGGCACGGCAGCACCACCGCTGGTCCGCTGCCCGTTTGCGG",
             ">TP3_hit_64",
             "TGCAGAAAACATGGAGATGGAGATGGCACGGCAGCACCACCGCTGGTCCGCTGCCCGTTTGCGG"),
             fas)

# now read the data
mydata <- readHMC(hmc, fastafile = fas)

# inspect the results
mydata
mydata$alleleDepth
mydata$alleleNucleotides
row.names(mydata$locTable)
}

\keyword{ file }

