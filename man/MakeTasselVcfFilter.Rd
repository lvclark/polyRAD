\name{MakeTasselVcfFilter}
\alias{MakeTasselVcfFilter}
\title{
Filter Lines of a VCF File By Call Rate and Allele Frequency
}
\description{
This function creates another function that can be used as a prefilter
by the function \code{filterVcf} in the package \pkg{VariantAnnotation}.
The user can set a minimum number of indiviuals with reads and a minimum
number of individuals with the minor allele (either the alternative or
reference allele).  The filter can be used to generate a smaller VCF file 
before reading with \code{\link{VCF2RADdata}}.
}
\usage{
MakeTasselVcfFilter(min.ind.with.reads = 200, min.ind.with.minor.allele = 10)
}
\arguments{
  \item{min.ind.with.reads}{
An integer indicating the minimum number of individuals that must have
reads in order for a marker to be retained.
}
  \item{min.ind.with.minor.allele}{
An integer indicating the minimum number of individuals that must have
the minor allele in order for a marker to be retained.
}
}
\details{
This function assumes the VCF file was output by the TASSEL GBSv2
pipeline.  This means that each genotype field begins with two digits
ranging from zero to three
separated by a forward slash to indicate the called genotype, followed
by a colon.
}
\value{
A function is returned.  The function takes as its only argument a character
vector representing a set of lines from a VCF file, with each line representing
one SNP.  The function returns a logical vector the same length as the
character vector, with \code{TRUE} if the SNP meets the threshold for call rate
and minor allele frequency, and \code{FALSE} if it does not.
}
\references{
\url{https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline}
}
\author{
Lindsay V. Clark
}

\examples{
require(VariantAnnotation)
# get the example VCF installed with polyRAD
exampleVCF <- system.file("extdata", "Msi01genes.vcf", package = "polyRAD")

# zip and index the file using Tabix
tempbgz <- tempfile(fileext = ".bgz") # temporary file for example
bgzip(exampleVCF, dest = tempbgz)
indexTabix(tempbgz, format = "vcf")

# make the filtering function
filterfun <- MakeTasselVcfFilter(300, 15)

# filter to a new file
tempvcf <- tempfile(fileext = ".vcf") # temporary file for example
filterVcf(tempbgz, destination = tempvcf, 
          prefilters = FilterRules(list(filterfun)))
}

\keyword{ file }
