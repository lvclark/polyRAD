\name{readProcessSamMulti}
\alias{readProcessSamMulti}
\title{
Import Preliminary Data to Determine Parameters for Isolocus Sorting
}
\description{
This function imports the files output by \code{process_sam_multi.py} to a
\code{"\link{RADdata}"} object so that \code{\link{HindHe}} can be run to
filter samples and determine optimal parameters for \code{process_isoloci.py}.
}
\usage{
readProcessSamMulti(alignfile,
                    depthfile = sub("align", "depth", alignfile),
                    expectedLoci = 1000,
                    min.ind.with.reads = 200,
                    min.ind.with.minor.allele = 10,
                    possiblePloidies = list(2), taxaPloidy = 2L,
                    contamRate = 0.001,
                    expectedAlleles = expectedLoci * 15,
                    maxLoci = expectedLoci)
}
\arguments{
  \item{alignfile}{
A file output by \code{process_sam_multi.py}, generally in the format
\code{prefix_1_align.csv}.
}
  \item{depthfile}{
A file output by \code{process_sam_multi.py}, generally in the format
\code{prefix_1_depth.csv}.
}
  \item{expectedLoci}{
The number of loci expected in the final object. The default, 1000, is fairly
small because this function is intended to be used for preliminary analysis only.
}
  \item{min.ind.with.reads}{
The minimum number of taxa with reads needed in order for a locus to be retained
in the output.
}
  \item{min.ind.with.minor.allele}{
The minimum number of taxa with the same minor allele needed in order for a
locus to be retained in the output.
}
  \item{possiblePloidies}{
A list indicating expected inheritance modes for markers.  See
\code{\link{RADdata}}.
}
  \item{taxaPloidy}{
A single integer, or an integer vector with one value per taxon, indicating
ploidy.  See \code{\link{RADdata}}.
}
  \item{contamRate}{
A number ranging from zero to one (although in practice probably less than 0.01)
indicating the expected sample cross-contamination rate.
}
  \item{expectedAlleles}{
The expected number of alleles in the dataset.
}
  \item{maxLoci}{
The maximum number of loci to import before ceasing to read the file.  Set to
\code{Inf} if you want to read the entire file.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \code{"RADdata"} object.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{readProcessIsoloci}}
}
\examples{
\dontrun{
myRAD <- readProcessSamMulti("mydata_2_align.csv")
}
}

\keyword{ file }
