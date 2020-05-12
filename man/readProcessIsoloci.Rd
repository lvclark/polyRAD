\name{readProcessIsoloci}
\alias{readProcessIsoloci}
\title{
Import Read Depth from Output of process_isoloci.py
}
\description{
After process_isoloci.py is used to assign RAD tags to alignment locations
within a highly duplicated genome, \code{readProcessIsoloci} imports the
resulting CSV to a \code{"\link{RADdata}"} object.
}
\usage{
readProcessIsoloci(sortedfile, min.ind.with.reads = 200,
                   min.ind.with.minor.allele = 10,
                   min.median.read.depth = 10,
                   possiblePloidies = list(2),
                   contamRate = 0.001,
                   nameFromTagStart = TRUE, mergeRareHap = TRUE)
}
\arguments{
  \item{sortedfile}{
File path to a CSV output by process_isoloci.py.
}
  \item{min.ind.with.reads}{
Minimum number of individuals with reads needed to retain a locus.
}
  \item{min.ind.with.minor.allele}{
Minimum number of individuals with reads in a minor allele needed to retain
a locus.
}
  \item{min.median.read.depth}{
Minimum median read depth across individuals (including individuals with depth 0)
needed to retain a locus.
}
  \item{possiblePloidies}{
A list indicating possible inheritance modes of loci.  See \code{\link{RADdata}}.
}
  \item{contamRate}{
Approximate rate of cross-contamination among samples.
}
  \item{nameFromTagStart}{
If \code{TRUE} loci will be named based on the alignment position and strand
of the RAD tag itself.  If \code{FALSE}, loci will be named based on the leftmost
position of the variable region of the RAD tag.  In either case,
\code{locTable$Pos} within the output will indicate the position of the variable
region of the tag.
}
  \item{mergeRareHap}{
Boolean indicating whether to run \code{\link{MergeRareHaplotypes}} after
building the \code{"RADdata"} object.
}
}
\details{
\code{\link{MergeIdenticalHaplotypes}} is used internally by this function to
merge alleles with identical sequence for the region shared by all tags, in
cases where tags vary in length within a locus.
}
\value{
A \code{"RADdata"} object containing read depth and alignment positions
from \code{sortedfile}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{readProcessSamMulti}}
}
%\examples{

%}

\keyword{ file }
