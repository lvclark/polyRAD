\name{EstimateContaminationRate}
\alias{EstimateContaminationRate}
\alias{EstimateConteminationRate.RADdata}
\title{
Estimate Sample Contamination Using Blanks
}
\description{
Based on mean read depth at blank and non-blank taxa, estimate sample
cross-contamination and add that information to the \code{"\link{RADdata}"}
object.
}
\usage{
EstimateContaminationRate(object, ...)
\method{EstimateContaminationRate}{RADdata}(object, multiplier = 1, \dots)
}
\arguments{
  \item{object}{
A \code{"RADdata"} object where \code{\link{SetBlankTaxa}} has already
been used to assign one or more taxa as blanks.
}
\item{multiplier}{
A single numeric value, or a named numeric vector with one value per blank
taxon in \code{object}, with names matching the blank taxa names.  Read depth
at blank taxa will be multiplied by this number when estimating sample 
cross-contamination.  See example below.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
This function estimates sample cross-contamination assuming that the only
source of contamination is from adapter or sample spill-over between wells
during library preparation, or contamination among the libraries themselves.
If you anticipate a higher rate of contamination during DNA extraction before
library preparation, you may wish to increase the value using
\code{\link{SetContamRate}}.

It is important to set the contamination rate to a reasonably accurate value
(\emph{i.e.} the right order of magnitude) in order for \pkg{polyRAD} to
be able to identify homozygotes that may otherwise appear heterozygous due
to contamination.
}
\value{
A \code{"RADdata"} object identical to \code{object} but with the 
\code{"contamRate"} attribute adjusted.
}

\author{
Lindsay V. Clark
}


\examples{
# dataset for this example
data(Msi01genes)

# give the name of the taxon that is blank
Msi01genes <- SetBlankTaxa(Msi01genes, "blank")

# Fifteen libraries were done; blank is pooled over all of them, and
# most other samples are pooled over two libraries.
mymult <- 2/15

# estimate the contamination rate
Msi01genes <- EstimateContaminationRate(Msi01genes, multiplier = mymult)
}

\keyword{ manip }
