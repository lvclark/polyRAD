\name{InbreedingFromHindHe}
\alias{InbreedingFromHindHe}
\title{
Estimate Inbreeding from Hind/He for a Given Ploidy
}
\description{
After running \code{\link{HindHe}} and examining the distribution of values
across taxa and loci, \code{InbreedingFromHindHe} can be used to estimate
the inbreeding statistic \eqn{F} from the median or mode value of
\eqn{H_{ind}/H_E}{Hind/He}.  The statistic estimated encompasses inbreeding
from all sources, including population structure, self-fertilization, and
preferential mating among relatives.  It is intended to be used as input to
the \code{process_isoloci.py} script.
}
\usage{
InbreedingFromHindHe(hindhe, ploidy)
}
\arguments{
  \item{hindhe}{
A value for \eqn{H_{ind}/H_E}{Hind/He}.  It should generally range from zero to
one.
}
  \item{ploidy}{
A single integer indicating the ploidy of the population.
}
}

\value{
A number indicating the inbreeding statistic \eqn{F}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Lindsay V. Clark
}

\seealso{
\code{\link{HindHe}}, \code{\link{readProcessSamMulti}}, \code{\link{readProcessIsoloci}}
}
\examples{
InbreedingFromHindHe(0.5, 2)
InbreedingFromHindHe(0.4, 2)
InbreedingFromHindHe(0.5, 4)
}
\keyword{ arith }
