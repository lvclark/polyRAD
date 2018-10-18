\name{TestOverdispersion}
\alias{TestOverdispersion}
\alias{TestOverdispersion.RADdata}
\title{
Test the Fit of Read Depth to Beta Binomial Distribution
}
\description{
This function is intended to help the user select a value to pass to the
\code{overdispersion} parameter of \code{\link{AddGenotypeLikelihood}},
generally via pipeline functions such as \code{\link{IterateHWE}} or
\code{\link{PipelineMapping2Parents}}.
}
\usage{
TestOverdispersion(object, ...)

\method{TestOverdispersion}{RADdata}(object, to_test = seq(6, 20, by = 2), \dots)
}
\arguments{
  \item{object}{
A \code{\link{RADdata}} object.  Genotype calling does not need to have
been performed, although for mapping populations it might be helpful to
have done a preliminary run of \code{PipelineMapping2Parents} without
linkage.
}
\item{to_test}{
A vector containing values to test.  These are values that will potentially
be used for the \code{overdispersion} argument of a pipeline function.  They
should all be positive numbers.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}
\details{
If no genotype calling has been performed, a single iteration under HWE using
default parameters will be done.  \code{object$ploidyChiSq} is then examined
to determine the most common/most likely inheritance mode for the whole
dataset.  The alleles that are examined are only those where this
inheritance mode has the lowest chi-squared value.

Within this inheritance mode and allele set, genotypes are selected where the
posterior probability of having a single copy of the allele is at least 0.95.
Read depth for these genotypes is then analyzed.  For each genotype, a 
probability is calculated for the read depth ratio to deviate from the expected
ratio by at least that much under the beta-binomial distribution.  This test
is performed for each overdispersion value provided in \code{to_test}.
}
\value{
A list of the same length as \code{to_test}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    UseMethod("TestOverdispersion", object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
