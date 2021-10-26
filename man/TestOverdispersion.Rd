\name{TestOverdispersion}
\alias{TestOverdispersion}
\alias{TestOverdispersion.RADdata}
\title{
Test the Fit of Read Depth to Beta-Binomial Distribution
}
\description{
This function is intended to help the user select a value to pass to the
\code{overdispersion} argument of \code{\link{AddGenotypeLikelihood}},
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
two-tailed probability is calculated for the read depth ratio to deviate from 
the expected ratio by at least that much under the beta-binomial distribution.
This test is performed for each overdispersion value provided in
\code{to_test}.
}
\value{
A list of the same length as \code{to_test} plus one.  The names of the list are
\code{to_test} converted to a character vector.  Each item in the list is a
vector of p-values, one per examined genotype, of the read depth ratio for
that genotype to deviate that much from the expected ratio.  The last item,
named "optimal", is a single number indicating the optimal value for the
overdispersion parameter based on the p-value distributions.  If the optimal
value was the minimum or maximum tested, \code{NA} is returned in the
\code{"optimal"} slot to encourage the user to test other values.
}

\author{
Lindsay V. Clark
}

\examples{
# dataset with overdispersion
data(Msi01genes)

# test several values for the overdispersion parameter
myP <- TestOverdispersion(Msi01genes, to_test = 8:10)

# view results as quantiles
sapply(myP[names(myP) != "optimal"],
       quantile, probs = c(0.01, 0.25, 0.5, 0.75, 0.99))
}

\keyword{ distribution }
\keyword{ methods }
