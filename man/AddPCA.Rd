\name{AddPCA}
\alias{AddPCA}
\title{
Perform Principal Components Analysis on \dQuote{RADdata} Object
}
\description{
This function uses read depth ratios or posterior genotype probabilities
(the latter preferentially) as input data for principal components analysis.
The PCA scores are then stored in the \code{$PCA} slot of the
\code{"\link{RADdata}"} object.
}
\usage{
AddPCA(object, ...)
\method{AddPCA}{RADdata}(object, nPcs = 20, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{nPcs}{
The number of principal component axes to return.
}
  \item{\dots}{
Additional arguments to be passed to the \code{pca} function from the
\pkg{pcaMethods} BioConductor package.
}
}
\details{
The PPCA (probabalistic PCA) method from \pkg{pcaMethods} is used,
due to the high missing data rate that is typical of genotyping-by-sequencing
datasets.
}
\value{
A \code{"RADdata"} object identical to the one passed to the function, but with
a matrix added to the \code{$PCA} slot.  This matrix contains PCA scores, with
taxa in rows, and PC axes in columns.
}

\author{
Lindsay V. Clark
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load data
data(exampleRAD)
# do PCA
exampleRAD <- AddPCA(exampleRAD, nPcs = 3)

plot(exampleRAD$PCA[,1], exampleRAD$PCA[,2])
}

\keyword{ methods }

