\name{reverseComplement}
\alias{reverseComplement}
\alias{reverseComplement,character-method}
\title{
Reverse Complement of DNA Sequence Stored as Character String
}
\description{
Whereas the \code{reverseComplement} function available in \pkg{Biostrings} only
functions on \code{XString} and \code{XStringSet} objects, the version in
\pkg{polyRAD} also works on character strings.  It is written as an S4 method
in order to avoid conflict with \pkg{Biostrings}.  It is primarily included for
internal use by \pkg{polyRAD}, but may be helpful at the user level as well.
}
\usage{
reverseComplement(x, ...)
}
\arguments{
  \item{x}{
A vector of character strings indicating DNA sequence using IUPAC codes.
}
  \item{\dots}{
Additional arguments (none implemented)
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A character vector.
}
\author{
Lindsay V. Clark
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{readDArTtag}} uses this function internally.
}
\examples{
reverseComplement(c("AAGT", "CCA"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
