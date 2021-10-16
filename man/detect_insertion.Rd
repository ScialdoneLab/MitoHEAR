\name{detect_insertion}
\alias{detect_insertion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
detect_insertion
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
detect_insertion(ref_sequence,different_sequence,length_comparison=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref_sequence}{
Character vector whose element are the bases of a DNA sequence to use as reference.
}
\item{different_sequence}{
Character vector whose element are the bases of a DNA sequence different from the reference.
}
\item{length_comparison}{
Integer number. Number of bases to consider for the comparison between the two DNA sequences in order to detect and remove insertions in the non-reference sequence.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Character vector of the *different_sequence* with length equal to *ref_sequence*, after having removed the insertions.
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
