\name{get_distribution}
\alias{get_distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get_distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_distribution(matrix,FUNCTION,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{A matrix with samples on the rows and bases on the columns.

}
\item{FUNCTION}{A character specifyng the function to be applied on each column of \emph{matrix}. The possible values are: \emph{mean},\emph{max},\emph{min},\emph{median} and \emph{sum}.

}
\item{index}{index returned by \emph{get_heteroplasmy}. Default=NULL.

}

}
\details{

}
\value{
It returns a vector with length equal to n_col of \emph{matrix} where each element contains the result of the operation defined by \emph{FUNCTION}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gabriele Lubatti <gabriele.lubatti@helmholtz-muenchen.de>
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
