\name{filter_bases}
\alias{filter_bases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filter_bases
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filter_bases(heteroplasmy_matrix,min_heteroplasmy,min_cells,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heteroplasmy_matrix}{
Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.
}
\item{min_heteroplasmy}{
Numeric value.
}
\item{min_cells}{
Numeric value.
}
\item{index}{
index returned by \emph{get_heteroplasmy}. Default=NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Character vector of bases that have an heteroplasmy greater than min_heteroplasmy in more than min_cells.
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
