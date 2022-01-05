\name{plot_heteroplasmy}
\alias{plot_heteroplasmy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_heteroplasmy
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_heteroplasmy(position,heteroplasmy_matrix,cluster,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{Character name of the base to plot.
}
\item{heteroplasmy_matrix}{Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.
}
\item{cluster}{Charachter vector specifying a partition of the samples, with length equal to n_row of heteroplasmy_matrix.
}
\item{index}{index returned by \emph{get_heteroplasmy}. Default=NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot object of the heteroplasmy level of a specific base across samples divided according to cluster.
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
