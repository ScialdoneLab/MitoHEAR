\name{plot_heteroplasmy_variability}
\alias{plot_heteroplasmy_variability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_heteroplasmy_variability
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_heteroplasmy_variability(heteroplasmy_matrix_sc,cluster_sc,threshold=0.1,frac=FALSE,index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heteroplasmy_matrix_sc}{
Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.
}
\item{cluster_sc}{
Charachter vector specifying a partition of the samples.
}
\item{threshold}{
Numeric value.
}
\item{frac}{
Logical. If FALSE the absolute number of cells that have at least one base with heteroplasmy above \emph{threshold} are shown divided by \emph{cluster_sc}. If TRUE, then the fraction of cells are shown.
}
\item{index}{
Fifth element returned by function \emph{get_heteroplasmy}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot2 object
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
