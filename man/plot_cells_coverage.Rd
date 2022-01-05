\name{plot_cells_coverage}
\alias{plot_cells_coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_cells_coverage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_cells_coverage(sum_matrix,cells_selected,cluster_ci,interactive=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sum_matrix}{
First element returned by the function \emph{get_heteroplasmy}.
}
\item{cells_selected}{
Character vector of cells for which the coverage is computed.
}
\item{cluster_ci}{
Character vector with partition information for cells specified in \emph{cells_selected}
}
\item{interactive}{
Logical. If TRUE an interactive plot is produced.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot2 object (if interactive=FALSE) or plotly object (if interactive=TRUE)
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
\link{https://plotly.com/r/}
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
