\name{Plot_distance_matrix}
\alias{Plot_distance_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_distance_matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_distance_matrix(dist_matrix_sc,cluster)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist_matrix_sc}{
Distance matrix returned by function \emph{clustering_dist_ang} (second element of the output).
}
\item{cluster}{
Character vector.Can be one of the two partitions returned by function \emph{clustering_dist_ang} (first element of the output).
}

}
\details{

}
\value{
Heatmap plot produced by function \emph{Heatmap}

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
\link{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}}

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
