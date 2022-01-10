\name{plot_correlation_bases}
\alias{plot_correlation_bases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_correlation_bases
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_correlation_bases(bases_example,index_ci,heteroplasmy_matrix_ci)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bases_example}{
Character vector specyfing the bases for which the Spearman correlation across samples is computed
}
\item{index_ci}{
Fifth element returned by \emph{get_heteroplasmy}.
}
\item{heteroplasmy_matrix_ci}{
Third element returned by \emph{get_heteroplasmy}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\link{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
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
