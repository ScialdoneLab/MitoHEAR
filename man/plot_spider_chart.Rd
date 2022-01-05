\name{plot_spider_chart}
\alias{plot_spider_chart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_spider_chart
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_spider_chart(name_base,cluster_ci,heteroplasmy_matrix_ci,index_ci)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name_base}{
Character name specyfing the base for which the spider chart is computed.
}
\item{cluster_ci}{
Character vector specyfing the partition for samples in \emph{heteroplasmy_matrix_ci}.
}
\item{heteroplasmy_matrix_ci}{
Third element returned by \emph{get_heteroplasmy}.
}
\item{index_ci}{
Fifth element returned by \emph{get_heteroplasmy}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
radarchart plot produced by function \emph{radarchart}.
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
\link{https://rdrr.io/cran/fmsb/man/radarchart.html}
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
