\name{plot_dpt}
\alias{plot_dpt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_dpt
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_dpt(position,heteroplasmy_matrix,cluster,time,gam_fit_result,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{Character name of the base to plot.
}
\item{heteroplasmy_matrix}{Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.
}

\item{cluster}{Charachter vector specifying a partition of the samples, with length equal to n_row of heteroplasmy_matrix.
}

\item{time}{Vector of diffusion pseudo time,with length equal to n_row of heteroplasmy_matrix.
}
\item{gam_fit}{Data frame returned by \emph{dpt_test}.
}
\item{index}{index returned by \emph{get_heteroplasmy}. Default=NULL.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot object of the heteroplasmy level of a specific base across samples and the GAM fitted curve. The title shows the adjusted p value (FDR) for the position obtained from \emph{get_heteroplasmy}.
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
\link{https://cran.r-project.org/web/packages/gam/gam.pdf}
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
