\name{Plot_distribution}
\alias{Plot_distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_distribution(quantity_counts_cell,x_name,title_name)
}

\arguments{
  \item{quantity_counts_cell}{numeric vector returned by \emph{get_distribution.}}
  \item{x_name}{Character name specifyng the xlab argument in ggplot.}
  \item{title_name}{Character name specifyng the ggtitle argument in ggplot..}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot density plot of the numeric vector \emph{quantity_counts_cell}
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
