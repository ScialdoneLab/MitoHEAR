\name{plot_coordinate_heteroplasmy}
\alias{plot_coordinate_heteroplasmy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_coordinate_heteroplasmy
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_coordinate_heteroplasmy(coordinate_dm,heteroplasmy_matrix_ci,index_ci,name_base)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinate_dm}{
Dataframe whit samples on the rows and coordinates names on the columns.
}

\item{heteroplasmy_matrix_ci}{
Third element returned by \emph{get_heteroplasmy}.
}

\item{index_ci}{
Fifth element returned by \emph{get_heteroplasmy}.
}

\item{name_base}{
Character name for the base for which we want to plot the expression level.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot2 object.
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
