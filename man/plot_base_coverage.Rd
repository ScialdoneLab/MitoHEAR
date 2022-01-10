\name{plot_base_coverage}
\alias{plot_base_coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_base_coverage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_base_coverage(sum_matrix,sum_matrix_qc,selected_cells,interactive=FALSE,text_size=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sum_matrix}{
First element returned by the function \emph{get_heteroplasmy}.
}
\item{sum_matrix_qc}{
Second element returned by the function \emph{get_heteroplasmy}.
}
\item{selected_cells}{
Character vector with cells used for plotting the coverage.
}
\item{interactive}{
Logical. If TRUE an interactive plot is produced.
}
\item{text_size}{
Integer value specifying the size of the text in the plot.
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
