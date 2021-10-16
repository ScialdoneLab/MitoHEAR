\name{Plot_boxplot}
\alias{Plot_boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot_boxplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_boxplot(distribution_1,distribution_2,label_1,label_2,name_x,name_y,name_title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distribution_1, distribution 2}{
Numeric vector
}
\item{label 1}{
Character vector of length equal to distribution_1
}
\item{label 1}{
Character vector of length equal to distribution_2
}
\item{name_x}{
Character name specifyng the xlab argument in ggplot.
}
\item{name_y}{
Character name specifyng the ylab argument in ggplot.
}
\item{name_title}{
Character name specifyng the ggtitle argument in ggplot.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot boxplot of the quantities specified by \emph{distribution_1} and \emph{distribution_2}, divided by the conditions denoted by \emph{label_1} and \emph{label_2}.
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
