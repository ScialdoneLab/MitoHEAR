\name{vi_comparison}
\alias{vi_comparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
vi_comparison
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vi_comparison(old_classification,new_classification,number_iter)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old_classification}{
Character vector. First column of the dataframe returned by function \emph{clustering_dist_ang} (first element of the output).
}
\item{new_classification}{
Character vector.Second column of the dataframe returned by function \emph{clustering_dist_ang} (first element of the output).
}

\item{number_iter}{
Integer value. Specify how many random partition are generated (starting from re-shuffle of labels in \emph{old_classification}).
}

}
\details{
We compute the variation of information (VI) between the partition provided by \emph{new_classification} and \emph{old_classification}.
The VI between a random partitions (obtained wth re-shuflle from original labels in \emph{old_classification}) and \emph{old_classification} is also computed.
A distribution of VI values from random partitions is built.
Finally, from the comparison with this distribution, an empirical p value is given to the VI of the unsupervised cluster analysis.
}
\value{

Empirical p value.

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
\link{https://www.rdocumentation.org/packages/mcclust/versions/1.0/topics/vi.dist}}

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
