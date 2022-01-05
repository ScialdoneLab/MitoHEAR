\name{get_wilcox_test}
\alias{get_wilcox_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get_wilcox_test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_wilcox_test(matrix,cluster,label_1,label_2,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.

}
\item{cluster}{Charachter vector specifying a partition of the samples, with length equal to n_row of matrix.
}
\item{label_1}{Character name of a first label included in cluster. It denotes the first group used for the Wilcoxon test
}
\item{label_2}{Character name of a second label included in cluster and different from label_1. it denotes the second group used for the Wilcoxon test.
}
\item{index}{index returned by \emph{get_heteroplasmy}. Default=NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
It returns a vector of length equal to n_row in matrix. Each element stands for a base and it contains the adjusted p-value (FDR), obtained in unpaired two-samples Wilcoxon test  from the comparison of the heteroplasmy between the label_1 and label_2 group.

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
\link{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test}
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
