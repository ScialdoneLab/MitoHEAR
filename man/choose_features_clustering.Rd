\name{choose_features_clustering}
\alias{choose_features_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
choose_features_clustering

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
choose_features_clustering(heteroplasmy_matrix, allele_matrix, cluster, top_pos, deepSplit_param,minClusterSize_param,min_value,threshold,index=NULL,max_frac=0.70)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heteroplasmy_matrix}{
Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.
}
\item{allele_matrix}{
Matrix with allele frequencies values obtained from \emph{get_heteroplasmy}.
}
\item{cluster}{
Charachter vector specifying a partition of the samples, with length equal to n_row of allele_matrix.
}
\item{top_pos}{
Numeric value. Number of bases sorted with decreasing values of distance variance (see section \emph{Details} below) among samples. If relevant_bases=NULL, then the bases for performing hierarchical clustering are the ones whose relative variance (variance of the base divided sum of variance among \emph{top_pos} bases) is above \emph{min_value}.
}



\item{deepSplit_param}{
Integer value between 0 and 4 for the deepSplit parameter of the function \emph{cutreeHybrid}. See section \emph{Details} below.

}


\item{minClusterSize_param}{
Integer value specifyng the minClusterSize parameter of the function \emph{cutreeHybrid}. See section \emph{Details} below.
}

\item{min_value}{
Numeric vector. For each value in the vector, the function \emph{clustering_dist_ang} is run with parameter \emph{min_value} equal to one element of the vector \emph{min_value}.
}

\item{threshold}{
Numeric value. If a base has heteroplasmy greater or equal to \emph{threshold} in more than \emph{max_frac} of cells, then the base is not considered for down stream analysis.
}



\item{index}{
index returned by \emph{get_heteroplasmy}.
}



\item{max_frac}{
Numeric value.If a base has heteroplasmy greater or equal to \emph{threshold} in more than \emph{max_frac} of cells, then the base is not considered for down stream analysis.
}

}

\value{
Clustree plot.

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
\link{https://cran.r-project.org/web/packages/clustree/clustree.pdf}
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
