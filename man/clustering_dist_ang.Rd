\name{clustering_dist_ang}
\alias{clustering_dist_ang}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
clustering_dist_ang
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clustering_dist_ang(heteroplasmy_matrix, allele_matrix, cluster, top_pos, deepSplit_param,minClusterSize_param,threshold=0.20,min_value,index=NULL,relevant_bases= NULL,max_frac=0.70)
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
\item{threshold}{
Numeric value. If a base has heteroplasmy greater or equal to \emph{threshold} in more than \emph{max_frac} of cells, then the base is not considered for down stream analysis.
}

\item{min_value}{
Numeric value. If relevant_bases=NULL, then the bases for performing hierarchical clustering are the ones whose relative variance (variance of the base divided sum of variance among \emph{top_pos} bases) is above \emph{min_value}.
}

\item{index}{
index returned by \emph{get_heteroplasmy}.
}



\item{relevant_bases}{
Character vector of bases to consider as features for performing hierarchical clustering on samples.Default=NULL.
}

\item{max_frac}{
Numeric value.If a base has heteroplasmy greater or equal to \emph{threshold} in more than \emph{max_frac} of cells, then the base is not considered for down stream analysis.
}

}
\details{
For each pair of samples and for each base, an angular distance matrix is computed based on the four allele frequencies. Then only the angular distances corresponding to the relevant_bases are kept. If relevant bases is NULL, then only the angular distances corresponding to the bases with relative distance variance among samples above \emph{min_value} are kept .
Finally the distance between each pair of samples is defined as the euclidean distance of the angular distances corresponding to the bases that pass the previous filtering step.
On this final distance matrix, a hierarchical clustering approach is performed using the function \emph{cutreeHybrid} of the package \emph{dynamicTreeCut}.


}
\value{
It returns a list with 4 elements:
\item{classe}{
Data frame with two columns and n_row equal to n_row in heteroplasmy_matrix. The first column is the old cluster annotation provided by cluster. The second columns is the new cluster annotation obtained with hierarichal clustering on distance matrix based on heteroplasmy values.

}

\item{res_ang_sqrt}{
Distance matrix based on heteroplasmy values as defined in the section \emph{Details}

}

\item{top_dist}{
Vector of bases used for hierarchical clustering. If \emph{relevant_bases} is not NULL, then top_dist=NULL

}

\item{common_idx}{
Vector of indeces of samples for which hierarchical clustering is performed. If \emph{index} is NULL, then common_idx=NULL

}

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
\link{https://www.rdocumentation.org/packages/dynamicTreeCut/versions/1.63-1/topics/cutreeHybrid}
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
