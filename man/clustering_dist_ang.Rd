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
clustering_dist_ang(heteroplasmy_matrix,allele_matrix,cluster,top_pos,deepSplit_param,minClusterSize_param,relevant_bases=NULL)
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
Charachter vector specifying a partion of the samples, with length equal to n_row of allele_matrix.
}
\item{top_pos}{
Numeric value. Top Number of bases with highest distance variance (see section \emph{Details} below) among samples. If relevant_bases=NULL, then these top number of bases are considered as features for performing hierarchical clustering on samples.
}
\item{deepSplit_param}{
Integer value between 0 and 4 for the deepSplit parameter of the function \emph{cutreeHybrid}. See section \emph{Details} below.

}
\item{minClusterSize_param}{
Integer value specifyng the minClusterSize parameter of the function \emph{cutreeHybrid}. See section \emph{Details} below.
}

\item{relevant_bases}{
Character vector of bases to consider as features for performing hierarchical clustering on samples.Default=NULL.
}

}
\details{
For each pair of samples and for each base, an angular distance matrix is computed based on the four allele frequencies. Then only the angular distances corresponding to the relevant_bases are kept. If relevant bases is NULL, then only the angular distances corresponding to the bases with highest distance variance among samples are kept.
Finally the distance between each pair of samples is defined as the euclidean distance of the angular distances corresponding to the bases that pass the previous filtering step.
On this final distance matrix, a hierarchical clustering approach is performed using the function \emph{cutreeHybrid} of the package \emph{dynamicTreeCut}.


}
\value{
It returns a list with 2 elements:
\item{classe}{
Data frame with two columns and n_row equal to n_row in heteroplasmy_matrix. The first column is the old cluster annotation provided by cluster. The second columns is the new cluster annotation obtained with hierarichal clustering on distance matrix based on heteroplasmy values.

}

\item{res_ang_sqrt}{
Distance matrix based on heteroplasmy values as defined in the section \emph{Details}

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
