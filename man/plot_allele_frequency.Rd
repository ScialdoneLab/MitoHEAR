\name{plot_allele_frequency}
\alias{plot_allele_frequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_allele_frequency
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_allele_frequency(position,heteroplasmy_matrix,allele_matrix,cluster,names_allele_qc,names_position_qc,size_text,title,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{Character name of the base to plot.

}
\item{heteroplasmy_matrix}{Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.

}
\item{allele_matrix}{Matrix with allele frequencies values obtained from \emph{get_heteroplasmy}.

}
\item{cluster}{Charachter vector specifying a partion of the samples, with length equal to n_row of allele_matrix.

}
\item{names_allele_qc}{Character vector with length equal to n_col of allele_matrix. Each element specifies the name of the base and the allele.
}

\item{names_position_qc}{Character vector with length equal to n_col of allele_matrix. Each element specifies the name of the base.
}

\item{size_tetx}{Character specifying the size of the text for ggplot
}

\item{title}{Character specifying the top argument in the function \emph{grid.arrange})
}
\item{index}{index returned by \emph{get_heteroplasmy}. Default=NULL.
}
}
\details{

}
\value{
\emph{grid.arrange} plot of allele frequencies of a specific base across samples divided according to cluster
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
{\link{https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html}
}}
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
