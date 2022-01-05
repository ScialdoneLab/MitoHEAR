\name{dpt_test}
\alias{dpt_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dpt_test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dpt_test(heteroplasmy_matrix,time,index=NULL,method="GAM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heteroplasmy_matrix}{Matrix with heteroplasmy values obtained from \emph{get_heteroplasmy}.
}
\item{time}{Vector of diffusion pseudo time,with length equal to n_row of heteroplasmy_matrix.
}
\item{index}{index returned by \emph{get_heteroplasmy}. Default=NULL.
}

\item{method}{Character name denoting the method to choose for assigning an adjusted p value to each of the bases. Can be one of GAM,pearson and spearman.
GAM: For each base, a GAM fit with formula z ~ lo(t) is performed between the heteroplasmy values (z) and the time (t). The p value from the table "Anova for Parametric Effects" is then assigned to the base.
pearson,spearman:for each base, a pearson or spearman correlation test is performed between the heteroplasmy values and the time . The p value obtained from the test is then assigned to the base.
In all the three possible  methods, all the p values are then corrected with the method FDR.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with 2 columns and number of rows equal to n_col in heteroplasmy_matrix. In the first column there are the names of the bases while in the second column there are the adjusted p value.
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
\link{https://www.rdocumentation.org/packages/gam/versions/1.20/topics/gam}
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
