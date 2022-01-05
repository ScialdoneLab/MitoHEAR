\name{heatmap_plot}
\alias{heatmap_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
heatmap_plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{


heatmap_plot(marker_plot, marker_plot_plot, new_classification,
                       old_classification,res_ang_sqrt, cluster_columns = F, cluster_rows = T, name_legend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker_plot}{
Character vector. Samples to select for showing in the heatmap.
}
\item{marker_plot_plot}{
Character vector. Selected samples will be labelled with the elements of this vector in the heatmap.
}
\item{new_classification}{
Character vector with length equal to selected samples,corresponding to the clustering analysis based on heteroplasmy values. See \emph{clustering_dist_ang} for more info.
}
\item{old_classification}{
Character vector with length equal to selected samples,corresponding to the partition available before clustering on heteroplasmy values.
}
\item{res_ang_sqrt
}{Distance matrix obtained from \emph{clustering_dist_ang}
.
}
\item{cluster_columns
}{Logical. Parameter for cluster_columns argument of the function \emph{Heatmap} in the package \emph{ComplexHeatmap} }
\item{cluster_rows
}{Logical. Parameter for cluster_rows argument of the function \emph{Heatmap}}

\item{name_legend
}{Character value.Parameter for name argument of the function \emph{Heatmap} }

}
\details{

}
\value{
Heatmap plot produced by function \emph{Heatmap}

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
\link{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}}

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
