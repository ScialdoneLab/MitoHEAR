\name{plot_coverage}
\alias{plot_coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_coverage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
plot_coverage(biomart_file,fastaFile,chr_name,heteroplasmy_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biomart_file}{Character string with full path to the txt file downloaded from BioMart (\link{https://m.ensembl.org/info/data/biomart/index.html}) . It must have the following five columns:Gene.stable.ID, Gene.name, Gene.start..bp., Gene.end..bp., Chromosome.scaffold.name
}
\item{fastaFile}{Character string with full path to the fasta file of the genomic region of interest. It should be the same file used in  \emph{get_raw_counts_allele}.
}

\item{chr_name}{Character specifyng the name of the chromosome of interest. It must be one of the names in the \emph{Chromosome.scaffold.name} column from the \emph{biomart_file}.
}

\item{heteroplasmy_matrix}{heteroplasmy_matrix as returned by \emph{get_heteroplasmy}.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plot as returned by \emph{karyoploteR} package.
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
\link{http://bioconductor.org/packages/release/bioc/html/karyoploteR.html}
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
