\name{get_raw_counts_allele}
\alias{get_raw_counts_allele}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get_raw_counts_allele
}
\description{
It is one the two main function of the \strong{MitoHEAR} package (together with \emph{get_heteroplasmy}).
The function allows to obtain a matrix of counts (n_row = number of sample, n_col= 4*number of bases) of the four alleles in each base, for every sample. It takes as input a vector of sorted bam files (one bam file for each sample) and a fasta file for the genomic region of interest. It is based on the \emph{pileup} function of the package Rsamtools.
}
\usage{
get_raw_counts_allele(bam_input,path_fasta,cell_names,cores_number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_input}{
Character vector of sorted bam files (full path). Each sample is defined by one bam file. For each bam file it is needed also the index bam file (.bai) at the same path.
}
\item{path_fasta}{
Character string with full path to the fasta file of the genomic region of interest.
}
\item{cell_names}{
Character vector of sample names.
}

\item{cores_number}{
Number of cores to use.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with three elements

\item{matrix_allele_counts }{Matrix of counts (n_row = number of sample, n_col= 4*number of bases) of the four alleles in each base, for every sample. The row names is equal to cell_names. }
\item{name_position_allele}{Character vector with length equal to n_col of matrix_allele_counts. Each element specifies the coordinate of genomic position for a base and the allele.}
\item{name_position}{Character vector with length equal to n_col of matrix_allele_counts. Each element specifies the coordinate of genomic position for a base.}

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
\link{https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/pileup}

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
