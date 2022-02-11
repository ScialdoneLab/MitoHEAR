#' get_distribution
#' @param FUNCTION A character specifying the function to be applied on each
#' column of \emph{matrix}. The possible values are:
#' \emph{mean},\emph{max},\emph{min},\emph{median} and \emph{sum}.
#' @param index index returned by \emph{get_heteroplasmy}.
#' @inheritParams plot_heteroplasmy
#' @return It returns a vector with length equal to n_col of \emph{matrix}
#' where each element contains the result of the operation defined by
#' \emph{FUNCTION}.
#' @author Gabriele Lubatti <gabriele.lubatti@@helmholtz-muenchen.de>
#' @export get_distribution
get_distribution <- function(heteroplasmy_matrix, FUNCTION, index=NULL) {
  base <- colnames(heteroplasmy_matrix)
  distribution <- rep(0, length(base))
   if (is.null(index)) {
    for (i in 1:length(base)) {
      Y <- as.vector(heteroplasmy_matrix[, base[i]])
      if (FUNCTION =='mean') {
        Y <- mean(Y)
        }
      if (FUNCTION =='max') {
        Y <- max(Y)
        }
      if (FUNCTION =='min') {
        Y <- min(Y)
        }
      if (FUNCTION =='median') {
        Y <- median(Y)
        }
      if (FUNCTION == 'sum') {
        Y <- sum(Y)
        }
      distribution[i] <- Y
      names(distribution)[i] <- base[i]
    }}

  else {
    for (i in 1:length(base)) {
      index_cell <- index[[which(names(index) ==  base[i])]]
      Y <- as.vector(heteroplasmy_matrix[index_cell, base[i]])
      if (FUNCTION ==  'mean') {
        Y <- mean(Y)
        }
      if (FUNCTION ==  'max') {
        Y <- max(Y)
        }
      if (FUNCTION ==  'min') {
        Y <- min(Y)
        }
      if (FUNCTION ==  'median') {
        Y <- median(Y)
        }
      if (FUNCTION ==  'sum') {
        Y <- sum(Y)
        }
      distribution[i] <- Y
      names(distribution)[i] <- base[i]
    }}
  return(distribution)
}








#' filter_bases
#' @param min_heteroplasmy Numeric value.
#' @param min_cells Numeric value.
#' @inheritParams plot_heteroplasmy
#' @return Character vector of bases that have an heteroplasmy greater than
#' \emph{min_heteroplasmy} in more than \emph{min_cells}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export filter_bases
filter_bases <- function(heteroplasmy_matrix, min_heteroplasmy, min_cells, index = NULL) {
  position = colnames(heteroplasmy_matrix)
  filter = rep(0, length(position))
   if (is.null(index) ) {
    for (i in 1:length(position) ) {
      Y = data.frame(t(heteroplasmy_matrix[, position[i]]))
      filter[i] = sum(Y>min_heteroplasmy)
    }
    heteroplasmy_matrix_select = heteroplasmy_matrix[, which(filter>min_cells)]
  }
  else{
    for (i in 1:length(position)) {
      Y = data.frame(t(heteroplasmy_matrix[as.numeric(index[[position[i]]]), position[i]]))
      filter[i] = sum(Y>min_heteroplasmy)
    }
    heteroplasmy_matrix_select = heteroplasmy_matrix[, which(filter>min_cells)]
    }
  return(colnames(heteroplasmy_matrix_select))
}








#' detect_insertion
#' @param ref_sequence Character vector whose elements are the bases of a DNA
#' sequence to use as reference.
#' @param different_sequence Character vector whose elements are the bases of a
#' DNA sequence different from the reference.
#' @param length_comparison Integer number. Number of bases to consider for the
#' comparison between the two DNA sequences in order to detect and remove
#' insertions in the non-reference sequence.
#' @return Character vector of the different_sequence with length equal to
#' ref_sequence, after having removed the insertions.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export detect_insertion
detect_insertion <- function(ref_sequence, different_sequence, length_comparison = 10) {
  max_length_insertion = length(different_sequence)-length(ref_sequence)
  for (i in 1:(length(ref_sequence)-length_comparison)) {
    if (length(different_sequence)>length(ref_sequence)) {
       if (different_sequence[i]!= ref_sequence[i]) {
        for(k in 1:max_length_insertion) {
           if (all(different_sequence[(i+k):(i+k+length_comparison)] ==  ref_sequence[i:(i+length_comparison)])) {
            different_sequence = different_sequence[-seq(i, (i+k-1))]
            break
          }
        }
      }}
    else{}
  }
  return(different_sequence)}


















