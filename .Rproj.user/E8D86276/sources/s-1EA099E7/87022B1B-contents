#' get_heteroplasmy
#'
#' It is one of the two main functions of the \strong{MitoHEAR} package
#' (together with \emph{get_raw_counts_allele}). It computes the allele
#' frequencies and the heteroplasmy matrix starting from the counts matrix
#' obtained with \emph{get_raw_counts_allele}.
#'
#' Starting from raw counts allele matrix, the function performed two
#' consequentially filtering steps. The first one is on the samples, keeping
#' only the ones that cover a number of bases above number_positions. The
#' second one is on the bases, defined by the parameter filtering. The
#' heteroplasmy for each sample-base pair is computed as \emph{1-max(f)}, where
#' \emph{f} are the frequencies of the four alleles.
#'
#' @param raw_counts_allele A raw counts matrix obtained from
#' \emph{get_raw_counts_allele}.
#' @param name_position_allele A character vector with elements specifyng the
#' genomic coordinate of the base and the allele (obtained from
#' \emph{get_raw_counts_allele}).
#' @param name_position A character vector with elements specifyng the genomic
#' coordinate of the base (obtained from \emph{get_raw_counts_allele}).
#' @param number_reads Integer specifying the minimun number of counts above
#' which we consider the base covered by the sample.
#' @param number_positions Integer specifying the minimun number of bases that
#' must be covered by the sample (with counts>number_reads), in order to keep
#' the sample for down-stream analysis.
#' @param filtering Numeric value equal to 1 or 2. If 1 then only the bases
#' that are covered by all the samples are kept for the downstream analysis. If
#' 2 then all the bases that are covered by more than 50\% of the the samples
#' in each cluster (specified by my.clusters) are kept for the down-stream
#' analysis. Default is 1.
#' @param my.clusters Charachter vector specifying a partition of the samples.
#' It is only used when filtering is equal to 2. Deafult is NULL
#' @return It returns a list with 5 elements:
#'
#' \item{sum_matrix}{A matrix (n_row=number of sample, n_col=number of bases)
#' with the counts for each sample/base, for all the initial samples and bases
#' included in the raw counts allele matrix.}
#'
#' \item{sum_matrix_qc}{A matrix (n_row=number of sample, n_col=number of
#' bases) with the counts for each sample/base, for all the samples and bases
#' that pass the two consequentially filtering steps.}
#'
#' \item{heteroplasmy_matrix}{A matrix with the same dimension of sum_matrix_qc
#' where each entry (i,j) is the heteroplasmy for sample i in base j.}
#'
#' \item{norm_counts_allele_filter_qc}{A matrix (n_row=number of sample,
#' n_col=4*number of bases) with allele frequencies, for all the samples and
#' bases that pass the two consequentially filtering steps.}
#'
#' \item{index}{If filtering is equal to 2: indices of all the samples that
#' cover a base, for all bases and samples that pass the two consequentially
#' filtering steps; if filtering equal to 1: NULL }
#' @author Gabriele Lubatti <gabriele.lubatti@@helmholtz-muenchen.de>
#' @examples
#'
#' load(system.file("extdata", "after_qc.Rda", package = "MitoHEAR"))
#' load(system.file("extdata", "output_SNP_mt.Rda", package = "MitoHEAR"))
#' row.names(after_qc)=after_qc$new_name
#' cells_fmk_epi=after_qc[(after_qc$condition=="Cell competition OFF")&
#' (after_qc$cluster==1|after_qc$cluster==3|after_qc$cluster==4),"new_name"]
#' after_qc_fmk_epi=after_qc[cells_fmk_epi,]
#' my.clusters=after_qc_fmk_epi$cluster
#'
#' matrix_allele_counts=output_SNP_mt[[1]]
#' name_position_allele=output_SNP_mt[[2]]
#' name_position=output_SNP_mt[[3]]
#' epiblast_ci=get_heteroplasmy(matrix_allele_counts[cells_fmk_epi,],
#' name_position_allele,name_position,number_reads=50,number_positions=2000,
#' filtering = 2,my.clusters)
#'
#'
#' sum_matrix=epiblast_ci[[1]]
#' sum_matrix_qc=epiblast_ci[[2]]
#' heteroplasmy_matrix_ci=epiblast_ci[[3]]
#' allele_matrix_ci=epiblast_ci[[4]]
#' index=epiblast_ci[[5]]
#'
#' @export get_heteroplasmy
get_heteroplasmy=function (raw_counts_allele, name_position_allele, name_position,
                           number_reads, number_positions, filtering = 1, my.clusters = NULL)
{
  raw_counts_allele_filter = raw_counts_allele
  sum_matrix = matrix(0, ncol = length(unique(name_position)),
                      nrow = length(row.names(raw_counts_allele_filter)))
  colnames(sum_matrix) = unique(name_position)
  row.names(sum_matrix) = row.names(raw_counts_allele_filter)
  colnames(raw_counts_allele_filter) = name_position
  for (i in unique(name_position)) {
    sum_matrix[, i] = apply(raw_counts_allele_filter[, colnames(raw_counts_allele_filter) %in%
                                                       i], 1, function(x) {
                                                         return(sum(x))
                                                       })
  }
  sum_matrix = as.matrix(sum_matrix)
  pos_cov = apply(sum_matrix, 1, function(x) {
    sum(x > number_reads)
  })
  pos_cov_20000 = pos_cov[pos_cov > number_positions]
  sum_matrix_qc = sum_matrix[pos_cov > number_positions, ]
  if (filtering == 1) {
    raw_counts_allele_update = apply(sum_matrix_qc, 2, function(x) {
      if (sum(x > number_reads) == length(row.names(sum_matrix_qc))) {
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    })
    sum_matrix_qc = sum_matrix_qc[, raw_counts_allele_update]
    index = NULL
  }
  if (filtering == 2) {
    my.clusters = my.clusters[pos_cov > number_positions]
    cluster_unique = unique(my.clusters)
    index = as.list(rep(0, length(unique(my.clusters))))
    for (i in 1:length(cluster_unique)) {
      index[[i]] = which(my.clusters == cluster_unique[i])
    }
    raw_counts_allele_update = apply(sum_matrix_qc, 2, function(x) {
      cell_cluster = as.list(rep(0, length(unique(my.clusters))))
      selection_base = as.list(rep(0, length(unique(my.clusters))))
      for (i in 1:length(cluster_unique)) {
        cell_cluster[[i]] = x[index[[i]]]
        if (sum(cell_cluster[[i]] > number_reads) > (0.5 *
                                                     length(cell_cluster[[i]]))) {
          selection_base[[i]] = TRUE
        }
        else {
          selection_base[[i]] = FALSE
        }
      }
      if (all(unlist(selection_base))) {
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    })
    sum_matrix_qc = sum_matrix_qc[, raw_counts_allele_update]
    index = apply(sum_matrix_qc, 2, function(x) {
      logic = which(x > number_reads)
      return((logic))
    })
  }
  raw_counts_allele_filter_qc = raw_counts_allele_filter[row.names(sum_matrix_qc),
                                                         colnames(raw_counts_allele_filter) %in% colnames(sum_matrix_qc)]
  norm_counts_allele_filter_qc = matrix(0, ncol = length(colnames(raw_counts_allele_filter_qc)),
                                        nrow = length(row.names(raw_counts_allele_filter_qc)))
  colnames(norm_counts_allele_filter_qc) = name_position[colnames(raw_counts_allele_filter) %in%
                                                           colnames(sum_matrix_qc)]
  row.names(norm_counts_allele_filter_qc) = row.names(sum_matrix_qc)
  colnames(raw_counts_allele_filter_qc) = name_position[colnames(raw_counts_allele_filter) %in%
                                                          colnames(sum_matrix_qc)]
  for (i in 1:length(colnames(raw_counts_allele_filter_qc))) {
    norm_counts_allele_filter_qc[, i] = raw_counts_allele_filter_qc[,
                                                                    i]/(sum_matrix_qc[, colnames(sum_matrix_qc) %in%
                                                                                        colnames(raw_counts_allele_filter_qc)[i]])
  }
  norm_counts_allele_filter_qc[is.na(norm_counts_allele_filter_qc)] = 0
  sum_allele = as.list(rep(0, length(colnames(sum_matrix_qc))))
  if (filtering == 1) {
    for (i in 1:length(colnames(sum_matrix_qc))) {
      sum_allele[[i]] = apply(norm_counts_allele_filter_qc[,
                                                           colnames(norm_counts_allele_filter_qc) %in% colnames(sum_matrix_qc)[i]],
                              1, sum)
    }
  }
  if (filtering == 2) {
    for (i in 1:length(colnames(sum_matrix_qc))) {
      index_cell = index[[which(names(index) == colnames(sum_matrix_qc)[i])]]
      sum_allele[[i]] = apply(norm_counts_allele_filter_qc[index_cell,
                                                           colnames(norm_counts_allele_filter_qc) %in% colnames(sum_matrix_qc)[i]],
                              1, sum)
    }
  }
  for (i in 1:length(colnames(sum_matrix_qc))) {
    if (all(sum_allele[[i]] > 0.99) != TRUE) {
      print(paste("Warning:sum of allele frequencies is not one for one or more cells in base ",
                  i, sep = ""))
    }
  }
  heteroplasmy_matrix = matrix(0, ncol = length(colnames(sum_matrix_qc)),
                               nrow = length(row.names(norm_counts_allele_filter_qc)))
  colnames(heteroplasmy_matrix) = colnames(sum_matrix_qc)
  row.names(heteroplasmy_matrix) = row.names(norm_counts_allele_filter_qc)
  colnames(norm_counts_allele_filter_qc) = name_position[name_position %in%
                                                           colnames(raw_counts_allele_filter_qc)]
  for (i in colnames(sum_matrix_qc)) {
    heteroplasmy_matrix[, i] = apply(norm_counts_allele_filter_qc[,
                                                                  colnames(norm_counts_allele_filter_qc) %in% i], 1,
                                     function(x) {
                                       x = as.numeric(as.vector(x))
                                       if (sum(x) == 0) {
                                         return(-1)
                                       }
                                       else {
                                         return(1 - max(x))
                                       }
                                     })
  }
  colnames(norm_counts_allele_filter_qc) = name_position_allele[name_position %in%
                                                                  colnames(raw_counts_allele_filter_qc)]
  return(list(sum_matrix, sum_matrix_qc, heteroplasmy_matrix,
              norm_counts_allele_filter_qc, index))
}

















