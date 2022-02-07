#' get_heteroplasmy
#'
#' It is one of the two main functions of the \strong{MitoHEAR} package
#' (together with \emph{get_raw_counts_allele}). It computes the allele
#' frequencies and the heteroplasmy matrix starting from the counts matrix
#' obtained with \emph{get_raw_counts_allele}.
#'
#' Starting from \emph{raw counts allele matrix}, the function performed two
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
#' must be covered by the sample (with counts>\emph{number_reads}), in order to keep
#' the sample for down-stream analysis.
#' @param filtering Numeric value equal to 1 or 2. If 1 then only the bases
#' that are covered by all the samples are kept for the downstream analysis. If
#' 2 then all the bases that are covered by more than 50\% of the the samples
#' in each cluster (specified by \emph{my.clusters}) are kept for the down-stream
#' analysis. Default is 1.
#' @param my.clusters Charachter vector specifying a partition of the samples.
#' It is only used when filtering is equal to 2. Deafult is NULL
#' @return It returns a list with 5 elements:
#'
#' \item{sum_matrix}{A matrix (n_row=number of sample, n_col=number of bases)
#' with the counts for each sample/base, for all the initial samples and bases
#' included in the \emph{raw counts allele matrix}.}
#'
#' \item{sum_matrix_qc}{A matrix (n_row=number of sample, n_col=number of
#' bases) with the counts for each sample/base, for all the samples and bases
#' that pass the two consequentially filtering steps.}
#'
#' \item{heteroplasmy_matrix}{A matrix with the same dimension of \emph{sum_matrix_qc}
#' where each entry (i,j) is the heteroplasmy for sample i in base j.}
#'
#' \item{allele_matrix}{A matrix (n_row=number of sample,
#' n_col=4*number of bases) with allele frequencies, for all the samples and
#' bases that pass the two consequentially filtering steps.}
#'
#' \item{index}{Indices of the samples that
#' cover a base, for all bases and samples that pass the two consequentially
#' filtering steps; if all the samples cover all the bases, then \emph{index} is NULL }
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
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
  if (length(unique(name_position))<2|length(row.names(raw_counts_allele))<2){stop(paste0('There are not at least 2 samples and at least 2 bases as input'))}

  sum_matrix = matrix(0, ncol = length(unique(name_position)), nrow = length(row.names(raw_counts_allele)))
  colnames(sum_matrix) = unique(name_position)
  row.names(sum_matrix) = row.names(raw_counts_allele)
  colnames(raw_counts_allele) = name_position
  for (i in unique(name_position)) {
    sum_matrix[, i] = apply(raw_counts_allele[, colnames(raw_counts_allele) %in%
                                                       i], 1, function(x) {
                                                         return(sum(x))
                                                       })
  }
  sum_matrix = as.matrix(sum_matrix)
  pos_cov = apply(sum_matrix, 1, function(x) {
    sum(x > number_reads)
  })
  if (sum(pos_cov > number_positions)<2){stop(paste0('There are not at least 2 samples that cover more at least 2 bases with more than ',number_reads," reads"))}

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
    if (sum(raw_counts_allele_update)<2){stop(paste0('There are no at least 2 bases that are covered by all samples with more than ',number_reads,' reads'))}

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
    if (sum(raw_counts_allele_update)<2){stop(paste0('There are not at least 2 bases that are covered by 50% of samples with more than ',number_reads,' reads'))}

    sum_matrix_qc = sum_matrix_qc[, raw_counts_allele_update]
    if (all(raw_counts_allele_update)){
      filtering=1
      index=NULL
    }
    else{
      index = apply(sum_matrix_qc, 2, function(x) {
        logic = which(x > number_reads)
        return((logic))
      })
    }

  }
  raw_counts_allele_qc = raw_counts_allele[row.names(sum_matrix_qc),
                                                         colnames(raw_counts_allele) %in% colnames(sum_matrix_qc)]
  allele_matrix = matrix(0, ncol = length(colnames(raw_counts_allele_qc)),
                                        nrow = length(row.names(raw_counts_allele_qc)))
  colnames(allele_matrix) = name_position[colnames(raw_counts_allele) %in%
                                                           colnames(sum_matrix_qc)]
  row.names(allele_matrix) = row.names(sum_matrix_qc)
  colnames(raw_counts_allele_qc) = name_position[colnames(raw_counts_allele) %in%
                                                          colnames(sum_matrix_qc)]
  for (i in 1:length(colnames(raw_counts_allele_qc))) {
    allele_matrix[, i] = raw_counts_allele_qc[,
                                                                    i]/(sum_matrix_qc[, colnames(sum_matrix_qc) %in%
                                                                                        colnames(raw_counts_allele_qc)[i]])
  }
  allele_matrix[is.na(allele_matrix)] = 0
  sum_allele = as.list(rep(0, length(colnames(sum_matrix_qc))))
  if (filtering == 1) {
    for (i in 1:length(colnames(sum_matrix_qc))) {
      sum_allele[[i]] = apply(allele_matrix[,
                                                           colnames(allele_matrix) %in% colnames(sum_matrix_qc)[i]],
                              1, sum)
    }
  }
  if (filtering == 2) {
    for (i in 1:length(colnames(sum_matrix_qc))) {
      index_cell = index[[which(names(index) == colnames(sum_matrix_qc)[i])]]
      sum_allele[[i]] = apply(allele_matrix[index_cell,
                                                           colnames(allele_matrix) %in% colnames(sum_matrix_qc)[i]],
                              1, sum)
    }
  }
  for (i in 1:length(colnames(sum_matrix_qc))) {
    if (all(sum_allele[[i]]==1) != TRUE) {

    warning(paste0('Sum of allele frequencies is not one for ',sum(sum_allele[[i]]!=1) ,' cells in base ',i))
    }
  }
  heteroplasmy_matrix = matrix(0, ncol = length(colnames(sum_matrix_qc)),
                               nrow = length(row.names(allele_matrix)))
  colnames(heteroplasmy_matrix) = colnames(sum_matrix_qc)
  row.names(heteroplasmy_matrix) = row.names(allele_matrix)
  colnames(allele_matrix) = name_position[name_position %in%
                                                           colnames(raw_counts_allele_qc)]
  for (i in colnames(sum_matrix_qc)) {
    heteroplasmy_matrix[, i] = apply(allele_matrix[,
                                                                  colnames(allele_matrix) %in% i], 1,
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
  colnames(allele_matrix) = name_position_allele[name_position %in%
                                                                  colnames(raw_counts_allele_qc)]
  return(list(sum_matrix, sum_matrix_qc, heteroplasmy_matrix,
              allele_matrix, index))
}

















