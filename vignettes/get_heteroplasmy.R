
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

















