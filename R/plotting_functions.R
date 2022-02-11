#' plot_heteroplasmy
#' @param position Character name of the base to plot.
#' @param heteroplasmy_matrix Third element returned by \emph{get_heteroplasmy}.
#' @param cluster Vector specifying a partition of the samples.
#' @param index Fifth element returned by \emph{get_heteroplasmy}.
#' @return ggplot object of the heteroplasmy level of a specific base across
#' samples divided according to cluster.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_heteroplasmy
plot_heteroplasmy <- function(position, heteroplasmy_matrix, cluster, index) {
  position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
  if (is.null(index)) {
    heteroplasmy <- heteroplasmy_matrix[, position]
  }
  else {
    cluster <- cluster[as.numeric(index[[position]])]
    heteroplasmy <- heteroplasmy_matrix[as.numeric(index[[position]]), position]
  }
  dati_scatter_1 <- data.frame(as.character(cluster), heteroplasmy)
  ggplot2::ggplot(data = dati_scatter_1, ggplot2::aes(x = factor(cluster), y = heteroplasmy, col = factor(cluster))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") +
    ggplot2::labs(y = "Heteroplasmy") + ggplot2::ggtitle(position)
}













#' plot_allele_frequency
#' @param allele_matrix Fourth element returned by \emph{get_heteroplasmy}.
#' @param names_allele_qc Character vector with length equal to n_col of
#' \emph{allele_matrix}. Each element specifies the name of the base and the allele.
#' @param names_position_qc Character vector with length equal to n_col of
#' \emph{allele_matrix}. Each element specifies the name of the base.
#' @param size_text Character specifying the size of the text for \emph{gridExtra} function
#' \emph{grid.arrange})
#' @return \emph{grid.arrange} plot of allele frequencies of a specific base
#' across samples divided according to cluster.
#' @inheritParams plot_heteroplasmy
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://cran.r-project.org/package=gridExtra}
#' @export plot_allele_frequency
plot_allele_frequency <- function(position, heteroplasmy_matrix, allele_matrix, cluster,
                                  names_allele_qc, names_position_qc, size_text, index)
{
  position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
  if (is.null(index)) {
    a <- position
    colnames(allele_matrix) <- names_position_qc
    allele <- colnames(allele_matrix)[names_allele_qc == a]
    allele_1 <- allele_matrix[, allele[1]]
    data_plot_1 <- data.frame(cluster, allele_1)
    allele_2 <- allele_matrix[, allele[2]]
    data_plot_2 <- data.frame(cluster, allele_2)
    allele_3 <- allele_matrix[, allele[3]]
    data_plot_3 <- data.frame(cluster, allele_3)
    allele_4 <- allele_matrix[, allele[4]]
    data_plot_4 <- data.frame(cluster, allele_4)
  }
  else {
    cluster <- cluster[as.numeric(index[[position]])]
    a <- position
    colnames(allele_matrix) <- names_position_qc
    allele <- colnames(allele_matrix)[names_allele_qc == a]
    allele_1 <- allele_matrix[as.numeric(index[[position]]),
                             allele[1]]
    data_plot_1 <- data.frame(cluster, allele_1)
    allele_2 <- allele_matrix[as.numeric(index[[position]]), allele[2]]
    data_plot_2 <- data.frame(cluster, allele_2)
    allele_3 <- allele_matrix[as.numeric(index[[position]]), allele[3]]
    data_plot_3 <- data.frame(cluster, allele_3)
    allele_4 <- allele_matrix[as.numeric(index[[position]]), allele[4]]
    data_plot_4 <- data.frame(cluster, allele_4)
  }
  plot_1 <- ggplot2::ggplot(data = data_plot_1, ggplot2::aes(x = factor(cluster), y = allele_1, col = factor(cluster))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggplot2::ggtitle(allele[1]) + ggplot2::ylim(0, 1) + ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))



  plot_2 <- ggplot2::ggplot(data = data_plot_2, ggplot2::aes(x = factor(cluster), y = allele_2, col = factor(cluster))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele Frequency")+ ggplot2::ggtitle(allele[2]) + ggplot2::ylim(0, 1) + ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))
  plot_3 = ggplot2::ggplot(data = data_plot_3, ggplot2::aes(x = factor(cluster), y = allele_3, col = factor(cluster) ))+ ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") +
    ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggplot2::ggtitle(allele[3]) + ggplot2::ylim(0, 1)  + ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))

  plot_4 <- ggplot2::ggplot(data = data_plot_4, ggplot2::aes(x = factor(cluster), y = allele_4, col = factor(cluster,
  ))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") +
    ggplot2::ggtitle(allele[4]) + ggplot2::ylim(0, 1)  + ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))
  gridExtra::grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 2, nrow = 2)
}


















#' plot_dpt
#' @param time Vector of diffusion pseudo time,with length equal to n_row of
#' \emph{heteroplasmy_matrix}.
#' @param gam_fit_result Data frame returned by \emph{dpt_test}.
#' @return ggplot object of the heteroplasmy level of a specific base across
#' samples and the GAM fitted curve. The title shows the adjusted p value (FDR)
#' for the position obtained from \emph{get_heteroplasmy}.
#' @inheritParams plot_heteroplasmy
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://cran.r-project.org/package=gam}
#' @export plot_dpt
plot_dpt <- function(position, heteroplasmy_matrix, cluster, time, gam_fit_result, index)
{
  if (is.null(index)) {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    mutation <- heteroplasmy_matrix[, position]
    data_plot_1 <- data.frame(time, mutation)
  }
  else {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    cluster <- cluster[as.numeric(index[[position]])]
    mutation <- heteroplasmy_matrix[as.numeric(index[[position]]), position]
    time <- time[as.numeric(index[[position]])]
    data_plot_1 <- data.frame(time, mutation)
  }
  plot <- ggplot2::ggplot(data_plot_1, ggplot2::aes(x = time, y = mutation)) +
    ggplot2::geom_point(ggplot2::aes(colour = (cluster))) + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Losing Score") + ggplot2::labs(y = "Heteroplasmy") +
    ggplot2::ggtitle(paste(as.character(position), "adjusted p value=",
                    round(gam_fit_result[position, 1], (floor(-log10(gam_fit_result[position, 1]))) + 1), sep = " ")) +
    ggplot2::geom_smooth(method = "loess", se = F, color = "black", fullrange = F)
  plot
}








#' plot_batch_epiblast
#' @noRd

plot_batch_epiblast <- function(position, heteroplasmy_matrix, batch, cluster, text_size, index)
{
  if (is.null(index)) {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    heteroplasmy <- heteroplasmy_matrix[, position]
    data_plot_1 <- data.frame(batch, heteroplasmy)
  }
  else {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    batch <- batch[as.numeric(index[[position]])]
    heteroplasmy <- heteroplasmy_matrix[as.numeric(index[[position]]), position]
    cluster <- cluster[as.numeric(index[[position]])]
    data_plot_1 <- data.frame(batch, heteroplasmy)
  }
  ggplot2::ggplot(data = data_plot_1, ggplot2::aes(x = batch, y = heteroplasmy,
                  col = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"), ordered = TRUE))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Batch") +
    ggplot2::labs(y = "Heteroplasmy") + ggplot2::ggtitle(position) +
    ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"), labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast")) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = text_size))
}











#' plot_batch
#' @param batch Vector of batch names,with length equal to n_row of
#' \emph{heteroplasmy_matrix}.
#' @param text_size Character specifying the size of the text for ggplot2.
#' @inheritParams plot_heteroplasmy
#' @return \emph{ggplot2} object of the heteroplasmy level of a specific base across
#' samples divided according to batch.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_batch
plot_batch <- function(position, heteroplasmy_matrix, batch, cluster, text_size, index)
{
  if (is.null(index)) {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    heteroplasmy <- heteroplasmy_matrix[, position]
    data_plot_1 <- data.frame(batch, heteroplasmy)
  }
  else {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    batch <- batch[as.numeric(index[[position]])]
    heteroplasmy <- heteroplasmy_matrix[as.numeric(index[[position]]), position]
    cluster <- cluster[as.numeric(index[[position]])]
    data_plot_1 <- data.frame(batch, heteroplasmy)
  }
  ggplot2::ggplot(data = data_plot_1, ggplot2::aes(x = batch, y = heteroplasmy, col = factor(cluster))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Batch") + ggplot2::labs(y = "Heteroplasmy") +
    ggplot2::ggtitle(position) + ggplot2::theme(axis.text = ggplot2::element_text(size = text_size))

}








#' plot_distribution
#' @param quantity_counts_cell Vector returned by
#' \emph{get_distribution}
#' @param name_x Character name specifyng the xlab argument in \emph{ggplot2}.
#' @param name_title Character name specifyng the ggtitle argument in \emph{ggplot2}.
#' @return \emph{ggplot2} density plot of the Vector quantity_counts_cell.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_distribution
plot_distribution <- function(quantity_counts_cell, name_x, name_title)
{
  data_plot <- data.frame(quantity_counts_cell)
  ggplot2::ggplot(data_plot, ggplot2::aes(x = quantity_counts_cell)) +
    ggplot2::geom_density(color = "darkblue", fill = "lightblue") +
    ggplot2::xlab(name_x) + ggplot2::ggtitle(name_title)
}








#' plot_condition
#' @param distribution_1,distribution_2 Numeric vector
#' @param label_1 Character vector of length equal to distribution_1
#' @param label_2 Character vector of length equal to distribution_2
#' @param name_y Character name specifyng the ylab argument in ggplot2.
#' @inheritParams plot_distribution
#' @return \emph{ggplot2} boxplot of the quantities specified by \emph{distribution_1}
#' and \emph{distribution_2}, separated by the conditions denoted by
#' \emph{label_1} and \emph{label_2}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_condition
plot_condition <- function(distribution_1, distribution_2, label_1, label_2, name_x, name_y, name_title)
{
  distribution_all <- c(distribution_1, distribution_2)
  label_all <- c(label_1, label_2)
  data_plot <- data.frame(distribution_all, label_all)
  res <- wilcox.test(distribution_all ~ label_all, data = data_plot, exact = FALSE)
  add_test = paste("Wilcoxon test:", round(res$p.value, 4), sep = " ")
  plot_gene <- ggplot2::ggplot(data_plot, ggplot2::aes(x = label_all, y = distribution_all, col = label_all)) +
    ggplot2::geom_boxplot() + ggplot2::xlab(name_x) + ggplot2::ylab(name_y) +
    ggplot2::labs(col = name_x) + ggplot2::ggtitle(paste(name_title, add_test, sep = " - "))
  return(plot_gene)
}






#' plot_heteroplasmy_epiblast
#' @noRd

plot_heteroplasmy_epiblast <- function(position, heteroplasmy_matrix, cluster, index)
{
  position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
  if (is.null(index)) {
    heteroplasmy <- heteroplasmy_matrix[, position]
  }
  else {
    cluster <- cluster[as.numeric(index[[position]])]
    heteroplasmy <- heteroplasmy_matrix[as.numeric(index[[position]]), position]
  }
  dati_scatter_1 <- data.frame(as.character(cluster), heteroplasmy)
  ggplot2::ggplot(data = dati_scatter_1, ggplot2::aes(x = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"), ordered = TRUE),
                                                      y = heteroplasmy,
                                                      col = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"), ordered = TRUE))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Heteroplasmy") +
    ggplot2::ggtitle(position) + ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"), labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"))
}








#' plot_allele_frequency_epiblast
#' @noRd
plot_allele_frequency_epiblast <- function(position, heteroplasmy_matrix, allele_matrix, cluster, names_allele_qc, names_position_qc, size_text, index) {
  position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
  if (is.null(index)) {
    a <- position
    colnames(allele_matrix) <- names_position_qc
    allele <- colnames(allele_matrix)[names_allele_qc == a]
    allele_1 <- allele_matrix[, allele[1]]
    data_plot_1 <- data.frame(cluster, allele_1)
    allele_2 <- allele_matrix[, allele[2]]
    data_plot_2 <- data.frame(cluster, allele_2)
    allele_3 <- allele_matrix[, allele[3]]
    data_plot_3 <- data.frame(cluster, allele_3)
    allele_4 <- allele_matrix[, allele[4]]
    data_plot_4 <- data.frame(cluster, allele_4)
  }
  else {
    cluster <- cluster[as.numeric(index[[position]])]
    a <- position
    colnames(allele_matrix) <- names_position_qc
    allele <- colnames(allele_matrix)[names_allele_qc == a]
    allele_1 <- allele_matrix[as.numeric(index[[position]]), allele[1]]
    data_plot_1 <- data.frame(cluster, allele_1)
    allele_2 <- allele_matrix[as.numeric(index[[position]]), allele[2]]
    data_plot_2 <- data.frame(cluster, allele_2)
    allele_3 <- allele_matrix[as.numeric(index[[position]]), allele[3]]
    data_plot_3 <- data.frame(cluster, allele_3)
    allele_4 <- allele_matrix[as.numeric(index[[position]]), allele[4]]
    data_plot_4 <- data.frame(cluster, allele_4)
  }
  plot_1 <- ggplot2::ggplot(data = data_plot_1, ggplot2::aes(x = factor(cluster,
                                                            levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE), y = allele_1, col = factor(cluster,
                                                            levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"), ordered = TRUE))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggplot2::ggtitle(allele[1]) + ggplot2::ylim(0, 1) + ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"), labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast")) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))



  plot_2 <- ggplot2::ggplot(data = data_plot_2, ggplot2::aes(x = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE), y = allele_2,
                                                            col = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) + ggplot2::geom_point() + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele Frequency")+ ggplot2::ggtitle(allele[2]) + ggplot2::ylim(0, 1) +
    ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"), labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast")) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))

  plot_3 <- ggplot2::ggplot(data = data_plot_3, ggplot2::aes(x = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE), y = allele_3,
                                                            col = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE))) +
    ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") +
    ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggplot2::ggtitle(allele[3]) + ggplot2::ylim(0, 1) +
    ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"), labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast")) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))

  plot_4 <- ggplot2::ggplot(data = data_plot_4, ggplot2::aes(x = factor(cluster,
                                                            levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE), y = allele_4, col = factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                            ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") +
    ggplot2::ggtitle(allele[4]) + ggplot2::ylim(0, 1) +
    ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"), labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast")) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = size_text))
  gridExtra::grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 2, nrow = 2)
}






#' plot_dpt_epiblast
#' @noRd
plot_dpt_epiblast <- function(position, heteroplasmy_matrix, cluster, time, gam_fit_result, index)
{
  if (is.null(index)) {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    cluster <- factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"), ordered = TRUE)
    mutation <- heteroplasmy_matrix[, position]
    data_plot_1 <- data.frame(time, mutation)
  }
  else {
    position <- colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) == position]
    cluster <- cluster[as.numeric(index[[position]])]
    cluster <- factor(cluster, levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"), ordered = TRUE)
    mutation <- heteroplasmy_matrix[as.numeric(index[[position]]), position]
    time <- time[as.numeric(index[[position]])]
    data_plot_1 <- data.frame(time, mutation)
  }
  plot <- ggplot2::ggplot(data_plot_1, ggplot2::aes(x = time, y = mutation)) +
    ggplot2::geom_point(ggplot2::aes(colour = (cluster))) + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Losing Score") + ggplot2::labs(y = "Heteroplasmy") +
    ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"),
                                labels = c("Winning Epiblast", "Intermediate", "Losing Epiblast")) +
    ggplot2::ggtitle(paste(as.character(position), "adjusted p value=",
                           round(gam_fit_result[position, 1], (floor(-log10(gam_fit_result[position, 1]))) + 1), sep = " ")) +
    ggplot2::geom_smooth(method = "loess", se = F, color = "black", fullrange = F)
  plot
}












#' plot_genome_coverage
#' @param biomart_file Character string with full path to the txt file
#' downloaded from BioMart \url{https://m.ensembl.org/info/data/biomart/index.html} . It must have
#' the following five columns:Gene.stable.ID, Gene.name, Gene.start..bp.,
#' Gene.end..bp., Chromosome.scaffold.name
#' @param path_fasta Character string with full path to the fasta file of the
#' genomic region of interest. It should be the same file used in
#' \emph{get_raw_counts_allele}.
#' @param chr_name Character specifyng the name of the chromosome of interest.
#' It must be one of the names in the \emph{Chromosome.scaffold.name} column
#' from the \emph{biomart_file}.
#' @inheritParams plot_heteroplasmy
#' @inheritParams get_raw_counts_allele
#' @return Plot as returned by \emph{karyoploteR} package.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{http://bioconductor.org/packages/release/bioc/html/karyoploteR.html}
#' @export plot_genome_coverage
plot_genome_coverage = function(biomart_file, path_fasta, chr_name, heteroplasmy_matrix) {
  if (!(requireNamespace("karyoploteR", quietly = TRUE) & requireNamespace("regioneR", quietly = TRUE) )) {
    stop("Package karyoploteR and regioneR needed for this function to work. Please install them: BiocManager::install('karyoploteR') and BiocManager::install('regioneR') ")
  }
  else{
  base_name <- rep(0, length(colnames(heteroplasmy_matrix)))
  for ( i in 1:length(colnames(heteroplasmy_matrix))) {
    base_name[i] <- strsplit(colnames(heteroplasmy_matrix), "_")[[i]][1]
  }
  colnames(heteroplasmy_matrix) <- base_name
  seq_name <- names(path_fasta)
  sequence <- paste(path_fasta)
  df_bg <- data.frame(seq_name, sequence)
  name_sequences <- rep(0, length(df_bg$seq_name))
  for (i in 1:length(df_bg$seq_name)) {
    name_sequences[i] <- strsplit(df_bg$seq_name[i], split = " ")[[1]][1]
  }

  df_bg$seq_name <- name_sequences

  length_sequence <- rep(0, length(df_bg$seq_name))
  for (q in 1:length(df_bg$seq_name)) {
    length_sequence[q] <- length(unlist(strsplit(df_bg$sequence[q], split = NULL)))

  }

  length_chr <- length_sequence[df_bg$seq_name == "MT"]

  mt_name <- biomart_file
  mt_name <- mt_name[mt_name$Chromosome.scaffold.name == chr_name, ]
  all_name <- mt_name$Gene.name
  p <- rep(0, length(all_name))
  for (j in 1:length(all_name)) {
    i <- which(mt_name$Gene.name == all_name[j])
    all_name_bases <- seq(mt_name$Gene.start..bp.[i], mt_name$Gene.end..bp.[i])
    all_name_bases <- as.character(all_name_bases)

    if (sum(colnames(heteroplasmy_matrix) %in% all_name_bases)>0) {
      p[i] <- all_name[i]
    }

  }
  mt_name_small <- mt_name[mt_name$Gene.name %in% p, ]

  start_plot <- c(mt_name_small$Gene.start..bp., mt_name_small$Gene.end..bp.)
  all_name_plot <- paste0(mt_name_small$Gene.name, "-start")
  all_name_plot_2 <- paste0(mt_name_small$Gene.name, "-end")
  names(start_plot) <- c(all_name_plot, all_name_plot_2)
  start_plot <- sort(start_plot)


  mt_name_new <- mt_name[!(mt_name$Gene.stable.ID %in% mt_name_small$Gene.stable.ID), ]

  start_plot_new <- c(mt_name_new$Gene.start..bp., mt_name_new$Gene.end..bp.)
  all_name_plot <- paste0(mt_name_new$Gene.name, "-start")
  all_name_plot_2 <- paste0(mt_name_new$Gene.name, "-end")
  names(start_plot_new) <- c(all_name_plot, all_name_plot_2)
  start_plot_new <- sort(start_plot_new)


  colnames(heteroplasmy_matrix) <- base_name
  all_base <- seq(1, length_chr)
  all_base <- as.character(all_base)

  name_all_base <- rep(0, length_chr)
  name_all_base[all_base %in% colnames(heteroplasmy_matrix)] <- "Covered"
  name_all_base[!all_base %in% colnames(heteroplasmy_matrix)] <- "Not_Covered"
  name_all_base <- paste0(name_all_base, all_base)
  colore <- rep(0, length_chr)
  colore[all_base %in% colnames(heteroplasmy_matrix)] <- "red"
  colore[!all_base %in% colnames(heteroplasmy_matrix)] <- "black"


  custom.genome <- regioneR::toGRanges(data.frame(chr = c(chr_name), start = c(1), end = c(length_chr)))

  kp <- karyoploteR::plotKaryotype(genome = custom.genome, plot.type = 2)

  df <- data.frame(chr = c(rep(chr_name, length(all_base))), start = all_base, end = all_base)
  karyoploteR::kpPlotRegions(kp, data = df, col = colore, border = colore, r0 = 0, r1 = 0.10)



  markers <- data.frame(chr = rep(chr_name, length(start_plot)), pos = start_plot, labels = names(start_plot))
  markers_new <- data.frame(chr = rep(chr_name, length(start_plot_new)), pos = start_plot_new, labels = names(start_plot_new))
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 1000, add.units = T, tick.len = 30)

  karyoploteR::kpPlotMarkers(kp, chr = markers$chr, x = markers$pos, labels = markers$labels, r0 = 0.10, r1 = 0.55, cex = 0.3, label.color = rep("red", length(markers_new)), label.margin = 1, line.color = rep("red", length(markers_new)))
  karyoploteR::kpPlotMarkers(kp, chr = markers_new$chr, x = markers_new$pos, labels = markers_new$labels, r0 = 0.10, r1 = 0.55, cex = 0.3, label.dist = 0.000001, data.panel = 2, marker.parts = c(0, 0.9, 0.1), label.margin = 0.0001)
  legend(x = "topleft", fill = c("red", "black"), legend = c("Covered", "Not covered"), cex = 0.5)

  }
}






#' plot_cells_coverage_epiblast
#' @noRd
plot_cells_coverage_epiblast = function(sum_matrix, cells_selected, cluster, interactive = FALSE) {
  cells_sum <- apply(sum_matrix[cells_selected, ], 1, sum)
  cells_sum <- cells_sum/1000
  cell_id <- names(cells_sum)
  df <- data.frame(cell_id, cells_sum, cluster)
  row.names(df) <- cell_id
  df <- df[order(cells_sum, decreasing = T), ]
  cell_id_plot <- factor(row.names(df), levels = row.names(df))
  cluster_plot_plot <- factor(df[, 3], levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"))
  if (interactive == FALSE) {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = cell_id_plot, y = df[, 2], fill = cluster_plot_plot)) +
      ggplot2::geom_bar(stat = "identity", ) +
      ggplot2::theme_minimal() + ggplot2::ylab("MT coverage per thousand reads") + ggplot2::xlab("Cells") +
      ggplot2::ggtitle("MT coverage per cell") + ggplot2::scale_fill_manual(values = c("#DD6400", "#0000FF", "#006400"),
                                                                          labels = c("Winner Epiblast", "Intermediate", "Loser Epiblast")) + ggplot2::labs(fill = "Category") +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + ggplot2::xlab("Cells")
    return(list(p))
  }
  if (interactive == TRUE) {
    if (! requireNamespace("plotly", quietly = TRUE)) {
      stop("Package plotly needed for interactive == TRUE. Please install it: install.packages('plotly') or set interactive == FALSE")
    }
    color_plot <- df[, 3]
    color_plot[df[, 3] == "Winner Epiblast"] <- "#DD6400"
    color_plot[df[, 3] == "Intermediate"] <- "#0000FF"
    color_plot[df[, 3] == "Loser Epiblast"] <- "#006400"




    fig <- plotly::plot_ly(
      x = (cell_id_plot),
      y = df[,2],
      name = "Cluster",
      type = "bar",
      marker = list(color = color_plot))
    fig <- fig %>% plotly::layout(title = "MT coverage per cell",
                                  xaxis = list(title = "Cells"),
                                  yaxis = list(title = "MT coverage per thousand reads"))

    return(list(fig))}

}










#' plot_cells_coverage
#' @param sum_matrix First element returned by the function
#' \emph{get_heteroplasmy}.
#' @param cells_selected Character vector of cells for which the coverage is
#' computed.
#' @param cluster Character vector with partition information for cells
#' specified in \emph{cells_selected}
#' @param interactive Logical. If TRUE an interactive plot is produced.
#' @return ggplot2 object (if \emph{interactive}=FALSE) or plotly object (if
#' \emph{interactive}=TRUE)
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://plotly.com/r/}
#' @export plot_cells_coverage
plot_cells_coverage = function(sum_matrix, cells_selected, cluster, interactive = FALSE) {
  cells_sum <- apply(sum_matrix[cells_selected, ], 1, sum)
  cells_sum <- cells_sum/1000
  cell_id <- names(cells_sum)
  df <- data.frame(cell_id, cells_sum, cluster)
  row.names(df) <- cell_id
  df <- df[order(cells_sum, decreasing = T), ]
  cell_id_plot <- factor(row.names(df), levels = row.names(df))
  if (interactive == FALSE) {
    p<-ggplot2::ggplot(data = df, ggplot2::aes(x = cell_id_plot, y = df[, 2], fill = df[, 3])) +
      ggplot2::geom_bar(stat = "identity", ) +
      ggplot2::theme_minimal() + ggplot2::ylab("MT coverage per thousand reads") + ggplot2::xlab("Cells") + ggplot2::ggtitle("MT coverage per cell") + ggplot2::labs(fill = "Category") +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + ggplot2::xlab("Cells")
    return(list(p))
  }
  if (interactive == TRUE) {

    if (! requireNamespace("plotly", quietly = TRUE)) {
      stop("Package plotly needed for interactive == TRUE. Please install it: install.packages('plotly') or set interactive == FALSE")
    }

    type <- length(levels(as.factor(df[, 3])))
    levels_type <- levels(as.factor(df[, 3]))
    color_plot <- rep(0, length(df[, 3]))
    for (i in 1:type) {
      color_plot[df[, 3] == levels_type[i]] <- gg_color_hue(type)[i]
    }
    fig <- plotly::plot_ly(
      x = (cell_id_plot),
      y = df[, 2],
      name = "Cluster",
      type = "bar",
      marker = list(color = color_plot ))
    fig <- fig %>% plotly::layout(title = "MT coverage per cell",
                                  xaxis = list(title = "Cells"),
                                  yaxis = list(title = "MT coverage per thousand reads"))

    return(list(fig))}


}







#' plot_base_coverage
#' @param sum_matrix_qc Second element returned by the function
#' \emph{get_heteroplasmy}.
#' @param selected_cells Character vector with cells used fro plotting the
#' coverage.
#' @param interactive Logical. If TRUE an interactive plot is produced.
#' @return \emph{ggplot2} object (if \emph{interactive}=FALSE) or plotly object (if
#' (if \emph{interactive}=TRUE)
#' @inheritParams plot_cells_coverage
#' @inheritParams plot_batch
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://plotly.com/r/}
#' @export plot_base_coverage
plot_base_coverage = function(sum_matrix, sum_matrix_qc, selected_cells, interactive = FALSE, text_size = 10) {
  condition_base <- rep("Removed", length(colnames(sum_matrix)))
  condition_base[colnames(sum_matrix) %in% colnames(sum_matrix_qc)] <- "Kept"
  if (length(selected_cells)>1) {
    base_mean <- apply(sum_matrix[selected_cells, ], 2, mean)}
  else{
    base_mean <- sum_matrix[selected_cells, ]

  }
  base_id <- names(base_mean)
  base_id <- substr(base_id, 0, nchar(base_id)-3)
  base_id <- as.numeric(base_id)
  df <- data.frame(base_id, base_mean, condition_base)
  row.names(df) <- base_id




  if (interactive == FALSE) {
    text_size <- 10

    p<-ggplot2::ggplot(data = df, ggplot2::aes(x = (df[, 1]), y = df[, 2], fill = df[, 3])) +
      ggplot2::geom_bar(stat = "identity", ) +
      ggplot2::theme_minimal() + ggplot2::ylab("Mean reads per base") + ggplot2::ggtitle("MT coverage per base") +
      ggplot2::labs(fill = "Category") + ggplot2::xlab("MT bases") +
      ggplot2::scale_x_continuous(breaks = round(seq(min((df[, 1])), max((df[, 1])), by = 1000), 1)) +
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                     axis.title.x = ggplot2::element_blank())


    return(list(p))
  }


  if (interactive == TRUE) {
    if (! requireNamespace("plotly", quietly = TRUE)) {
      stop("Package plotly needed for interactive == TRUE. Please install it: install.packages('plotly') or set interactive == FALSE")
    }
    text_size <- 10

    type <- length(levels(as.factor(df[, 3])))
    levels_type <- levels(as.factor(df[, 3]))
    color_plot <- rep(0, length(df[, 3]))
    for (i in 1:type) {
      color_plot[df[, 3] == levels_type[i]] <- gg_color_hue(type)[i]
    }
    fig <- plotly::plot_ly(
      x  =  (df[, 1]),
      y  = df[,2],
      name = "Cluster",
      type = "bar",
      marker = list(color = color_plot ))

    fig <- fig %>% plotly::layout(title = "MT coverage per base",
                                  xaxis = list(title = "MT bases"),
                                  yaxis = list(title = "Mean reads per base"))


    return(list(fig))
  }
}









#' plot_correlation_bases
#' @param bases_vector Character vector specyfing the bases for which the
#' spearman correlation across samples is computed.
#' @inheritParams plot_heteroplasmy
#' @return Heatmap plot produced by function \emph{Heatmap}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#' @export plot_correlation_bases
plot_correlation_bases = function(bases_vector, index, heteroplasmy_matrix) {
  if (is.null(index)) {

    dat <- heteroplasmy_matrix[, bases_vector]

    common_idx <- seq(1, length(row.names(heteroplasmy_matrix)))
  }

  if (!is.null(index)) {
  common_idx <- rep(list(0), length(bases_vector))
  for (i in 1:length(bases_vector)) {

    common_idx[[i]] <- as.numeric(index[[bases_vector[i]]])
  }


  common_idx <- Reduce(intersect, common_idx)

  dat <- heteroplasmy_matrix[common_idx, bases_vector]}


  correlation <- cor(dat, method = "spearman")




  ht21 <- ComplexHeatmap::Heatmap(correlation
                                 ,cluster_rows = TRUE
                                 , col = circlize::colorRamp2(c(0, round(max(correlation))), c("white", "red"))
                                 , name = "Spearman correlation"
                                 , column_title = paste0("Spearman correlation with ", length(common_idx), " cells")
                                 , cluster_columns = TRUE

                                 , row_names_gp = gpar(fontsize = 8)
                                 ,column_names_gp = gpar(fontsize = 8)
                                 , show_column_names = T
                                 , show_row_names = T

  )

  ComplexHeatmap::draw(ht21)



}







#' plot_spider_chart
#' @param name_base Character name specyfing the base.
#' @inheritParams plot_heteroplasmy
#' @return radarchart plot produced by function \emph{radarchart}.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://rdrr.io/cran/fmsb/man/radarchart.html}
#' @export plot_spider_chart
plot_spider_chart = function(name_base, cluster, heteroplasmy_matrix, index) {

  if (! requireNamespace("fmsb", quietly = TRUE)) {
    stop("Package fmsb needed for this function. Please install it: install.packages('fmsb')")
  }

  levels_cluster <- levels(factor(cluster))
  score_base <- rep(1, length(levels_cluster))
  if (is.null(index)) {
    cluster_small <- cluster
    heteroplasmy_matrix_small <- heteroplasmy_matrix
  }
  if (!is.null(index)) {
  cluster_small <- cluster[as.numeric(index[[name_base]])]
  heteroplasmy_matrix_small <- heteroplasmy_matrix[as.numeric(index[[name_base]]), ]}
  for (i in 1:length(levels_cluster)) {
    score_base[i] <- mean(heteroplasmy_matrix_small[cluster_small == levels_cluster[i], name_base])
    score_base[i] <- round(score_base[i], 3)
  }

  score_base_max <- rep(0, length(levels_cluster))
  for (i in 1:length(levels_cluster)) {
    score_base_max[i] <- max(apply(heteroplasmy_matrix_small[cluster_small == levels_cluster[i], ], 2, mean))
  }

  score_base_max <- max(as.vector(score_base_max))


  dat <- data.frame(score_base)
  row.names(dat) <- levels_cluster
  colnames(dat) <- name_base
  dat <- t(dat)

  max_min <- data.frame(rep(0, length(colnames(dat))), rep(score_base_max, length(colnames(dat))))
  max_min <- t(max_min)
  rownames(max_min) <- c("Max", "Min")


  df <- rbind(max_min, dat)



  student1_data <- df[c("Min", "Max", name_base), ]


  build_radarchart <- function(data, color ,
                               vlabels = colnames(data), vlcex = 0.7,
                               caxislabels = NULL, title = NULL) {

    fmsb::radarchart(
      data, axistype = 0,

      pcol = color, pfcol = color, plwd = 2, plty = 1,

      cglcol = "grey", cglty = 1, cglwd = 0.8,

      axislabcol = "grey",

      vlcex = vlcex, vlabels = vlabels,
      caxislabels = caxislabels, title = title
    )
  }

  op <- par(mar = c(1, 2, 2, 1))

  student1_data <- df[c("Min", "Max", name_base), ]
  student1_data <- as.data.frame(student1_data)
  build_radarchart(student1_data, color = gg_color_hue(1)[1], caxislabels = c(0, 5, 10, 15, 20), title = name_base)
  par(op)


}






#' plot_coordinate_heteroplasmy
#' @param coordinate_dm Dataframe whit samples on the rows and coordinates
#' names on the columns.
#' @inheritParams plot_heteroplasmy
#' @inheritParams plot_spider_chart
#' @return \emph{ggplot2} object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_coordinate_heteroplasmy
plot_coordinate_heteroplasmy = function(coordinate_dm, heteroplasmy_matrix, index, name_base) {
  if (is.null(index)) {
    base_col <- heteroplasmy_matrix[, name_base]
    row.names(coordinate_dm) <- row.names(heteroplasmy_matrix)
    dat <- cbind(coordinate_dm, base_col)
    dat <- dat[order(dat[, 3], decreasing = F), ]


    ggplot2::ggplot(dat, ggplot2::aes(coordinate_dm[, 1],  coordinate_dm[, 1])) +
      ggplot2::geom_point(ggplot2::aes(colour = dat[, 3] )) +
      ggplot2::scale_colour_gradient(low = "white", high = "midnight blue") +
      ggplot2::labs(col = "Heteroplasmy") + ggplot2::ggtitle(paste0("Heteroplasmy level of ", name_base))

  }

   if (!is.null(index)) {
  base_col <- heteroplasmy_matrix[as.numeric(index[[name_base]]), name_base]
  coordinate_dm <- coordinate_dm[as.numeric(index[[name_base]]), ]
  dat <- cbind(coordinate_dm, base_col)



  ggplot2::ggplot(dat, ggplot2::aes(coordinate_dm[, 1], coordinate_dm[, 2])) +
    ggplot2::geom_point(ggplot2::aes(colour = dat[, 3] )) +
    ggplot2::scale_colour_gradient(low = "white", high = "midnight blue") +
    ggplot2::labs(col = "Heteroplasmy") + ggplot2::ggtitle(paste0("Heteroplasmy level of ", name_base))
  }
}







#' plot_coordinate_cluster_epiblast
#' @noRd
plot_coordinate_cluster_epiblast = function(coordinate_dm, cluster) {
  Cluster_col <- cluster
  Cluster <- as.factor(Cluster_col)
  dat <- cbind(coordinate_dm, Cluster)

  cord_clu <- ggplot2::ggplot(dat, ggplot2::aes(x = coordinate_dm[, 1], y = coordinate_dm[, 2], color =  Cluster)) + ggplot2::geom_point()


  cord_clu + ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"),
                                      labels = c("Normal (winner) Epiblast", "Intermediate",  "Loser Epiblast"))+
    ggplot2::labs(x = "DC1", y = "DC2") + ggplot2::ggtitle("Diffusion Map")}







#' plot_coordinate_cluster
#' @inheritParams plot_coordinate_heteroplasmy
#' @inheritParams plot_heteroplasmy
#' @return \emph{ggplot2} object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_coordinate_cluster
plot_coordinate_cluster = function(coordinate_dm, cluster) {
  Cluster_col <- cluster
  Cluster <- as.factor(Cluster_col)
  dat <- cbind(coordinate_dm,Cluster)

  cord_clu <- ggplot2::ggplot(dat, ggplot2::aes(x = coordinate_dm[, 1], y = coordinate_dm[, 2], color = Cluster)) + ggplot2::geom_point()


  cord_clu + ggplot2::labs(x = "DC1", y = "DC2") + ggplot2::ggtitle("Diffusion Map")}









#' plot_heteroplasmy_variability
#' @param threshold Numeric value.
#' @param frac Logical. If FALSE the absolute number of cells that have at
#' least one base with heteroplasmy above \emph{threshold} are shown separated by
#' \emph{cluster}. If TRUE, then the fraction of cells are shown.
#' @inheritParams plot_heteroplasmy
#' @return \emph{ggplot2} object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @export plot_heteroplasmy_variability
plot_heteroplasmy_variability = function(heteroplasmy_matrix, cluster, threshold = 0.1, frac = FALSE, index) {
  if (is.null(index)) {
  number_het <- apply(heteroplasmy_matrix, 1, function(x) {
    threshold <- threshold
    x <- x[x>threshold]
    return(length(x))
  })
  }

  if (!is.null(index)) {
    number_het <- rep(0, length(row.names(heteroplasmy_matrix)))

    for ( i in 1:length(row.names(heteroplasmy_matrix))) {
      base_covered <- rep(0, length(colnames(heteroplasmy_matrix)))
      for (j in 1:length(colnames(heteroplasmy_matrix))) {
      if (sum(index[[colnames(heteroplasmy_matrix)[j]]] %in% i)>0) {
        base_covered[j] <- names(index)[j]
      }}
      base_covered <- base_covered[base_covered!=0]
      x <- heteroplasmy_matrix[i, base_covered]
      x <- x[x>threshold]
      number_het[i] <- length(x)
    }
  }

   if (frac == FALSE) {
    stages <- names(table(cluster[number_het>0]))
    number <- as.vector((table(cluster[number_het>0])))
    df <- data.frame(stage = stages, number = number)
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = stages, y = number)) +
      ggplot2::geom_bar(stat = "identity", color = "blue", fill = "blue") + ggplot2::ggtitle(paste("Number of cells with H>", threshold, " across stages", sep = "")) + ggplot2::ylab("Number of cells") + ggplot2::xlab("Stage")

  }

   if (frac == TRUE) {
    stages <- names(table(cluster[number_het>0]))
    number <- as.vector((table(cluster[number_het>0])))/as.vector((table(cluster[(cluster) %in% stages])))
    df <- data.frame(stage = stages,
                     number = number)
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = stages, y = number)) +
      ggplot2::geom_bar(stat = "identity", color = "blue", fill = "blue") +
      ggplot2::ggtitle(paste("Fraction of cells with H>", threshold, " across stages", sep = "")) +
      ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Stage")

  }
  return(list(p))}




