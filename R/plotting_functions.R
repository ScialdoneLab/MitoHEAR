

plot_heteroplasmy=function (position, heteroplasmy_matrix, cluster, index = NULL)
{
  position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                             position]
  if (is.null(index)) {
    heteroplasmy = heteroplasmy_matrix[, position]
  }
  else {
    cluster = cluster[as.numeric(index[[position]])]
    heteroplasmy = heteroplasmy_matrix[as.numeric(index[[position]]),
                                       position]
  }
  dati_scatter_1 = data.frame(as.character(cluster), heteroplasmy)
  ggplot2::ggplot(data = dati_scatter_1, aes(x = factor(cluster), y = heteroplasmy,
                                             col = factor(cluster))) + ggplot2::geom_boxplot(width = 0.9,
                                                                                             alpha = 0.8) + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") +
    labs(y = "Heteroplasmy") + ggplot2::ggtitle(position)
}









plot_allele_frequency=function (position, heteroplasmy_matrix, allele_matrix, cluster,
                                names_allele_qc, names_position_qc, size_text, index = NULL)
{
  position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                             position]
  if (is.null(index)) {
    a = position
    colnames(allele_matrix) = names_position_qc
    allele = colnames(allele_matrix)[names_allele_qc == a]
    allele_1 = allele_matrix[, allele[1]]
    data_plot_1 = data.frame(cluster, allele_1)
    allele_2 = allele_matrix[, allele[2]]
    data_plot_2 = data.frame(cluster, allele_2)
    allele_3 = allele_matrix[, allele[3]]
    data_plot_3 = data.frame(cluster, allele_3)
    allele_4 = allele_matrix[, allele[4]]
    data_plot_4 = data.frame(cluster, allele_4)
  }
  else {
    cluster = cluster[as.numeric(index[[position]])]
    a = position
    colnames(allele_matrix) = names_position_qc
    allele = colnames(allele_matrix)[names_allele_qc == a]
    allele_1 = allele_matrix[as.numeric(index[[position]]),
                             allele[1]]
    data_plot_1 = data.frame(cluster, allele_1)
    allele_2 = allele_matrix[as.numeric(index[[position]]),
                             allele[2]]
    data_plot_2 = data.frame(cluster, allele_2)
    allele_3 = allele_matrix[as.numeric(index[[position]]),
                             allele[3]]
    data_plot_3 = data.frame(cluster, allele_3)
    allele_4 = allele_matrix[as.numeric(index[[position]]),
                             allele[4]]
    data_plot_4 = data.frame(cluster, allele_4)
  }
  plot_1 = ggplot2::ggplot(data = data_plot_1, aes(x = factor(cluster), y = allele_1, col = factor(cluster)))+ ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggtitle(allele[1]) + ggplot2::ylim(0, 1) + ggplot2::theme(axis.text = element_text(size = size_text))



  plot_2 = ggplot2::ggplot(data = data_plot_2, aes(x = factor(cluster), y = allele_2, col = factor(cluster))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele Frequency")+ ggplot2::ggtitle(allele[2]) + ylim(0, 1) + ggplot2::theme(axis.text = element_text(size = size_text))
  plot_3 = ggplot2::ggplot(data = data_plot_3, aes(x = factor(cluster), y = allele_3, col = factor(cluster) ))+ ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") +
    ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggplot2::ggtitle(allele[3]) + ggplot2::ylim(0, 1)  + ggplot2::theme(axis.text = element_text(size = size_text))

  plot_4 = ggplot2::ggplot(data = data_plot_4, aes(x = factor(cluster), y = allele_4, col = factor(cluster,
  ))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") +
    ggplot2::ggtitle(allele[4]) + ggplot2::ylim(0, 1)  + ggplot2::theme(axis.text = element_text(size = size_text))
  gridExtra::grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 2, nrow = 2)
}














plot_dpt=function (position, heteroplasmy_matrix, cluster, time, gam_fit_result,
                            index = NULL)
{
  if (is.null(index)) {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    mutation = heteroplasmy_matrix[, position]
    data_plot_1 = data.frame(time, mutation)
  }
  else {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    cluster = cluster[as.numeric(index[[position]])]
    mutation = heteroplasmy_matrix[as.numeric(index[[position]]),
                                   position]
    time = time[as.numeric(index[[position]])]
    data_plot_1 = data.frame(time, mutation)
  }
  plot <- ggplot2::ggplot(data_plot_1, aes(x = time, y = mutation)) +
    ggplot2::geom_point(aes(colour = (cluster))) + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Losing Score") + ggplot2::labs(y = "Heteroplasmy") +
    ggplot2::ggtitle(paste(as.character(position), "adjusted p value=",
                           round(gam_fit_result[position, 1], (floor(-log10(gam_fit_result[position,
                                                                                           1]))) + 1), sep = " ")) + ggplot2::geom_smooth(method = "loess",
                                                                                                                                          se = F, color = "black", fullrange = F)
  plot
}




plot_batch_epiblast=function (position, heteroplasmy_matrix, batch, colour, cluster,
                              text_size, index = NULL)
{
  if (is.null(index)) {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    heteroplasmy = heteroplasmy_matrix[, position]
    data_plot_1 = data.frame(batch, heteroplasmy)
  }
  else {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    batch = batch[as.numeric(index[[position]])]
    heteroplasmy = heteroplasmy_matrix[as.numeric(index[[position]]),
                                       position]
    colour = colour[as.numeric(index[[position]])]
    cluster = cluster[as.numeric(index[[position]])]
    data_plot_1 = data.frame(batch, heteroplasmy)
  }
  ggplot2::ggplot(data = data_plot_1, aes(x = batch, y = heteroplasmy,
                                          col = factor(cluster, levels = c("Winner Epiblast", "Intermediate",
                                                                           "Loser Epiblast"), ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9,
                                                                                                                                        alpha = 0.8) + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Batch") +
    ggplot2::labs(y = "Heteroplasmy") + ggplot2::ggtitle(position) + ggplot2::scale_color_manual(values = c("#DD6400",
                                                                                                            "#0000FF", "#006400"), labels = c("Winner Epiblast",
                                                                                                                                              "Intermediate", "Loser Epiblast")) + ggplot2::theme(axis.text = element_text(size = text_size))
}







plot_batch=function (position, heteroplasmy_matrix, batch, colour, cluster,
                     text_size, index = NULL)
{
  if (is.null(index)) {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    heteroplasmy = heteroplasmy_matrix[, position]
    data_plot_1 = data.frame(batch, heteroplasmy)
  }
  else {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    batch = batch[as.numeric(index[[position]])]
    heteroplasmy = heteroplasmy_matrix[as.numeric(index[[position]]),
                                       position]
    colour = colour[as.numeric(index[[position]])]
    cluster = cluster[as.numeric(index[[position]])]
    data_plot_1 = data.frame(batch, heteroplasmy)
  }
  ggplot2::ggplot(data = data_plot_1, aes(x = batch, y = heteroplasmy,
                                          col = factor(cluster))) + ggplot2::geom_boxplot(width = 0.9,
                                                                                          alpha = 0.8) + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Batch") +
    ggplot2::labs(y = "Heteroplasmy") + ggplot2::ggtitle(position)+ ggplot2::theme(axis.text = element_text(size = text_size))

}




Plot_distribution=function (quantity_counts_cell, x_name, title_name)
{
  data_plot = data.frame(quantity_counts_cell)
  ggplot2::ggplot(data_plot, aes(x = quantity_counts_cell)) + ggplot2::geom_density(color = "darkblue",
                                                                                    fill = "lightblue") + xlab(x_name) + ggplot2::ggtitle(title_name)
}




Plot_boxplot=function (distribution_1, distribution_2, label_1, label_2, name_x,
                       name_y, name_title)
{
  distribution_all = c(distribution_1, distribution_2)
  label_all = c(label_1, label_2)
  data_plot = data.frame(distribution_all, label_all)
  res <- wilcox.test(distribution_all ~ label_all, data = data_plot,
                     exact = FALSE)
  add_test = paste("Wilcoxon test:", round(res$p.value, 4),
                   sep = " ")
  plot_gene = ggplot2::ggplot(data_plot, aes(x = label_all, y = distribution_all,
                                             col = label_all)) + ggplot2::geom_boxplot() + ggplot2::xlab(name_x) + ggplot2::ylab(name_y) +
    labs(col = name_x) + ggplot2::ggtitle(paste(name_title, add_test,
                                                sep = " - "))
  return(plot_gene)
}


plot_heteroplasmy_epiblast=function (position, heteroplasmy_matrix, cluster, index = NULL)
{
  position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                             position]
  if (is.null(index)) {
    heteroplasmy = heteroplasmy_matrix[, position]
  }
  else {
    cluster = cluster[as.numeric(index[[position]])]
    heteroplasmy = heteroplasmy_matrix[as.numeric(index[[position]]),
                                       position]
  }
  dati_scatter_1 = data.frame(as.character(cluster), heteroplasmy)
  ggplot2::ggplot(data = dati_scatter_1, aes(x = factor(cluster, levels = c("Winner Epiblast",
                                                                            "Intermediate", "Loser Epiblast"), ordered = TRUE), y = heteroplasmy,
                                             col = factor(cluster, levels = c("Winner Epiblast", "Intermediate",
                                                                              "Loser Epiblast"), ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9,
                                                                                                                                           alpha = 0.8) + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") +
    labs(y = "Heteroplasmy") + ggplot2::ggtitle(position) + ggplot2::scale_color_manual(values = c("#DD6400",
                                                                                                   "#0000FF", "#006400"), labels = c("Winner Epiblast",
                                                                                                                                     "Intermediate", "Loser Epiblast"))
}




plot_allele_frequency_epiblast=function (position, heteroplasmy_matrix, allele_matrix, cluster,
                                         names_allele_qc, names_position_qc, size_text, index = NULL)
{
  position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                             position]
  if (is.null(index)) {
    a = position
    colnames(allele_matrix) = names_position_qc
    allele = colnames(allele_matrix)[names_allele_qc == a]
    allele_1 = allele_matrix[, allele[1]]
    data_plot_1 = data.frame(cluster, allele_1)
    allele_2 = allele_matrix[, allele[2]]
    data_plot_2 = data.frame(cluster, allele_2)
    allele_3 = allele_matrix[, allele[3]]
    data_plot_3 = data.frame(cluster, allele_3)
    allele_4 = allele_matrix[, allele[4]]
    data_plot_4 = data.frame(cluster, allele_4)
  }
  else {
    cluster = cluster[as.numeric(index[[position]])]
    a = position
    colnames(allele_matrix) = names_position_qc
    allele = colnames(allele_matrix)[names_allele_qc == a]
    allele_1 = allele_matrix[as.numeric(index[[position]]),
                             allele[1]]
    data_plot_1 = data.frame(cluster, allele_1)
    allele_2 = allele_matrix[as.numeric(index[[position]]),
                             allele[2]]
    data_plot_2 = data.frame(cluster, allele_2)
    allele_3 = allele_matrix[as.numeric(index[[position]]),
                             allele[3]]
    data_plot_3 = data.frame(cluster, allele_3)
    allele_4 = allele_matrix[as.numeric(index[[position]]),
                             allele[4]]
    data_plot_4 = data.frame(cluster, allele_4)
  }
  plot_1 = ggplot2::ggplot(data = data_plot_1, aes(x = factor(cluster,
                                                              levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                              ordered = TRUE), y = allele_1, col = factor(cluster,
                                                                                                          levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                                                                          ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggtitle(allele[1]) + ggplot2::ylim(0, 1) + ggplot2::scale_color_manual(values = c("#DD6400",
                                                                                      "#0000FF", "#006400"), labels = c("Winner Epiblast",
                                                                                                                        "Intermediate", "Loser Epiblast")) + ggplot2::theme(axis.text = element_text(size = size_text))



  plot_2 = ggplot2::ggplot(data = data_plot_2, aes(x = factor(cluster,
                                                              levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                              ordered = TRUE), y = allele_2, col = factor(cluster,
                                                                                                          levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                                                                          ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele Frequency")+ ggplot2::ggtitle(allele[2]) + ylim(0, 1) + ggplot2::scale_color_manual(values = c("#DD6400",
                                                                                                                                                                                                                     "#0000FF", "#006400"), labels = c("Winner Epiblast",
                                                                                                                                                                                                                                                       "Intermediate", "Loser Epiblast")) + ggplot2::theme(axis.text = element_text(size = size_text))
  plot_3 = ggplot2::ggplot(data = data_plot_3, aes(x = factor(cluster,
                                                              levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                              ordered = TRUE), y = allele_3, col = factor(cluster,
                                                                                                          levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                                                                          ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") +
    ggplot2::labs(col = "Cluster") + ggplot2::labs(x = "Cluster") + ggplot2::labs(y = "Allele frequency") +
    ggplot2::ggtitle(allele[3]) + ggplot2::ylim(0, 1) + ggplot2::scale_color_manual(values = c("#DD6400",
                                                                                               "#0000FF", "#006400"), labels = c("Winner Epiblast",
                                                                                                                                 "Intermediate", "Loser Epiblast")) + ggplot2::theme(axis.text = element_text(size = size_text))

  plot_4 = ggplot2::ggplot(data = data_plot_4, aes(x = factor(cluster,
                                                              levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                              ordered = TRUE), y = allele_4, col = factor(cluster,
                                                                                                          levels = c("Winner Epiblast", "Intermediate", "Loser Epiblast"),
                                                                                                          ordered = TRUE))) + ggplot2::geom_boxplot(width = 0.9, alpha = 0.8) +
    ggplot2::geom_point() + ggplot2::labs(col = "Cluster") + ggplot2::labs(y = "Allele Frequency") + ggplot2::labs(x = "Cluster") +
    ggplot2::ggtitle(allele[4]) + ggplot2::ylim(0, 1) + ggplot2::scale_color_manual(values = c("#DD6400",
                                                                                               "#0000FF", "#006400"), labels = c("Winner Epiblast",
                                                                                                                                 "Intermediate", "Loser Epiblast")) + ggplot2::theme(axis.text = element_text(size = size_text))
  gridExtra::grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 2, nrow = 2)
}


plot_dpt_epiblast=function (position, heteroplasmy_matrix, cluster, time, gam_fit_result,
                                     index = NULL)
{
  if (is.null(index)) {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    cluster = factor(cluster, levels = c("Winner Epiblast",
                                         "Intermediate", "Loser Epiblast"), ordered = TRUE)
    mutation = heteroplasmy_matrix[, position]
    data_plot_1 = data.frame(time, mutation)
  }
  else {
    position = colnames(heteroplasmy_matrix)[colnames(heteroplasmy_matrix) ==
                                               position]
    cluster = cluster[as.numeric(index[[position]])]
    cluster = factor(cluster, levels = c("Winner Epiblast",
                                         "Intermediate", "Loser Epiblast"), ordered = TRUE)
    mutation = heteroplasmy_matrix[as.numeric(index[[position]]),
                                   position]
    time = time[as.numeric(index[[position]])]
    data_plot_1 = data.frame(time, mutation)
  }
  plot <- ggplot2::ggplot(data_plot_1, aes(x = time, y = mutation)) +
    ggplot2::geom_point(aes(colour = (cluster))) + ggplot2::labs(col = "Cluster") +
    ggplot2::labs(x = "Losing Score") + ggplot2::labs(y = "Heteroplasmy") +
    ggplot2::scale_color_manual(values = c("#DD6400", "#0000FF", "#006400"),
                                labels = c("Winning Epiblast", "Intermediate", "Losing Epiblast")) +
    ggplot2::ggtitle(paste(as.character(position), "adjusted p value=",
                           round(gam_fit_result[position, 1], (floor(-log10(gam_fit_result[position,
                                                                                           1]))) + 1), sep = " ")) + ggplot2::geom_smooth(method = "loess",
                                                                                                                                          se = F, color = "black", fullrange = F)
  plot
}








plot_coverage=function(biomart_file,fastaFile,chr_name,heteroplasmy_matrix){
  base_name=rep(0,length(colnames(heteroplasmy_matrix)))
  for ( i in 1:length(colnames(heteroplasmy_matrix_ci))){
    base_name[i]=strsplit(colnames(heteroplasmy_matrix_ci),"_")[[i]][1]
  }
  colnames(heteroplasmy_matrix)=base_name
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df_bg <- data.frame(seq_name, sequence)
  name_sequences=rep(0,length(df_bg$seq_name))
  for (i in 1:length(df_bg$seq_name)){
    name_sequences[i]=strsplit(df_bg$seq_name[i],split=" ")[[1]][1]
  }

  df_bg$seq_name=name_sequences

  length_sequence=rep(0,length(df_bg$seq_name))
  for (q in 1:length(df_bg$seq_name)){
    length_sequence[q]=length(unlist(strsplit(df_bg$sequence[q],split=NULL)))

  }

  length_chr=length_sequence[df_bg$seq_name=="MT"]

  mt_name=biomart_file
  mt_name=mt_name[mt_name$Chromosome.scaffold.name==chr_name,]
  nomi=mt_name$Gene.name
  p=rep(0,length(nomi))
  for (j in 1:length(nomi)){
    i=which(mt_name$Gene.name==nomi[j])
    bo=seq(mt_name$Gene.start..bp.[i],mt_name$Gene.end..bp.[i])
    bo=as.character(bo)

    if (sum(colnames(heteroplasmy_matrix)%in%bo)>0){
      p[i]=nomi[i]
    }

  }
  mt_name_small=mt_name[mt_name$Gene.name%in%p,]

  start_plot=c(mt_name_small$Gene.start..bp.,mt_name_small$Gene.end..bp.)
  nomi_plot=paste0(mt_name_small$Gene.name,"-start")
  nomi_plot_2=paste0(mt_name_small$Gene.name,"-end")
  names(start_plot)=c(nomi_plot,nomi_plot_2)
  start_plot=sort(start_plot)


  mt_name_new=mt_name[!(mt_name$Gene.stable.ID%in%mt_name_small$Gene.stable.ID),]

  start_plot_new=c(mt_name_new$Gene.start..bp.,mt_name_new$Gene.end..bp.)
  nomi_plot=paste0(mt_name_new$Gene.name,"-start")
  nomi_plot_2=paste0(mt_name_new$Gene.name,"-end")
  names(start_plot_new)=c(nomi_plot,nomi_plot_2)
  start_plot_new=sort(start_plot_new)


  colnames(heteroplasmy_matrix)=base_name
  all_base=seq(1,length_chr)
  all_base=as.character(all_base)

  name_all_base=rep(0,length_chr)
  name_all_base[all_base%in%colnames(heteroplasmy_matrix)]="Covered"
  name_all_base[!all_base%in%colnames(heteroplasmy_matrix)]="Not_Covered"
  name_all_base=paste0(name_all_base,all_base)
  colore=rep(0,length_chr)
  colore[all_base%in%colnames(heteroplasmy_matrix)]="red"
  colore[!all_base%in%colnames(heteroplasmy_matrix)]="black"


  custom.genome <- regioneR::toGRanges(data.frame(chr=c(chr_name), start=c(1), end=c(length_chr)))

  kp <- karyoploteR::plotKaryotype(genome = custom.genome,plot.type = 2)

  df <- data.frame(chr=c(rep(chr_name,length(all_base))), start=all_base, end=all_base)
  karyoploteR::kpPlotRegions(kp, data=df, col=colore, border=colore, r0=0, r1=0.10)



  markers <- data.frame(chr=rep(chr_name, length(start_plot)), pos=start_plot, labels=names(start_plot))
  markers_new <- data.frame(chr=rep(chr_name, length(start_plot_new)), pos=start_plot_new, labels=names(start_plot_new))
  karyoploteR::kpAddBaseNumbers(kp,tick.dist=1000,minor.tick.dist=1000,add.units=T,tick.len = 30)

  karyoploteR::kpPlotMarkers(kp, chr=markers$chr, x=markers$pos, labels=markers$labels,r0=0.10,r1=0.55,cex=0.3,label.color = rep("red",length(markers_new)),label.margin = 1,line.color = rep("red",length(markers_new)))
  karyoploteR::kpPlotMarkers(kp, chr=markers_new$chr, x=markers_new$pos, labels=markers_new$labels,r0=0.10,r1=0.55,cex=0.3,label.dist = 0.000001,data.panel=2,marker.parts = c(0, 0.9, 0.1),label.margin = 0.0001)

  legend(x = "topleft", fill = c("red", "black"), legend = c("Covered", "Not covered"),cex=0.5)

}





