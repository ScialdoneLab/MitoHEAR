dpt_test=function (heteroplasmy_matrix, time, index = NULL, method = "GAM")
{
  if (is.null(index)) {
    position = colnames(heteroplasmy_matrix)
    p_value = rep(0, length(position))
    for (i in 1:length(position)) {
      Y = data.frame(t(heteroplasmy_matrix[, position[i]]))
      t = time
      if (method != "GAM") {
        p_value[i] = cor.test(as.numeric(Y), t, method = method)$p.value
      }
      if (method == "GAM") {
        gam.res <- apply(Y, 1, function(z) {
          z = as.numeric(as.vector(z))
          d <- data.frame(z = z, t = t)
          tmp <- gam::gam(z ~ gam::lo(t), data = d)
          p <- summary(tmp)[4][[1]][1, 5]
          f <- fitted(tmp)
          c(p, f)
        })
        genes.table <- data.frame(genes.names = rownames(Y))
        genes.table$pvals <- gam.res[1, ]
        genes.table$FDR <- p.adjust(genes.table$pvals,
                                    method = "fdr")
        genes.table <- genes.table[order(genes.table$FDR),
                                   ]
        genes.table$genes.names <- as.character(genes.table$genes.names)
        row.names(genes.table) = genes.table$genes.names
        results.gam.tot = genes.table
        gam.fitted <- gam.res[-1, ]
        results_gam_nomi = results.gam.tot$genes.names
        gam_fitted = gam.fitted
        p_value[i] = results.gam.tot$pvals
      }
    }
  }
  else {
    position = colnames(heteroplasmy_matrix)
    p_value = rep(0, length(position))
    for (i in 1:length(position)) {
      Y = data.frame(t(heteroplasmy_matrix[as.numeric(index[[position[i]]]),
                                           position[i]]))
      t = time[as.numeric(index[[position[i]]])]
      if (method != "GAM") {
        p_value[i] = cor.test(as.numeric(Y), t, method = method)$p.value
      }
      if (method == "GAM") {
        gam.res <- apply(Y, 1, function(z) {
          z = as.numeric(as.vector(z))
          d <- data.frame(z = z, t = t)
          tmp <- gam::gam(z ~ gam::lo(t), data = d)
          p <- summary(tmp)[4][[1]][1, 5]
          f <- fitted(tmp)
          c(p, f)
        })
        genes.table <- data.frame(genes.names = rownames(Y))
        genes.table$pvals <- gam.res[1, ]
        genes.table$FDR <- p.adjust(genes.table$pvals,
                                    method = "fdr")
        genes.table <- genes.table[order(genes.table$FDR),
                                   ]
        genes.table$genes.names <- as.character(genes.table$genes.names)
        row.names(genes.table) = genes.table$genes.names
        results.gam.tot = genes.table
        gam.fitted <- gam.res[-1, ]
        results_gam_nomi = results.gam.tot$genes.names
        gam_fitted = gam.fitted
        p_value[i] = results.gam.tot$pvals
      }
    }
  }
  fdr_update = p.adjust(p_value, method = "fdr")
  position_fdr = data.frame(fdr_update, colnames(heteroplasmy_matrix))
  colnames(position_fdr) = c("FDR_value", "Position")
  row.names(position_fdr) = colnames(heteroplasmy_matrix)
  position_fdr = position_fdr[order(position_fdr$FDR_value),
                              ]
  return(position_fdr)
}

get_wilcox_test=function(matrix,cluster,label_1,label_2,index=NULL){
  base=colnames(matrix)
  distribution=rep(0,length(base))
  if(is.null(index)){
    for (i in 1:length(base)){
      Y=as.vector(matrix[,base[i]])
      Y_label_1=Y[cluster==label_1]
      Y_label_2=Y[cluster==label_2]
      Y_label_all=c(Y_label_1,Y_label_2)
      name_label_1=rep(label_1,length(Y_label_1))
      name_label_2=rep(label_2,length(Y_label_2))
      name_label_all=c(name_label_1,name_label_2)
      data_test=data.frame(Y_label_all,name_label_all)
      colnames(data_test)=c("Heteroplasmy","Cluster")

      res <- wilcox.test( Heteroplasmy~ Cluster, data = data_test)
      distribution[i]=res$p.value
      if(is.nan(res$p.value)&all(Y_label_all==0)){
        distribution[i]=1
      }
      names(distribution)[i]=base[i]




    }}



  else{

    for (i in 1:length(base)){
      index_cell=index[[which(names(index)==base[i])]]

      cluster_index=cluster[index_cell]
      Y=as.vector(matrix[index_cell,base[i]])
      Y_label_1=Y[cluster_index==label_1]
      Y_label_2=Y[cluster_index==label_2]
      Y_label_all=c(Y_label_1,Y_label_2)
      name_label_1=rep(label_1,length(Y_label_1))
      name_label_2=rep(label_2,length(Y_label_2))
      name_label_all=c(name_label_1,name_label_2)
      data_test=data.frame(Y_label_all,name_label_all)
      colnames(data_test)=c("Heteroplasmy","Cluster")

      res <- wilcox.test( Heteroplasmy~ Cluster, data = data_test)
      distribution[i]=res$p.value
      if(is.nan(res$p.value)&all(Y_label_all==0)){
        distribution[i]=1
      }
      names(distribution)[i]=base[i]

    }}
  distribution_adjusted=p.adjust(distribution,method="fdr")
  return(distribution_adjusted)
}









