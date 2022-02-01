#' dpt_test
#' @param time Vector of diffusion pseudo time.
#' @param index index returned by \emph{get_heteroplasmy}.
#' @param method Character name denoting the method to choose for assigning an
#' adjusted p value to each of the bases. Can be one of GAM,pearson and
#' spearman. GAM: For each base, a GAM fit with formula z ~ lo(t) is performed
#' between the heteroplasmy values (z) and the time (t). The p value from the
#' table "Anova for Parametric Effects" is then assigned to the base.
#' pearson,spearman:for each base, a pearson or spearman correlation test is
#' performed between the heteroplasmy values and the time . The p value
#' obtained from the test is then assigned to the base. In all the three
#' possible methods, all the p values are then corrected with the method FDR.
#' @inheritParams plot_heteroplasmy
#' @return A data frame with 2 columns and number of rows equal to n_col in
#' \emph{heteroplasmy_matrix}. In the first column there are the names of the bases
#' while in the second column there are the adjusted p value.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/gam/versions/1.20/topics/gam}
#' @export dpt_test
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
        if (! requireNamespace("gam", quietly = TRUE)) {
          stop("Package gam needed for method==gam. Please install it: install.packages('gam')")
        }
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
        if (! requireNamespace("gam", quietly = TRUE)) {
          stop("Package gam needed for method==gam. Please install it: install.packages('gam')")
        }
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





#' get_wilcox_test
#' @param label_1 Character name of a first label included in cluster. It
#' denotes the first group used for the Wilcoxon test
#' @param label_2 Character name of a second label included in cluster and
#' different from label_1. it denotes the second group used for the Wilcoxon
#' test.
#' @inheritParams plot_heteroplasmy
#' @return It returns a vector of length equal to n_row in matrix. Each element
#' stands for a base and it contains the adjusted p-value (FDR), obtained in
#' unpaired two-samples Wilcoxon test from the comparison of the heteroplasmy
#' between the label_1 and label_2 group.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test}
#' @export get_wilcox_test
get_wilcox_test=function(heteroplasmy_matrix,cluster,label_1,label_2,index=NULL){
  base=colnames(heteroplasmy_matrix)
  distribution=rep(0,length(base))
  if(is.null(index)){
    for (i in 1:length(base)){
      Y=as.vector(heteroplasmy_matrix[,base[i]])
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
      Y=as.vector(heteroplasmy_matrix[index_cell,base[i]])
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









