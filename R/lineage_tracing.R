

clustering_dist_ang=function (heteroplasmy_matrix, allele_matrix, cluster, top_pos, deepSplit_param,minClusterSize_param,threshold=0.20,min_value,index=NULL,relevant_bases= NULL,max_frac=0.70)
{




  if (is.null(index)){


    entropia_matrix=heteroplasmy_matrix
    res_ang = rep(0, length(colnames(entropia_matrix)))
    res_ang = as.list(res_ang)
    names(res_ang) = colnames(entropia_matrix)
    j = 1
    for (i in 1:length(colnames(entropia_matrix))) {
      allele_matrix_2_1 = allele_matrix[, (j:(j + 3))]
      dist_ang_pos = rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] = dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])]=0
      j = j + 3 + 1
    }
    res_ang_square = lapply(res_ang, FUN = function(x) {
      return(x^2)
    })
    if (!is.null(relevant_bases)) {
      res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                           relevant_bases)]
    }

    if (is.null(relevant_bases)) {
      var_dist = lapply(res_ang_square, FUN = function(x) {
        return(var(as.vector(as.dist(x))))
      })

      mean_dist = lapply(res_ang_square, FUN = function(x) {
        return(mean(as.vector(as.dist(x))))
      })


      var_dist = as.numeric(var_dist)



      names(var_dist) = names(res_ang)

      var_dist = var_dist[order(var_dist, decreasing = T)]


      mean_max=apply(entropia_matrix[,names(var_dist)],2,function(x){
        #x=mean(x)
        #x=x[x>=threshold]

        if (sum(x>=threshold)>max_frac*length(x)){return(FALSE)
        }


        else(return(TRUE))
      })

      var_dist=var_dist[as.vector(mean_max)]




      if (length(var_dist)>=top_pos){
        var_dist_top_all = var_dist[1:top_pos]}
      if (length(var_dist)<top_pos){
        var_dist_top_all = var_dist
      }
      var_dist_sum=rep(0,length(var_dist_top_all))
      for ( i in 1:length(var_dist_top_all)){
        var_dist_sum[i]=var_dist_top_all[i]/sum(var_dist_top_all)
      }

      if (var_dist_sum[1]>min_value){
        var_dist_top=var_dist_top_all[var_dist_sum>min_value]
      }

      if (var_dist_sum[1]<=min_value){print(paste0("No variance above ",min_value,".","Using all features."))
        var_dist_top=var_dist_top_all
      }



      top_dist = names(var_dist_top)
      res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                           top_dist)]
    }



    res_ang_sum = Reduce("+", res_ang_sel)
    res_ang_sqrt = sqrt(res_ang_sum)
    dist_ang_sqrt = as.dist(res_ang_sqrt)
    my.tree <- hclust(dist_ang_sqrt)
    my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label
    length(unique(my.clusters))
    classe <- data.frame(old_classification = cluster, new_classification = my.clusters)
    if (!is.null(relevant_bases)){
      top_dist=NULL}
    common_idx=NULL
    return(list(classe, res_ang_sqrt,top_dist,common_idx))
  }


  if (!is.null(index)){

    entropia_matrix=heteroplasmy_matrix
    res_ang = rep(0, length(colnames(entropia_matrix)))
    res_ang = as.list(res_ang)
    names(res_ang) = colnames(entropia_matrix)
    j = 1
    for (i in 1:length(colnames(entropia_matrix))) {
      allele_matrix_2_1 = allele_matrix[index[[colnames(entropia_matrix)[i]]], (j:(j + 3))]
      dist_ang_pos = rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] = dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])]=0
      j = j + 3 + 1
    }
    res_ang_square = lapply(res_ang, FUN = function(x) {
      return(x^2)
    })
    if (!is.null(relevant_bases)) {
      res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                           relevant_bases)]
    }

    if (is.null(relevant_bases)) {
      var_dist = lapply(res_ang_square, FUN = function(x) {
        return(var(as.vector(as.dist(x))))
      })

      mean_dist = lapply(res_ang_square, FUN = function(x) {
        return(mean(as.vector(as.dist(x))))
      })


      var_dist = as.numeric(var_dist)



      names(var_dist) = names(res_ang)

      var_dist = var_dist[order(var_dist, decreasing = T)]



      mean_max=rep(0,length(names(var_dist)))
      for (i in 1:length(names(var_dist))){
        x=entropia_matrix[index[[names(var_dist)[i]]],names(var_dist)[i]]
        if (sum(x>=threshold)>max_frac*length(x)){mean_max[i]="FALSE"
        }


        else(mean_max[i]="TRUE")
      }




      var_dist=var_dist[as.logical(mean_max)]




      if (length(var_dist)>=top_pos){
        var_dist_top_all = var_dist[1:top_pos]}

      if (length(var_dist)<top_pos){
        var_dist_top_all = var_dist}

      var_dist_sum=rep(0,length(var_dist_top_all))
      for ( i in 1:length(var_dist_top_all)){
        var_dist_sum[i]=var_dist_top_all[i]/sum(var_dist_top_all)
      }


      if (var_dist_sum[1]>min_value){
        var_dist_top=var_dist_top_all[var_dist_sum>min_value]
      }

      if (var_dist_sum[1]<=min_value){print(paste0("No variance above ",min_value,".","Using all features."))
        var_dist_top=var_dist_top_all
      }






      top_dist = names(var_dist_top)
      res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                           top_dist)]
    }

    bases_example=names(res_ang_sel)
    common_idx=rep(list(0),length(bases_example))


    for (i in 1:length(bases_example)){

      common_idx[[i]]=as.numeric(index[[bases_example[i]]])

    }
    common_idx=Reduce(intersect, common_idx)
    common_cells=row.names(entropia_matrix)[common_idx]
    res_ang_sel_com=res_ang_sel
    for (i in 1:length(res_ang_sel)){
      res_ang_sel_com[[i]]=res_ang_sel[[i]][common_cells,common_cells]
    }
    res_ang_sel=res_ang_sel_com
    res_ang_sum = Reduce("+", res_ang_sel)
    res_ang_sqrt = sqrt(res_ang_sum)
    dist_ang_sqrt = as.dist(res_ang_sqrt)
    my.tree <- hclust(dist_ang_sqrt)
    my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label
    length(unique(my.clusters))
    classe <- data.frame(old_classification = cluster[common_idx], new_classification = my.clusters)
    if (!is.null(relevant_bases)){
      top_dist=NULL}
    return(list(classe, res_ang_sqrt,top_dist,common_idx))
  }

}



choose_features_clustering=function (heteroplasmy_matrix, allele_matrix, cluster, top_pos, deepSplit_param,
                                     minClusterSize_param,min_value,threshold,index=NULL,max_frac=0.70)
{




entropia_matrix=heteroplasmy_matrix
  if (is.null(index)){
    res_ang = rep(0, length(colnames(entropia_matrix)))
    res_ang = as.list(res_ang)
    names(res_ang) = colnames(entropia_matrix)
    j = 1
    for (i in 1:length(colnames(entropia_matrix))) {
      allele_matrix_2_1 = allele_matrix[, (j:(j + 3))]
      dist_ang_pos = rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] = dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])]=0
      j = j + 3 + 1
    }
    res_ang_square = lapply(res_ang, FUN = function(x) {
      return(x^2)
    })



    var_dist = lapply(res_ang_square, FUN = function(x) {
      return(var(as.vector(as.dist(x))))
    })

    mean_dist = lapply(res_ang_square, FUN = function(x) {
      return(mean(as.vector(as.dist(x))))
    })

    var_dist = as.numeric(var_dist)

    names(var_dist) = names(res_ang)

    var_dist = var_dist[order(var_dist, decreasing = T)]


    mean_max=apply(entropia_matrix[,names(var_dist)],2,function(x){


      if (sum(x>=threshold)>max_frac*length(x)){return(FALSE)
      }


      else(return(TRUE))
    })

    var_dist=var_dist[as.vector(mean_max)]


    number_pos=rep(list(0),length(min_value))
    for (i in 1:length(min_value)) {

      if(length(var_dist)>=top_pos){

        var_dist_top_all = var_dist[1:top_pos]
        var_dist_sum=rep(0,length(var_dist_top_all))
        for ( j in 1:length(var_dist_top_all)){
          var_dist_sum[j]=var_dist_top_all[j]/sum(var_dist_top_all)
        }
        if(min_value[i]>=var_dist_sum[1]){print(paste0("No variance above ",min_value[i],".","Using all features"))
          var_dist_top=var_dist_top_all }

        if(min_value[i]<var_dist_sum[1]){
          var_dist_top=var_dist_top_all[var_dist_sum>min_value[i]]


        }}
      if(length(var_dist)<top_pos)
      { var_dist_top_all = var_dist
      var_dist_sum=rep(0,length(var_dist_top_all))
      for ( j in 1:length(var_dist_top_all)){
        var_dist_sum[j]=var_dist_top_all[j]/sum(var_dist_top_all)
      }
      if(min_value[i]>=var_dist_sum[1]){print(paste0("No variance above ",min_value[i],".","Using all features"))
        var_dist_top=var_dist_top_all  }

      if(min_value[i]<var_dist_sum[1]){
        var_dist_top=var_dist_top_all[var_dist_sum>min_value[i]]


      }
      }

      top_dist = names(var_dist_top)
      res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                           top_dist)]
      res_ang_sum = Reduce("+", res_ang_sel)
      res_ang_sqrt = sqrt(res_ang_sum)
      dist_ang_sqrt = as.dist(res_ang_sqrt)
      my.tree <- hclust(dist_ang_sqrt)
      my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                  deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label






      number_pos[[i]]=my.clusters

      names(number_pos[[i]])=rep(paste0("N",i),length(number_pos[[i]]))
    }
    df <- data.frame(matrix(unlist(number_pos), ncol=length(number_pos), byrow=FALSE),stringsAsFactors=FALSE)
    colnames(df)=unique(names(unlist(number_pos)))
    clustree::clustree(df, prefix = "N",node_size = 10,layout = "sugiyama") +
      guides(edge_colour = FALSE, edge_alpha = FALSE) +
      theme(legend.position = "bottom")



  }


  if(!is.null(index)){

    res_ang = rep(0, length(colnames(entropia_matrix)))
    res_ang = as.list(res_ang)
    names(res_ang) = colnames(entropia_matrix)
    j = 1
    for (i in 1:length(colnames(entropia_matrix))) {
      allele_matrix_2_1 = allele_matrix[index[[colnames(entropia_matrix)[i]]], (j:(j + 3))]
      dist_ang_pos = rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] = dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])]=0
      j = j + 3 + 1
    }
    res_ang_square = lapply(res_ang, FUN = function(x) {
      return(x^2)
    })



    var_dist = lapply(res_ang_square, FUN = function(x) {
      return(var(as.vector(as.dist(x))))
    })

    mean_dist = lapply(res_ang_square, FUN = function(x) {
      return(mean(as.vector(as.dist(x))))
    })

    var_dist = as.numeric(var_dist)

    names(var_dist) = names(res_ang)

    var_dist = var_dist[order(var_dist, decreasing = T)]


    mean_max=rep(0,length(names(var_dist)))
    for (i in 1:length(names(var_dist))){
      x=entropia_matrix[index[[names(var_dist)[i]]],names(var_dist)[i]]
      if (sum(x>=threshold)>max_frac*length(x)){mean_max[i]="FALSE"
      }


      else(mean_max[i]="TRUE")
    }


    var_dist=var_dist[as.logical(mean_max)]


    number_pos=rep(list(0),length(min_value))
    number_idx=rep(list(0),length(min_value))
    for (i in 1:length(min_value)) {

      if(length(var_dist)>=top_pos){

        var_dist_top_all = var_dist[1:top_pos]
        var_dist_sum=rep(0,length(var_dist_top_all))
        for ( j in 1:length(var_dist_top_all)){
          var_dist_sum[j]=var_dist_top_all[j]/sum(var_dist_top_all)
        }
        if(min_value[i]>=var_dist_sum[1]){print(paste0("No variance above ",min_value[i],".","Using all features"))
          var_dist_top=var_dist_top_all }

        if(min_value[i]<var_dist_sum[1]){
          var_dist_top=var_dist_top_all[var_dist_sum>min_value[i]]


        }}
      if(length(var_dist)<top_pos)
      { var_dist_top_all = var_dist
      var_dist_sum=rep(0,length(var_dist_top_all))
      for ( j in 1:length(var_dist_top_all)){
        var_dist_sum[j]=var_dist_top_all[j]/sum(var_dist_top_all)
      }
      if(min_value[i]>=var_dist_sum[1]){print(paste0("No variance above ",min_value[i],".","Using all features"))
        var_dist_top=var_dist_top_all  }

      if(min_value[i]<var_dist_sum[1]){
        var_dist_top=var_dist_top_all[var_dist_sum>min_value[i]]


      }
      }

      top_dist = names(var_dist_top)
      res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                           top_dist)]



      bases_example=names(res_ang_sel)
      common_idx=rep(list(0),length(bases_example))


      for (l in 1:length(bases_example)){

        common_idx[[l]]=as.numeric(index[[bases_example[l]]])

      }
      common_idx=Reduce(intersect, common_idx)
      common_cells=row.names(entropia_matrix)[common_idx]
      res_ang_sel_com=res_ang_sel
      for (k in 1:length(res_ang_sel)){
        res_ang_sel_com[[k]]=res_ang_sel[[k]][common_cells,common_cells]
      }
      res_ang_sel=res_ang_sel_com
      res_ang_sum = Reduce("+", res_ang_sel)
      res_ang_sqrt = sqrt(res_ang_sum)
      dist_ang_sqrt = as.dist(res_ang_sqrt)
      my.tree <- hclust(dist_ang_sqrt)
      my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                  deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label


      number_pos[[i]]=my.clusters

      names(number_pos[[i]])=rep(paste0("N",i),length(number_pos[[i]]))

      number_idx[[i]]=common_idx

      names(number_idx[[i]])=common_cells
    }








    common_idx=Reduce(intersect, number_idx)
    for (i in 1:length(number_pos)){

      number_pos[[i]]=number_pos[[i]][number_idx[[i]]%in%common_idx]
    }

    df <- data.frame(matrix(unlist(number_pos), ncol=length(number_pos), byrow=FALSE),stringsAsFactors=FALSE)
    colnames(df)=unique(names(unlist(number_pos)))
    clustree::clustree(df, prefix = "N",node_size = 10,layout = "sugiyama") +
      guides(edge_colour = FALSE, edge_alpha = FALSE) +
      theme(legend.position = "bottom")



  }





}


heatmap_plot=function (marker_plot, marker_plot_plot, new_classification,
                       old_classification, res_ang_sqrt, cluster_columns = F, cluster_rows = T,
                       name_legend)
{
  norm_es=res_ang_sqrt
  cluster_unique = unique(new_classification)
  cluster_unique = sort(cluster_unique, decreasing = F)
  index = as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    index[[i]] = which(new_classification == cluster_unique[i])
  }
  condition_unique = as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    condition_unique[[i]] = old_classification[index[[i]]]
  }
  cell_all_day = colnames(norm_es)
  cell_unique = as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    cell_unique[[i]] = cell_all_day[index[[i]]]
  }
  cluster_unique_2 = as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    cluster_unique_2[[i]] = new_classification[index[[i]]]
  }
  cluster_finale = unlist(cluster_unique_2)
  good_finale = unlist(cell_unique)
  medium_finale = unlist(condition_unique)



  norm_es_plot = norm_es[good_finale, good_finale]

  color_cluster = rep(0, length(unique(new_classification)))
  for (i in 1:length(color_cluster)) {
    color_cluster[i] = gg_color_hue(length(color_cluster))[i]
  }
  names(color_cluster) = as.character(cluster_unique)
  color_condition = rep(0, length(unique(old_classification)))
  for (i in 1:length(color_condition)) {
    color_condition[i] = gg_color_hue(length(color_condition))[i]
  }
  names(color_condition) = as.character(unique(old_classification))
  Condition_factor = factor(old_classification)
  names(color_condition) = as.character(levels(Condition_factor))
  data_heatmap = data.frame(New_classification = cluster_finale,
                            Old_classification = medium_finale)
  haCol1 = ComplexHeatmap::HeatmapAnnotation(df = data_heatmap, col = list(New_classification = color_cluster,
                                                                           Old_classification = color_condition), show_legend = T)
  ht21 = ComplexHeatmap::Heatmap(as.matrix(norm_es_plot), cluster_rows =cluster_rows,
                                 col = circlize::colorRamp2(c(0, (max(norm_es_plot))), c("white",
                                                                                         "red")), name = name_legend, cluster_columns = cluster_columns,
                                 top_annotation = haCol1, row_names_gp = grid::gpar(fontsize = 6),
                                 show_column_names = F, show_row_names = T,clustering_method_rows = "complete",clustering_method_columns="complete")
  ComplexHeatmap::draw(ht21)
}



Plot_distance_matrix=function(dist_matrix_sc,cluster){

  color_cluster = rep(0, length(unique(cluster)))
  for (i in 1:length(color_cluster)) {
    color_cluster[i] = gg_color_hue(length(color_cluster))[i]
  }
  names(color_cluster) = as.character(unique(cluster))



  data_heatmap=data.frame(cluster=cluster)

  haCol1 = ComplexHeatmap::HeatmapAnnotation(df = data_heatmap, col = list(cluster = color_cluster), show_legend = T,show_annotation_name=F)
  ht21 = ComplexHeatmap::Heatmap(dist_matrix_sc, cluster_rows =T,
                                 col = circlize::colorRamp2(c(0, (max(dist_matrix_sc))), c("white",
                                                                                           "red")), name = "Euclidean Distance", cluster_columns = T,
                                 top_annotation = haCol1, row_names_gp = grid::gpar(fontsize = 6),
                                 show_column_names = F, show_row_names = F)
  ComplexHeatmap::draw(ht21)
}

vi_comparison=function(old_classification,new_classification,number_iter){
  random_vi=rep(list(0),number_iter)
  for (i in 1:number_iter){
    set.seed(i)

    random_clas=sample(old_classification,length(old_classification))
    random_vi[[i]]=mcclust::vi.dist(old_classification,random_clas)

  }
  final_random_vi=unlist(random_vi)
  unsupervised_vi=mcclust::vi.dist(old_classification,new_classification)
  empirical_pvalue=sum(final_random_vi<=unsupervised_vi)/length(final_random_vi)
  print(paste0("Empirical p value:",empirical_pvalue))
  return(empirical_pvalue)
}



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}






