

clustering_dist_ang=function (entropia_matrix, allele_matrix, cluster, top_pos, deepSplit_param,
                              minClusterSize_param, relevant_bases = NULL)
{
  res_ang = rep(0, length(colnames(entropia_matrix)))
  res_ang = as.list(res_ang)
  names(res_ang) = colnames(entropia_matrix)
  j = 1
  for (i in 1:length(colnames(entropia_matrix))) {
    allele_matrix_2_1 = allele_matrix[, (j:(j + 3))]
    dist_ang_pos = rdist::pdist(allele_matrix_2_1, metric = "angular")
    res_ang[[i]] = dist_ang_pos
    j = j + 3 + 1
  }
  res_ang_square = lapply(res_ang, FUN = function(x) {
    return(x^2)
  })
  if (!is.null(relevant_bases)) {
    res_ang_sel = res_ang_square[which(names(res_ang) %in%
                                         relevant_bases)]
  }
  else {
    var_dist = lapply(res_ang_square, FUN = function(x) {
      return(var(as.vector(as.dist(x))))
    })
    var_dist = as.numeric(var_dist)
    names(var_dist) = names(res_ang)
    var_dist = var_dist[order(var_dist, decreasing = T)]
    var_dist_top = var_dist[1:top_pos]
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
  return(list(classe, res_ang_sqrt))
}



heatmap_plot=function (marker_plot, marker_plot_plot, new_classification,
                       old_classification, norm_es, cluster_columns = F, cluster_rows = T,
                       name_legend)
{
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
  norm_es_plot = norm_es[marker_plot, good_finale]
  row.names(norm_es_plot) <- marker_plot_plot
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
  ht21 = ComplexHeatmap::Heatmap(as.matrix(norm_es_plot), cluster_rows = cluster_rows,
                                 col = circlize::colorRamp2(c(0, (max(norm_es_plot))), c("white",
                                                                                         "red")), name = name_legend, cluster_columns = cluster_columns,
                                 top_annotation = haCol1, row_names_gp = grid::gpar(fontsize = 6),
                                 show_column_names = F, show_row_names = T)
  ComplexHeatmap::draw(ht21)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}






