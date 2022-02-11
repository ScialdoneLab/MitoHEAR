#'clustering_angular_distance
#'
#' For each pair of samples and for each base, an angular distance matrix is
#' computed based on the four allele frequencies. Then only the angular
#' distances corresponding to the relevant_bases are kept. If relevant bases is
#' NULL, then only the angular distances corresponding to the bases with
#' relative distance variance among samples above \emph{min_value} are kept .
#' Finally the distance between each pair of samples is defined as the
#' euclidean distance of the angular distances corresponding to the bases that
#' pass the previous filtering step. On this final distance matrix, a
#' hierarchical clustering approach is performed using the function
#' \emph{cutreeHybrid} of the package \emph{dynamicTreeCut}.
#' @param top_pos Numeric value. Number of bases sorted with decreasing values
#' of distance variance (see section \emph{Details} below) among samples. If
#' \emph{relevant_bases}=NULL, then the bases for performing hierarchical clustering
#' are the ones whose relative variance (variance of the base divided sum of
#' variance among \emph{top_pos} bases) is above \emph{min_value}.
#' @param deepSplit_param Integer value between 0 and 4 for the \emph{deepSplit}
#' parameter of the function \emph{cutreeHybrid}. See section \emph{Details}
#' below.
#' @param minClusterSize_param Integer value specifyng the \emph{minClusterSize}
#' parameter of the function \emph{cutreeHybrid}. See section \emph{Details}
#' below.
#' @param threshold Numeric value. If a base has heteroplasmy greater or equal
#' to \emph{threshold} in more than \emph{max_frac} of cells, then the base is
#' not considered for down stream analysis.
#' @param min_value Numeric value. If \emph{relevant_bases}=NULL, then the bases for
#' performing hierarchical clustering are the ones whose relative variance
#' (variance of the base divided sum of variance among \emph{top_pos} bases) is
#' above \emph{min_value}.
#' @param relevant_bases Character vector of bases to consider as features for
#' performing hierarchical clustering on samples.Default=NULL.
#' @param max_frac Numeric value.If a base has heteroplasmy greater or equal to
#' \emph{threshold} in more than \emph{max_frac} of cells, then the base is not
#' considered for down stream analysis.
#' @inheritParams plot_heteroplasmy
#' @inheritParams plot_allele_frequency
#' @return It returns a list with 4 elements: \item{classification}{ Dataframe with
#' two columns and n_row equal to n_row in heteroplasmy_matrix. The first
#' column is the old cluster annotation provided by cluster. The second columns
#' is the new cluster annotation obtained with hierarichal clustering on
#' distance matrix based on heteroplasmy values.
#'
#' }
#'
#' \item{dist_ang_matrix}{ Distance matrix based on heteroplasmy values as defined
#' in the section \emph{Details}
#'
#' }
#'
#' \item{top_bases_dist}{ Vector of bases used for hierarchical clustering. If
#' \emph{relevant_bases} is not NULL, then \emph{top_bases_dist}=NULL
#'
#' }
#'
#' \item{common_idx}{ Vector of indeces of samples for which hierarchical
#' clustering is performed. If \emph{index} is NULL, then \emph{common_idx}=NULL
#'
#' }
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/dynamicTreeCut/versions/1.63-1/topics/cutreeHybrid}
#' @export clustering_angular_distance
clustering_angular_distance <- function (heteroplasmy_matrix, allele_matrix, cluster, top_pos, deepSplit_param, minClusterSize_param, threshold = 0.20, min_value, index, relevant_bases =  NULL, max_frac = 0.70) {
  if (is.null(index)) {



    res_ang <- rep(0, length(colnames(heteroplasmy_matrix)))
    res_ang <- as.list(res_ang)
    names(res_ang) <- colnames(heteroplasmy_matrix)
    j <- 1
    for (i in 1:length(colnames(heteroplasmy_matrix))) {
      allele_matrix_2_1 <- allele_matrix[, (j:(j + 3))]
      dist_ang_pos <- rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] <- dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])]<-0
      j <- j + 4
    }
    res_ang_square <- lapply(res_ang, FUN = function(x) {
      return(x^2)
    })
    if (!is.null(relevant_bases)) {
      res_ang_sel <- res_ang_square[which(names(res_ang) %in% relevant_bases)]
    }

    if (is.null(relevant_bases)) {
      var_dist <- lapply(res_ang_square, FUN = function(x) {
        distance_vector <- as.vector(as.dist(x))
        if (length(distance_vector) == 1){
          return(0)
          }
        else{
          return(var(distance_vector))
        }
      })

      mean_dist <- lapply(res_ang_square, FUN = function(x) {
        return(mean(as.vector(as.dist(x))))
      })


      var_dist <- as.numeric(var_dist)



      names(var_dist) <- names(res_ang)

      var_dist <- var_dist[order(var_dist, decreasing = T)]


      mean_max <- apply(heteroplasmy_matrix[, names(var_dist)], 2, function(x){


        if (sum(x >= threshold) > max_frac*length(x)) {
          return(FALSE)
        }


        else(return(TRUE))
      })

      if (all(!as.vector(mean_max))) {
        percentage_frac <- max_frac*100
        stop(paste0("All the bases have heteroplasmy greater than ", threshold, " in more than", percentage_frac, "% of cells"))
      }

      var_dist <- var_dist[as.vector(mean_max)]




      if (length(var_dist) >= top_pos) {
        var_dist_top_all <- var_dist[1:top_pos]}
      if (length(var_dist) < top_pos) {
        var_dist_top_all <- var_dist
        warning(paste0("Less than ", top_pos, " bases present. All the bases will be consider for downstream analysis."))
      }
      var_dist_sum <- rep(0, length(var_dist_top_all))

      if (sum(var_dist_top_all) == 0) {
        stop("All the bases has variance 0 across samples")
      }
      for ( i in 1:length(var_dist_top_all)) {
        var_dist_sum[i] <- var_dist_top_all[i]/sum(var_dist_top_all)
      }

      if (var_dist_sum[2] > min_value) {
        var_dist_top <- var_dist_top_all[var_dist_sum>min_value]
      }

      if (var_dist_sum[2] <= min_value){
        warning(paste0("There are not at least two bases with relative variance above  ", min_value, ".", "All the bases will be used."))
        var_dist_top <- var_dist_top_all
      }



      top_bases_dist <- names(var_dist_top)
      res_ang_sel <- res_ang_square[which(names(res_ang) %in% top_bases_dist)]
    }



    res_ang_sum <- Reduce("+", res_ang_sel)
    dist_ang_matrix <- sqrt(res_ang_sum)
    dist_ang_sqrt <- as.dist(dist_ang_matrix)
    my.tree <- hclust(dist_ang_sqrt)
    my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt), deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label
    classification <- data.frame(old_classification = cluster, new_classification = my.clusters)
    if (!is.null(relevant_bases)) {
      top_bases_dist <- NULL}
    common_idx <- NULL
    return(list(classification, dist_ang_matrix, top_bases_dist, common_idx))
  }


  if (!is.null(index)) {


    res_ang <- rep(0, length(colnames(heteroplasmy_matrix)))
    res_ang <- as.list(res_ang)
    names(res_ang) <- colnames(heteroplasmy_matrix)
    j <- 1
    for (i in 1:length(colnames(heteroplasmy_matrix))) {
      allele_matrix_2_1 <- allele_matrix[index[[colnames(heteroplasmy_matrix)[i]]], (j:(j + 3))]
      dist_ang_pos <- rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] <- dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])] <- 0
      j <- j + 4
    }
    res_ang_square <- lapply(res_ang, FUN = function(x) {
      return(x^2)
    })
    if (!is.null(relevant_bases)) {
      res_ang_sel <- res_ang_square[which(names(res_ang) %in% relevant_bases)]
    }

    if (is.null(relevant_bases)) {
      var_dist <- lapply(res_ang_square, FUN = function(x) {
        distance_vector <- as.vector(as.dist(x))
        if (length(distance_vector) == 1){
          return(0)
          }
        else{
          return(var(distance_vector))
        }

      })

      mean_dist <- lapply(res_ang_square, FUN = function(x) {
        return(mean(as.vector(as.dist(x))))
      })


      var_dist <- as.numeric(var_dist)



      names(var_dist) <- names(res_ang)

      var_dist <- var_dist[order(var_dist, decreasing = T)]



      mean_max <- rep(0, length(names(var_dist)))
      for (i in 1:length(names(var_dist))) {
        x <- heteroplasmy_matrix[index[[names(var_dist)[i]]], names(var_dist)[i]]
        if (sum(x >= threshold) > max_frac*length(x)){
          mean_max[i] <- "FALSE"
        }


        else(mean_max[i] <- "TRUE")
      }


      if (all(!as.logical(mean_max))){
        percentage_frac <- max_frac*100
        stop(paste0("All the bases have heteroplasmy greater than ", threshold, " in more than", percentage_frac, "% of cells"))
      }


      var_dist <- var_dist[as.logical(mean_max)]




      if (length(var_dist) >= top_pos){
        var_dist_top_all <- var_dist[1:top_pos]}

      if (length(var_dist) < top_pos){
        var_dist_top_all <- var_dist
        warning(paste0("Less than ", top_pos, " bases present. All the bases will be consider for downstream analysis."))
        }

      var_dist_sum <- rep(0, length(var_dist_top_all))

      if (sum(var_dist_top_all) == 0) {
        stop("All the bases has variance 0 across samples")
      }
      for ( i in 1:length(var_dist_top_all)) {
        var_dist_sum[i] <- var_dist_top_all[i]/sum(var_dist_top_all)
      }


      if (var_dist_sum[2] > min_value) {
        var_dist_top <- var_dist_top_all[var_dist_sum>min_value]
      }

      if (var_dist_sum[2] <= min_value) {
        var_dist_top <- var_dist_top_all
        warning(paste0("There are not at least two bases with relative variance above  ", min_value, ".", "All the bases will be used."))
      }






      top_bases_dist <- names(var_dist_top)
      res_ang_sel <- res_ang_square[which(names(res_ang) %in% top_bases_dist)]
    }

    bases_example <- names(res_ang_sel)
    common_idx <- rep(list(0), length(bases_example))


    for (i in 1:length(bases_example)) {

      common_idx[[i]] <- as.numeric(index[[bases_example[i]]])

    }
    common_idx <- Reduce(intersect, common_idx)
    if (length(common_idx) < 2) {
      stop("There are not at least 2 samples that cover all the filtered bases")
    }
    common_cells <- row.names(heteroplasmy_matrix)[common_idx]
    res_ang_sel_com <- res_ang_sel
    for (i in 1:length(res_ang_sel)){
      res_ang_sel_com[[i]] <- res_ang_sel[[i]][common_cells, common_cells]
    }
    res_ang_sel <- res_ang_sel_com
    res_ang_sum <- Reduce("+", res_ang_sel)
    dist_ang_matrix <- sqrt(res_ang_sum)
    dist_ang_sqrt <- as.dist(dist_ang_matrix)
    my.tree <- hclust(dist_ang_sqrt)
    my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label
    classification <- data.frame(old_classification = cluster[common_idx], new_classification = my.clusters)
    if (!is.null(relevant_bases)){
      top_bases_dist <- NULL
      }
    return(list(classification, dist_ang_matrix, top_bases_dist, common_idx))
  }

}







#' choose_features_clustering
#' @param min_value_vector Numeric vector. For each value in the vector, the function
#' \emph{clustering_angular_distance} is run with parameter \emph{min_value} equal to
#' one element of the vector \emph{min_value_vector}.
#' @inheritParams plot_heteroplasmy
#' @inheritParams plot_allele_frequency
#' @inheritParams clustering_angular_distance
#' @return Clustree plot.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://cran.r-project.org/package=clustree}
#' @export choose_features_clustering
choose_features_clustering <- function (heteroplasmy_matrix, allele_matrix, cluster, top_pos, deepSplit_param,
                                     minClusterSize_param, min_value_vector, threshold = 0.20, index, max_frac = 0.70)
{

  if (!(requireNamespace("clustree", quietly = TRUE))) {
    stop("Package clustree needed for this function to work. Please install it:install.packages('clustree')")
  }





  if (is.null(index)) {
    res_ang <- rep(0, length(colnames(heteroplasmy_matrix)))
    res_ang <- as.list(res_ang)
    names(res_ang) <- colnames(heteroplasmy_matrix)
    j <- 1
    for (i in 1:length(colnames(heteroplasmy_matrix))) {
      allele_matrix_2_1 <- allele_matrix[, (j:(j + 3))]
      dist_ang_pos <- rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] <- dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])] <- 0
      j <- j + 4
    }
    res_ang_square <- lapply(res_ang, FUN = function(x) {
      return(x^2)
    })



    var_dist <- lapply(res_ang_square, FUN = function(x) {
      distance_vector <- as.vector(as.dist(x))
      if (length(distance_vector) == 1){
        return(0)
        }
      else{
        return(var(distance_vector))
      }

    })

    mean_dist <- lapply(res_ang_square, FUN = function(x) {
      return(mean(as.vector(as.dist(x))))
    })

    var_dist <- as.numeric(var_dist)

    names(var_dist) <- names(res_ang)

    var_dist <- var_dist[order(var_dist, decreasing = T)]


    mean_max <- apply(heteroplasmy_matrix[, names(var_dist)], 2, function(x){


      if (sum(x >= threshold) > max_frac*length(x)) {
        return(FALSE)
      }


      else(return(TRUE))
    })

    if (all(!as.vector(mean_max))) {
      percentage_frac <- max_frac*100
      stop(paste0("All the bases have heteroplasmy greater than ", threshold, " in more than", percentage_frac, "% of cells"))
    }
    var_dist <- var_dist[as.vector(mean_max)]


    number_pos <- rep(list(0), length(min_value_vector))
    for (i in 1:length(min_value_vector)) {

      if(length(var_dist) >= top_pos){

        var_dist_top_all <- var_dist[1:top_pos]
        var_dist_sum <- rep(0, length(var_dist_top_all))
        if (sum(var_dist_top_all) == 0){
          stop("All the bases has variance 0 across samples")
        }
        for ( j in 1:length(var_dist_top_all)){
          var_dist_sum[j] <- var_dist_top_all[j] / sum(var_dist_top_all)
        }
        if (min_value_vector[i] >= var_dist_sum[2]) {
          warning(paste0("There are not at least two bases with relative variance above  ", min_value_vector[i], ".", "All the bases will be used."))
          var_dist_top <- var_dist_top_all }

        if (min_value_vector[i] < var_dist_sum[2]) {
          var_dist_top <- var_dist_top_all[var_dist_sum>min_value_vector[i]]


        }}
      if (length(var_dist) < top_pos) {
        var_dist_top_all <- var_dist
        warning(paste0("Less than ", top_pos, " bases present. All the bases will be consider for downstream analysis."))
        var_dist_sum <- rep(0, length(var_dist_top_all))

      if (sum(var_dist_top_all) == 0){
        stop("All the bases has variance 0 across samples")
      }

      for ( j in 1:length(var_dist_top_all)) {
        var_dist_sum[j] <- var_dist_top_all[j] / sum(var_dist_top_all)
      }



      if (min_value_vector[i] >= var_dist_sum[2]) {
        warning(paste0("There are not at least two bases with relative variance above  ", min_value_vector[i], ".", "All the bases will be used."))
        var_dist_top <- var_dist_top_all
        }

      if (min_value_vector[i] < var_dist_sum[2]) {
        var_dist_top <- var_dist_top_all[var_dist_sum > min_value_vector[i]]


      }
      }

      top_bases_dist <- names(var_dist_top)
      res_ang_sel <- res_ang_square[which(names(res_ang) %in% top_bases_dist)]
      res_ang_sum <- Reduce("+", res_ang_sel)
      dist_ang_matrix <- sqrt(res_ang_sum)
      dist_ang_sqrt <- as.dist(dist_ang_matrix)
      my.tree <- hclust(dist_ang_sqrt)
      my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                  deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label






      number_pos[[i]] <- my.clusters

      names(number_pos[[i]]) <- rep(paste0("N", i), length(number_pos[[i]]))
    }
    df <- data.frame(matrix(unlist(number_pos), ncol = length(number_pos), byrow = FALSE), stringsAsFactors = FALSE)
    colnames(df) <- unique(names(unlist(number_pos)))
    clustree::clustree(df, prefix = "N", node_size = 10, layout = "sugiyama") +
      ggplot2::guides(edge_colour = FALSE, edge_alpha = FALSE) +
      ggplot2::theme(legend.position = "bottom")



  }


  if (!is.null(index)){

    res_ang <- rep(0, length(colnames(heteroplasmy_matrix)))
    res_ang <- as.list(res_ang)
    names(res_ang) <- colnames(heteroplasmy_matrix)
    j <- 1
    for (i in 1:length(colnames(heteroplasmy_matrix))) {
      allele_matrix_2_1 <- allele_matrix[index[[colnames(heteroplasmy_matrix)[i]]], (j:(j + 3))]
      dist_ang_pos <- rdist::pdist(allele_matrix_2_1, metric = "angular")
      res_ang[[i]] <- dist_ang_pos
      res_ang[[i]][is.na(res_ang[[i]])] <- 0
      j <- j + 4
    }
    res_ang_square <- lapply(res_ang, FUN = function(x) {
      return(x^2)
    })



    var_dist <- lapply(res_ang_square, FUN = function(x) {
      distance_vector <- as.vector(as.dist(x))
      if (length(distance_vector) == 1) {
        return(0)
        }
      else{
        return(var(distance_vector))
      }
    })

    mean_dist <- lapply(res_ang_square, FUN = function(x) {
      return(mean(as.vector(as.dist(x))))
    })

    var_dist <- as.numeric(var_dist)

    names(var_dist) <- names(res_ang)

    var_dist <- var_dist[order(var_dist, decreasing = T)]


    mean_max <- rep(0, length(names(var_dist)))
    for (i in 1:length(names(var_dist))){
      x <- heteroplasmy_matrix[index[[names(var_dist)[i]]],names(var_dist)[i]]
      if (sum(x >= threshold) > max_frac*length(x)) {
        mean_max[i] <- "FALSE"
      }


      else(mean_max[i] <- "TRUE")
    }

    if (all(!as.logical(mean_max))) {
      percentage_frac <- max_frac*100
      stop(paste0("All the bases have heteroplasmy greater than ", threshold, " in more than", percentage_frac, "% of cells"))
    }
    var_dist <- var_dist[as.logical(mean_max)]


    number_pos <- rep(list(0), length(min_value_vector))
    number_idx <- rep(list(0), length(min_value_vector))
    for (i in 1:length(min_value_vector)) {

      if (length(var_dist) >= top_pos) {

        var_dist_top_all <- var_dist[1:top_pos]
        var_dist_sum <- rep(0, length(var_dist_top_all))

        if (sum(var_dist_top_all) ==0 ) {
          stop("All the bases has variance 0 across samples")
        }

        for ( j in 1:length(var_dist_top_all)){
          var_dist_sum[j] <- var_dist_top_all[j] / sum(var_dist_top_all)
        }
        if (min_value_vector[i] >= var_dist_sum[2]) {
          warning(paste0("There are not at least two bases with relative variance above  ", min_value_vector[i], ".", "All the bases will be used."))
          var_dist_top <- var_dist_top_all }

        if (min_value_vector[i] < var_dist_sum[2]) {
          var_dist_top <- var_dist_top_all[var_dist_sum>min_value_vector[i]]


        }}
      if (length(var_dist)<top_pos) {
        var_dist_top_all <- var_dist
        warning(paste0("Less than ", top_pos, " bases present. All the bases will be consider for downstream analysis."))
        var_dist_sum <- rep(0, length(var_dist_top_all))
      if (sum(var_dist_top_all) == 0){
        stop("All the bases has variance 0 across samples")
      }

      for ( j in 1:length(var_dist_top_all)) {
        var_dist_sum[j] <- var_dist_top_all[j] / sum(var_dist_top_all)
      }
      if (min_value_vector[i] >= var_dist_sum[2]) {
        warning(paste0("There are not at least two bases with relative variance above ", min_value_vector[i], ".", "All the bases will be used."))
        var_dist_top <- var_dist_top_all  }

      if (min_value_vector[i] < var_dist_sum[2]) {
        var_dist_top <- var_dist_top_all[var_dist_sum > min_value_vector[i]]


      }
      }

      top_bases_dist <- names(var_dist_top)
      res_ang_sel <- res_ang_square[which(names(res_ang) %in% top_bases_dist)]



      bases_example <- names(res_ang_sel)
      common_idx <- rep(list(0), length(bases_example))


      for (l in 1:length(bases_example)) {

        common_idx[[l]] <- as.numeric(index[[bases_example[l]]])

      }
      common_idx <- Reduce(intersect, common_idx)
      if (length(common_idx) < 2){
      stop("There are not at least 2 samples that cover all the filtered bases")
        }
      common_cells <- row.names(heteroplasmy_matrix)[common_idx]
      res_ang_sel_com <- res_ang_sel
      for (k in 1:length(res_ang_sel)) {
        res_ang_sel_com[[k]] <- res_ang_sel[[k]][common_cells, common_cells]
      }
      res_ang_sel <- res_ang_sel_com
      res_ang_sum <- Reduce("+", res_ang_sel)
      dist_ang_matrix <- sqrt(res_ang_sum)
      dist_ang_sqrt <- as.dist(dist_ang_matrix)
      my.tree <- hclust(dist_ang_sqrt)
      my.clusters <- dynamicTreeCut::cutreeHybrid(my.tree, distM = as.matrix(dist_ang_sqrt),
                                                  deepSplit = deepSplit_param, minClusterSize = minClusterSize_param)$label


      number_pos[[i]] <- my.clusters

      names(number_pos[[i]]) <- rep(paste0("N", i), length(number_pos[[i]]))

      number_idx[[i]] <- common_idx

      names(number_idx[[i]]) <- common_cells
    }








    common_idx <- Reduce(intersect,number_idx)
    for (i in 1:length(number_pos)) {

      number_pos[[i]] <- number_pos[[i]][number_idx[[i]]%in%common_idx]
    }

    df <- data.frame(matrix(unlist(number_pos), ncol = length(number_pos), byrow = FALSE), stringsAsFactors = FALSE)
    colnames(df) <- unique(names(unlist(number_pos)))
    clustree::clustree(df, prefix = "N", node_size = 10, layout = "sugiyama") +
      ggplot2::guides(edge_colour = FALSE, edge_alpha = FALSE) +
      ggplot2::theme(legend.position = "bottom")



  }





}






#' plot_heatmap
#' @param new_classification Character vector.Second column of the dataframe
#' returned by function \emph{clustering_angular_distance} (first element of the
#' output).
#' @param old_classification Character vector. First column of the dataframe
#' returned by function \emph{clustering_angular_distance} (first element of the
#' output).
#' @param dist_ang_matrix Distance matrix obtained from \emph{clustering_angular_distance}
#' (second element of the output).
#' @param cluster_columns Logical. Parameter for cluster_columns argument of
#' the function \emph{Heatmap} in the package \emph{ComplexHeatmap}
#' @param cluster_rows Logical. Parameter for cluster_rows argument of the
#' function \emph{Heatmap}
#' @param name_legend Character value.Parameter for name argument of the
#' function \emph{Heatmap}
#' @return Heatmap plot produced by function \emph{Heatmap}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#' @export plot_heatmap
plot_heatmap <- function(new_classification,old_classification, dist_ang_matrix, cluster_columns = F, cluster_rows = T, name_legend)
{

  cluster_unique <- unique(new_classification)
  cluster_unique <- sort(cluster_unique, decreasing = F)
  index <- as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    index[[i]] <- which(new_classification == cluster_unique[i])
  }
  condition_unique <- as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    condition_unique[[i]] <- old_classification[index[[i]]]
  }
  cell_all_day <- colnames(dist_ang_matrix)
  cell_unique <- as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    cell_unique[[i]] <- cell_all_day[index[[i]]]
  }
  cluster_unique_2 <- as.list(rep(0, length(unique(new_classification))))
  for (i in 1:length(cluster_unique)) {
    cluster_unique_2[[i]] <- new_classification[index[[i]]]
  }
  cluster_finale <- unlist(cluster_unique_2)
  good_finale <- unlist(cell_unique)
  medium_finale <- unlist(condition_unique)



  dist_ang_matrix_plot <- dist_ang_matrix[good_finale, good_finale]

  color_cluster <- rep(0, length(unique(new_classification)))
  for (i in 1:length(color_cluster)) {
    color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
  }
  names(color_cluster) <- as.character(cluster_unique)
  color_condition <- rep(0, length(unique(old_classification)))
  for (i in 1:length(color_condition)) {
    color_condition[i] <- gg_color_hue(length(color_condition))[i]
  }
  names(color_condition) <- as.character(unique(old_classification))
  Condition_factor <- factor(old_classification)
  names(color_condition) <- as.character(levels(Condition_factor))
  data_heatmap <- data.frame(New_classification = cluster_finale, Old_classification = medium_finale)
  haCol1 <- ComplexHeatmap::HeatmapAnnotation(df = data_heatmap, col = list(New_classification = color_cluster, Old_classification = color_condition), show_legend = T)
  ht21 <- ComplexHeatmap::Heatmap(as.matrix(dist_ang_matrix_plot), cluster_rows = cluster_rows,
                                 col = circlize::colorRamp2(c(0, (max(dist_ang_matrix_plot))), c("white", "red")), name = name_legend, cluster_columns = cluster_columns,
                                 top_annotation = haCol1, row_names_gp = gpar(fontsize = 6),
                                 show_column_names = F, show_row_names = F, clustering_method_rows = "complete", clustering_method_columns = "complete")
  ComplexHeatmap::draw(ht21)
}







#' plot_distance_matrix
#' @param cluster Vector.Can be one of the two partitions returned by
#' function \emph{clustering_angular_distance} (first element of the output).
#' @inheritParams plot_heatmap
#' @return Heatmap plot produced by function \emph{Heatmap}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#' @export plot_distance_matrix
plot_distance_matrix <- function(dist_ang_matrix, cluster){

  color_cluster <- rep(0, length(unique(cluster)))
  for (i in 1:length(color_cluster)) {
    color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
  }
  names(color_cluster) <- as.character(unique(cluster))



  data_heatmap <- data.frame(cluster = cluster)

  haCol1 <- ComplexHeatmap::HeatmapAnnotation(df = data_heatmap, col = list(cluster = color_cluster), show_legend = T, show_annotation_name = F)
  ht21 <- ComplexHeatmap::Heatmap(dist_ang_matrix, cluster_rows = T,
                                 col = circlize::colorRamp2(c(0, (max(dist_ang_matrix))), c("white","red")), name = "Euclidean Distance", cluster_columns = T,
                                 top_annotation = haCol1, row_names_gp = gpar(fontsize = 6),
                                 show_column_names = F, show_row_names = F)
  ComplexHeatmap::draw(ht21)
}





#' vi_comparison
#' We compute the variation of information (VI) between the partition provided
#' by \emph{new_classification} and \emph{old_classification}. The VI between a
#' random partitions (obtained wth re-shuflle from original labels in
#' \emph{old_classification}) and \emph{old_classification} is also computed. A
#' distribution of VI values from random partitions is built. Finally, from the
#' comparison with this distribution, an empirical p value is given to the VI
#' of the unsupervised cluster analysis.
#' @param number_iter Integer value. Specify how many random partition are
#' generated (starting from re-shuffle of labels in \emph{old_classification}).
#' @inheritParams plot_heatmap
#' @return
#'
#' Empirical p value.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/mcclust/versions/1.0/topics/vi.dist}
#' @export vi_comparison
vi_comparison <- function(old_classification, new_classification, number_iter){
  random_vi <- rep(list(0), number_iter)
  for (i in 1:number_iter) {
    set.seed(i)

    random_clas <- sample(old_classification, length(old_classification))
    random_vi[[i]] <- mcclust::vi.dist(old_classification, random_clas)

  }
  final_random_vi <- unlist(random_vi)
  unsupervised_vi <- mcclust::vi.dist(old_classification, new_classification)
  empirical_pvalue <- sum(final_random_vi <= unsupervised_vi) / length(final_random_vi)
  message(paste0("Empirical p value:", empirical_pvalue))
  return(empirical_pvalue)
}














