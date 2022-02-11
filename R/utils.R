#' gg_color_hue
#' @noRd
#' @importFrom grDevices hcl
#' @importFrom graphics legend par
#' @importFrom stats as.dist cor cor.test fitted hclust median p.adjust var wilcox.test
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom grid gpar
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

