#' This function is utilized to determine the optimal number of clusters for k-means clustering on multi-sample omics data.
#' @name cncalc
#' @param data A time-series omics matrix.
#' @export cncalc
#' @importFrom fpc pamk
#' @importFrom ggprism theme_prism
#' @importFrom factoextra fviz_nbclust
cncalc <- function(data){
  df <- data.frame(round(t(apply(data, 1, scale)), 2)) %>%
    setNames(colnames(df))
  pamk.best <- pamk(df)
  print(pamk.best$nc)
  fviz_nbclust(df, kmeans, method = "silhouette")+
    theme_prism(border = TRUE)
}
