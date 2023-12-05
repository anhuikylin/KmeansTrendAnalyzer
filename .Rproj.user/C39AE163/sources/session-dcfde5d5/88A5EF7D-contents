#' This function is designed for computing k-means clustering and creating visualizations for multi-sample omics data.
#' @name KmeansR
#' @param data A time-series omics matrix.
#' @param centers Number of cluster.
#' @param table Logical value,If TRUE is selected, a kmeans_result.csv will be output.The content involves performing z-score normalization on each row of the data frame, along with information about different clusters.
#' @param angle The angle of rotation for x-axis labels.
#' @param box logical. Should a border be drawn around the plot?
#' @export KmeansR
#' @import magrittr
#' @import tidyverse
#' @import utils
#' @import stats
#' @import patchwork
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#' @import conflicted
#' @import magrittr
#' @importFrom grDevices rainbow
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr count
#' @importFrom dplyr right_join
#' @importFrom dplyr rename
#' @importFrom stats setNames
#' @importFrom patchwork wrap_plots
#' @importFrom reshape2 melt
#' @importFrom purrr map2
#' @importFrom conflicted conflict_prefer
#' @importFrom ggprism theme_prism
#' @examples
#' library(KmeansTrendAnalyzer)
#' library(tidyverse)
#' data(mango_FI_DAM)
#' df <- mango_FI_DAM %>%
#'   select(Index,`TA-1`:`C-3`) %>%
#'   column_to_rownames("Index") %>%
#'   mutate(TA = rowMeans(select(., 1:3)),
#'          TB = rowMeans(select(., 4:6)),
#'          TC = rowMeans(select(., 7:9)),
#'          A = rowMeans(select(., 10:12)),
#'          B = rowMeans(select(., 13:15)),
#'          C = rowMeans(select(., 16:18))) %>%
#'   select(!contains("-"))
#'cncalc(df)
#' KmeansR(df,centers = 3,angle = 0,box = TRUE)
KmeansR <- function(data,centers,table = FALSE,angle = 90,box = FALSE) {
  data_scale <- data.frame(round(t(apply(data, 1, scale)), 2))
  colnames(data_scale) <- colnames(data)
  cl <- kmeans(data_scale,centers = centers)
  data_new <- data.frame("index" = rownames(data_scale), "Cluster" = cl$cluster, data_scale) %>%
    as_tibble() %>%
    mutate("Cluster" = paste0("Cluster", Cluster))
  if (table == TRUE) {
    write.csv(data_new,"./kmeans_result.csv")
    write.csv(as.data.frame(cl$centers),"./centre_line.csv")
  }
  data_new <- data.frame("index" = rownames(data_scale), Cluster = cl$cluster, data_scale) %>%
    as_tibble() %>%
    mutate(Cluster = paste0("Cluster", Cluster)) %>%
    count(Cluster) %>%
    rename(Cluster_Count = n) %>%
    right_join(data_new, by = "Cluster") %>%
    mutate(Cluster = paste0(Cluster,":",Cluster_Count)) %>%
    select(-Cluster_Count)
  data_new = melt(data_new) %>% as_tibble()
  centers_line <- melt(cl$centers)
  centers_line <- split(centers_line, centers_line$Var1)
  # 定义绘图数据
  plot_data <- split(data_new, data_new$Cluster)
  # Generate a palette of distinct colors
  num_clusters <- length(unique(data_new$Cluster))  # 获取聚类集群的数量
  colors <- rainbow(num_clusters)  # 使用RColorBrewer生成一组颜色
  # Create a named vector of colors
  color_vector <- setNames(colors, unique(data_new$Cluster))
  # Modify the plotting code to use the color_vector
  plots <- map2(plot_data, centers_line, function(df, centers) {
    ggplot(df, aes(x = variable, y = value, group = index, color = Cluster)) +
      geom_line(show.legend = FALSE) +
      labs(x = "", y = "Standardised value") +
      labs(title = df$Cluster)+
      scale_color_manual(values = color_vector) +  # 使用自定义的颜色向量
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(panel.grid = element_blank()) +
      theme(axis.text = element_text(colour = 'black'))+
      theme(text=element_text(size=11,  family="serif"))+
      geom_line(data = centers, aes(x = Var2,
                                    y = value,
                                    group = factor(Var1)),
                col = "black", linewidth = 2)+
      theme_prism(border = as.logical(box))+
      theme(axis.text.x = element_text(angle = as.numeric(angle), vjust = 0.5))
  })
  patchwork_plot <- wrap_plots(plots)
  print(patchwork_plot)
}
