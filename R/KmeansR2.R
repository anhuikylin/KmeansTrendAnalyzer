#' This function is designed for computing k-means clustering and creating visualizations for multi-sample omics data.
#' @name KmeansR
#' @param data A time-series omics matrix.
#' @param centers Number of cluster.
#' @param table Logical value,If TRUE is selected, a kmeans_result.csv will be output.The content involves performing z-score normalization on each row of the data frame, along with information about different clusters.
#' @param angle The angle of rotation for x-axis labels.
#' @param box logical. Should a border be drawn around the plot?
#' @param label_size Label Size on Bar Chart
#' @param legend.position Legend Position on Bar Chart. Specify the desired location of legends, choosing from options such as "none," "left," "right," "bottom," "top," or a two-element numeric vector.
#' @export KmeansR2
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
#' data <- mango_FI_DAM %>%
#'   select(Index,`TA-1`:`C-3`,Class = Class.I) %>%
#'   column_to_rownames("Index") %>%
#'   mutate(TA = rowMeans(select(., 1:3)),
#'          TB = rowMeans(select(., 4:6)),
#'          TC = rowMeans(select(., 7:9)),
#'          A = rowMeans(select(., 10:12)),
#'          B = rowMeans(select(., 13:15)),
#'          C = rowMeans(select(., 16:18))) %>%
#'   select(!contains("-")) %>%
#'   mutate(Colour = case_when(
#'     Class == "Amino acids and derivatives" ~ "#FF0000",
#'     Class == "Lipids" ~ "#FFFF00",
#'     Class == "Nucleotides and derivatives" ~ "#00FF00",
#'     Class == "Organic acids" ~ "#00FFFF",
#'     Class == "Phenolic acids" ~ "#0000FF",
#'     Class == "Others" ~ "#FF00FF"
#'   ))
#' df_mango <- data
#' KmeansR2(df_mango)
#' # ggsave("df_mango.pdf",width = 20,height = 10)

KmeansR2 <- function(data,centers = 6,table = TRUE,angle = 90,
                     box = FALSE,label_size = 3,
                     legend.position = "right") {
  Class <- data %>%
    rownames_to_column(var = "Index") %>%
    select(Index,Class)
  Colour <- data %>%
    rownames_to_column(var = "Index") %>%
    select(Class,Colour) %>%
    distinct(Class, .keep_all = TRUE)
  data <- data %>%
    select(-c(Class,Colour))
  data_scale <- data.frame(round(t(apply(data, 1, scale)), 2)) %>%
    mutate_all(~as.numeric(.))
  colnames(data_scale) <- colnames(data)
  cl <- kmeans(data_scale,centers = centers)###############
  data_new <- data.frame("Index" = rownames(data_scale), "Cluster" = cl$cluster, data_scale) %>%
    as_tibble() %>%
    mutate("Cluster" = paste0("Cluster", Cluster))
  data_new <- data.frame("Index" = rownames(data_scale), Cluster = cl$cluster, data_scale) %>%
    as_tibble() %>%
    mutate(Cluster = paste0("Cluster", Cluster)) %>%
    count(Cluster) %>%
    rename(Cluster_Count = n) %>%
    right_join(data_new, by = "Cluster") %>%
    mutate(Cluster = paste0(Cluster,":",Cluster_Count)) %>%
    select(-Cluster_Count)%>%
    left_join(Class,by = "Index") %>%
    relocate(Class,.after = Cluster) %>%
    relocate(Index,.before = Cluster)
  if (table == TRUE) {
    write.csv(data_new %>%
                left_join(Colour,by = "Class") %>%
                relocate(Colour,.after = Cluster),
              "./kmeans_result.csv")
    write.csv(as.data.frame(cl$centers),"./centre_line.csv")
  }
  data_new = reshape2::melt(data_new) %>% as_tibble()
  ##################
  centers_line <- reshape2::melt(cl$centers)
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
    p1 <- ggplot(df, aes(x = variable, y = value, group = Index, color = Cluster)) +
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
    ##################
    df_bar <- df %>%
      count(Class) %>%
      left_join(Colour,by = "Class") %>%
      arrange(desc(Class)) %>%
      mutate(
        Class2 = "Class",
        n = n,
        all = sum(n),
        decimal = n/all,
        percentage = paste0(round(decimal*100, 2), "%"),
        position = cumsum(decimal) - decimal/2
      ) %>%
      arrange(Class)
    # df_bar$Class <- factor(df_bar$Class,levels=df_bar$Class)
    # df_bar$Colour <- factor(df_bar$Colour,levels=df_bar$Colour)
    p2 <- ggplot(df_bar, aes(x = Class2, y = decimal, fill = Class))+
      geom_bar(stat = "identity", position = "stack")+
      geom_text(aes(label = percentage, y = position),size = label_size)+
      ggprism::theme_prism()+
      labs(x = "",y = "")+
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank())+
      scale_fill_manual(values = df_bar$Colour,
                        labels = df_bar$Class)
    (p1 | p2)+
      plot_layout(widths = c(2, 0.5),guides='collect') &
      theme(legend.position=legend.position,
            legend.key.size = unit(5, "pt"))
  })
  patchwork_plot <- wrap_plots(plots,)
  print(patchwork_plot)
  ##################
  return(plot_data)
}
