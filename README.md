# KmeansTrendAnalyzer

![R](https://img.shields.io/badge/R-4.3.0-greenyellow) ![MacOS version](https://img.shields.io/badge/license-MIT-salmon)

#### Introduce
The `cncalc` function in the KmeansTrendAnalyzer package is employed for determining the optimal number of clusters, while the `KmeansR` function is utilized to analyze the relative abundance trends of genes/proteins/metabolites across different groups. This involves standardizing the relative abundance of genes/proteins/metabolites by their average values within each group and subsequently conducting K-means clustering analysis. Following trend analysis, further gene/protein/metabolite functional analysis, such as GO and KEGG enrichment analysis, can be performed on the identified clusters of interest. This process continues until relevant genes of interest are unearthed.

#### Installation tutorial

##### Github
```
devtools::install_github("anhuikylin/KmeansTrendAnalyzer")
```

#### Examples


```
library(KmeansTrendAnalyzer)
library(tidyverse)
data(mango_FI_DAM)
df <- mango_FI_DAM %>%
  select(Index,`TA-1`:`C-3`) %>%
  column_to_rownames("Index") %>%
  mutate(TA = rowMeans(select(., 1:3)),
         TB = rowMeans(select(., 4:6)),
         TC = rowMeans(select(., 7:9)),
         A = rowMeans(select(., 10:12)),
         B = rowMeans(select(., 13:15)),
         C = rowMeans(select(., 16:18))) %>%
  select(!contains("-"))
cncalc(df)
KmeansR(df,centers = 3,angle = 0,box = TRUE)
```
![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/e975c613-4c04-4c82-8f03-d5ace0c958fb)

![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/dab8bd59-33b3-49be-96c4-123643214fd9)


```
library(KmeansTrendAnalyzer)
library(tidyverse)
data <- mango_FI_DAM %>%
  select(Index,`TA-1`:`C-3`,Class = Class.I) %>%
  column_to_rownames("Index") %>%
  mutate(TA = rowMeans(select(., 1:3)),
         TB = rowMeans(select(., 4:6)),
         TC = rowMeans(select(., 7:9)),
         A = rowMeans(select(., 10:12)),
         B = rowMeans(select(., 13:15)),
         C = rowMeans(select(., 16:18))) %>%
  select(!contains("-")) %>% 
  mutate(Colour = case_when(
    Class == "Amino acids and derivatives" ~ "#FF0000",
    Class == "Lipids" ~ "#FFFF00",
    Class == "Nucleotides and derivatives" ~ "#00FF00",
    Class == "Organic acids" ~ "#00FFFF",
    Class == "Phenolic acids" ~ "#0000FF",
    Class == "Others" ~ "#FF00FF"
  ))
df_mango <- data
KmeansR2(df_mango)
# ggsave("df_mango.pdf",width = 20,height = 10)

```
![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/8cc5980d-b464-45ee-b62e-6cb9a63d9369)

