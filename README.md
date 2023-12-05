# 🍠KmeansTrendAnalyzer

#### 🦪Introduce
The `cncalc` function in the KmeansTrendAnalyzer package is employed for determining the optimal number of clusters, while the `KmeansR` function is utilized to analyze the relative abundance trends of genes/proteins/metabolites across different groups. This involves standardizing the relative abundance of genes/proteins/metabolites by their average values within each group and subsequently conducting K-means clustering analysis. Following trend analysis, further gene/protein/metabolite functional analysis, such as GO and KEGG enrichment analysis, can be performed on the identified clusters of interest. This process continues until relevant genes of interest are unearthed.

#### 🍣Installation tutorial

##### 🍜Github
```
devtools::install_github("anhuikylin/KmeansTrendAnalyzer")
```

#### 🦐Examples


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
![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/b46bea95-f240-4677-9490-558d1d2558fa)

![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/67059930-87cd-4384-a5d5-9b7f919e9091)
