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
[![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/b46bea95-f240-4677-9490-558d1d2558fa)](https://private-user-images.githubusercontent.com/103125590/287997863-096a1d80-d19a-46a3-ba2b-2ce1a6602efe.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDE3NzA4NjUsIm5iZiI6MTcwMTc3MDU2NSwicGF0aCI6Ii8xMDMxMjU1OTAvMjg3OTk3ODYzLTA5NmExZDgwLWQxOWEtNDZhMy1iYTJiLTJjZTFhNjYwMmVmZS5wbmc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBSVdOSllBWDRDU1ZFSDUzQSUyRjIwMjMxMjA1JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDIzMTIwNVQxMDAyNDVaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT0yNjNkOGZhZDI2NWI1MGMwMzdmMWZmOTk1MTZjODA4ODFmNGMyZGUzZmM3NWRlYTMzNGIzYjU3MzY5OTU0ZDdkJlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZhY3Rvcl9pZD0wJmtleV9pZD0wJnJlcG9faWQ9MCJ9.Ey1izrH7xFvEGdclgpuUAVRpk_YtJ0HU5H50rVNLLf0)

[![image](https://github.com/anhuikylin/KmeansTrendAnalyzer/assets/103125590/67059930-87cd-4384-a5d5-9b7f919e9091)](https://private-user-images.githubusercontent.com/103125590/287997964-ba276975-f35a-4054-b9a3-eb80dd7c994d.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDE3NzA4NjUsIm5iZiI6MTcwMTc3MDU2NSwicGF0aCI6Ii8xMDMxMjU1OTAvMjg3OTk3OTY0LWJhMjc2OTc1LWYzNWEtNDA1NC1iOWEzLWViODBkZDdjOTk0ZC5wbmc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBSVdOSllBWDRDU1ZFSDUzQSUyRjIwMjMxMjA1JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDIzMTIwNVQxMDAyNDVaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT02OTBlMGE0N2E3NzBkZTc3MmFiMmUwZDI5NjljYzU4YTc3YjAxNWEwZmJmYjdmYzM3ODNlY2E0NGExNmQ1MmVkJlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZhY3Rvcl9pZD0wJmtleV9pZD0wJnJlcG9faWQ9MCJ9.jbeU1IJ72vgiYrZXSDx101AyVH2ZXBGXGlWnpvFDttc)


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
https://private-user-images.githubusercontent.com/103125590/287998243-64081669-9ecf-46d5-a24a-78750080c453.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDE3NzA4NjUsIm5iZiI6MTcwMTc3MDU2NSwicGF0aCI6Ii8xMDMxMjU1OTAvMjg3OTk4MjQzLTY0MDgxNjY5LTllY2YtNDZkNS1hMjRhLTc4NzUwMDgwYzQ1My5wbmc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBSVdOSllBWDRDU1ZFSDUzQSUyRjIwMjMxMjA1JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDIzMTIwNVQxMDAyNDVaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT0zZGI4ZmM1NTY5NjQ2NzllMjIwNmNjNGVkNjA2OTk3YTBiNTkyMzY2ZjJlODcxNWU1ZTQ3ZmZhZDdiMTI5OTk0JlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZhY3Rvcl9pZD0wJmtleV9pZD0wJnJlcG9faWQ9MCJ9.XZxoizmAAmijRow4DPf2MLsHjhtTwJelPJobc1I7t6c
