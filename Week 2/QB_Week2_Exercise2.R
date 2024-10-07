library(tidyverse)
library(ggthemes)
library(ggplot2)


df <- read.table("~/Documents/Quant bio/qbb2024-answers/Week 2/snp_counts.txt", header = FALSE)
headers <- c("MAF", "Feature", "Enrichment")
colnames(df) <- headers

#view(df)

df2 <- df %>%
  dplyr::mutate(Log2_enrichment = log2(Enrichment)) %>%
  dplyr::mutate(Feature = gsub("_chr1.bed", "", Feature)) %>% # Mutate Feature column to just have feature name
  dplyr::mutate(MAF = gsub("chr1_snps_", "", MAF)) %>%
  dplyr::mutate(MAF = gsub(".bed", "", MAF))

#view(df2)


ggplot(data = df2, aes(x = MAF, y = Enrichment, color = Feature, group = Feature)) +
  geom_point() +
  geom_line() +
  ggtitle("MAF Log2 Enrichment by Feature")
ggsave(filename = "snp_enrichments.pdf")





  
  