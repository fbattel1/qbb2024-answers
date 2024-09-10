library(tidyverse)

expression <- read_tsv("qbb2024-answers/day4-lunch/dicts_expr.tsv") #read in gene expression data 

expression <- expression %>%
  mutate(Tissue_Data = paste0(Tissue, " ", GeneID)) %>% # Add new category with tissue and gene 
  mutate(Log2_Expr = log2(Expr + 1)) #Log transform data

ggplot(data = expression, 
       mapping = aes(x = Tissue_Data, y = Log2_Expr)) +
  geom_violin() +
  labs(x = "Tissue Type + Gene", y = "Log 2 Gene Expression") + 
  coord_flip() 

# Given the tissue specificity and high expression level, I am surprised at the variability in expression in each gene. 
# Expression of individual genes varies the least in pancreas genes. This may be because they are more tightly regulated in the pancreas.

