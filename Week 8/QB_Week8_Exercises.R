### Exercise 1

# Load packages

#BiocManager::install("zellkonverter")

library(zellkonverter)
library(scuttle)
library(scater)
library(scran)
library(ggplot2)


# Load and inspect data
gut <- readH5AD("v2_fca_biohub_gut_10x_raw.h5ad")
head(gut)

assayNames(gut) <- "counts"
head(gut)

gut <- logNormCounts(gut)
gut

## Q1
# 11788 genes are quantitated
# There are 13407 cells in the data set
# PCA, TSNE, and UMAp are present 

## Q2
colData(gut)
# There are 39 columns

colnames(colData(gut))
# Tissue, n_counts, and n_genes
# tissue which tells the type of tissue the sample is from
  # This can help ID clusters
# n_counts which tells the number of reads in each sample
  # This can inform which samples contribute the most to the dataset 
# n_genes which tells the number of genes in each sample
  # This can provide information about cells that contribute the most genes

plotReducedDim(gut, "X_umap", colour_by = "broad_annotation")


### Exercise 2

# Create matrix
assay(gut)

# Create vector 
genecounts <- rowSums(assay(gut))

# Summarize stats
summary(genecounts)

## Q3
# The mean is 3185 and the median is 254.
# This suggests that most of the data is low with a few higher values bringing the average up 

head(sort(genecounts, decreasing=TRUE ))
# lncRNA:Hsromega, pre-rRNA:CR4585, and lncRNA:roX1 have the highest expression
# They are all not mRNA 


## Q4a
cellcounts <- colSums(assay(gut))
summary(cellcounts)
hist(cellcounts)

# The mean is 3622
# Those with much higher counts are likely very transcriptionally active
# Perhaps they are actively growing/dividing or responding to a stimulus 


## Q4b
cellsdetected <- colSums(assay(gut)>0)
summary(cellsdetected)
hist(cellsdetected)

# The mean is 1059 
# This represents the fraction where expression is detected, removing 0 values from the statistics

## Q5

# Grep to make vector of mito gene names 
mito <- grep("^mt:", rownames(gut), value = TRUE)
head(mito)

# Create df 
df <- perCellQCMetrics(gut, subsets = list(Mito = mito))

# Convert to df and check summary
df_as_df <- as.data.frame(df)
summary(df_as_df)

# Bind metrics to cell metadata 
colData(gut) <- cbind( colData(gut), df )

# Plot mito subsets against broad_annotation 
broad_vs_mito <- plotColData(gut, x = "broad_annotation", y = "subsets_Mito_percent") + 
  theme( axis.text.x=element_text(angle=90))
#ggsave(broad_vs_mito, filename = "broad_vs_mito.png", width = 8, height = 6)

# Cells that are highly active or frequently turning over (ex. epithelial cell, gut cell) may have more mitochondrial reads
# This is because mitochondria are more active here 



### Exercise 3 

## Q6a
# New vector for cells of interst 
coi <- colData(gut)$broad_annotation == "epithelial cell"

# New single cell exp. object
epi <- gut[, coi]

# Plot epi 
epi_plot <- plotReducedDim(epi, "X_umap", colour_by = "annotation")
#ggsave(epi_plot, filename = "epi.png", width = 8, height = 6)


# Create new list with pairwise comparisons
marker.info <-scoreMarkers(epi, colData(epi)$annotation)

# ID top marker genes 
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])


## Q6b

# The six top marker genes in the anterior midgut are Mal-A6, Men-b, vnd, betaTry, Mal-A1, and Nhe2
# Mal-A6 - maltase, carbohydrate metabolic process
# Men-b - malic enzyme, glucose homeostasis 
# vnd - ventral nervous system function 
# betaTry - digestive enzyme 
# Mal-A1 - maltase, carbohydrate metabolic process
# Nhe2 - Na/H hydrogen exchanger, increases pH

# This region appears to specialize in metabolizing carbohydrates 


# Plot expression of top marker gene 
marker_genes <- plotExpression(epi, features = "Mal-A6", x = "annotation") +
  theme( axis.text.x=element_text(angle=90))
ggsave(marker_genes, filename = "marker_genes.png", width = 8, height = 10)


## Q7
# Repeat for somatic precursor cells 
coi_spc <- colData(gut)$broad_annotation == "somatic precursor cell"

# New single cell exp. object
som_pre <- gut[, coi_spc]

# Plot 
#spc_plot <- plotReducedDim(som_pre, "X_umap", colour_by = "annotation")
#ggsave(epi_plot, filename = "epi.png", width = 8, height = 6)


# Create new list with pairwise comparisons
marker.info.spc <-scoreMarkers(som_pre, colData(som_pre)$annotation)

# ID top marker genes 
chosen_spc <- marker.info.spc[["intestinal stem cell"]]
ordered_spc <- chosen_spc[order(chosen_spc$mean.AUC, decreasing=TRUE),]
head(ordered_spc[,1:4])

# Create vector with names of top six genes
goi <- rownames(ordered_spc)[1:6]

# Plot
spc_plot <- plotExpression(som_pre, features = goi, x = "annotation") +
  theme( axis.text.x=element_text(angle=90))
ggsave(spc_plot, filename = "spc_plot.png", width = 6, height = 8)


# Enteroblasts and intestinal stem cells have more similar expression 
# DI looks specific for intestinal stem cells 


