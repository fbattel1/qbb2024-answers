library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(hexbin)
library(ggfortify)


## Q3.1

# Load file and parse
data = readr::read_tsv('salmon.merged.gene_counts.tsv') %>%
  column_to_rownames(var="gene_name")%>%
  select(-gene_id) %>%
  mutate_if(is.numeric, as.integer) %>%
  filter(rowSums(across()) > 100)
#view(data)

# Pull out narrow region samples. 21 columns
narrow = data %>% select("A1_Rep1":"P2-4_Rep3")
#view(narrow)



## Q3.2 

# Create metadata tibble - ORIGINAL
#narrow_metadata = tibble(tissue=as.factor(c("A1_Rep1", "A1_Rep2", "A1_Rep3", 
 #                                           "A2-3_Rep1", "A2-3_Rep2", "A2-3_Rep3", 
  #                                          "Cu_Rep1", "Cu_Rep2", "Cu_Rep3", 
   #                                         "LFC-Fe_Rep1", "LFC-Fe_Rep2", "LFC-Fe_Rep3", 
    #                                        "Fe_Rep1", "Fe_Rep2", "Fe_Rep3",
     #                                       "P1_Rep1", "P1_Rep2", "P1_Rep3", 
      #                                      "P2-4_Rep1", "P2-4_Rep2", "P2-4_Rep3")),
       #                 rep=as.factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)))
#view(narrow_metadata)

# Metadata tibble - CORRECTED
narrow_metadata = tibble(tissue=as.factor(c("A1", "A1", "A1", 
                                            "A2-3", "A2-3", "A2-3", 
                                            "Cu", "Cu", "Cu", 
                                            "LFC-Fe", "LFC-Fe", "LFC-Fe", 
                                            "Fe", "Fe", "Fe",
                                            "P1", "P1", "P1", 
                                            "P2-4", "P2-4", "P2-4")),
                         rep=as.factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)))
view(narrow_metadata)

# Create object
narrowdata = DESeqDataSetFromMatrix(countData=as.matrix(narrow), colData=narrow_metadata, design=~tissue)

# Plot variance by average - no vst
meanSdPlot(assay(narrowdata))

# Batch-correct 
narrowVstdata = vst(narrowdata)

# Plot data without batch-effects
meanSdPlot(assay(narrowVstdata))



## Q3.3

# PCA analysis and plot
narrowPcaData = plotPCA(narrowVstdata,intgroup=c("rep","tissue"), returnData=TRUE)

ggplot(narrowPcaData, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)
ggsave("PCA.png")



## Q3.4

# Convert to matrix
narrowVstdata = as.matrix(assay(narrowVstdata))

# Replicate means
combined = narrowVstdata[,seq(1, 21, 3)] # want only ever third sample, takes 1st replicate from each 
combined = combined + narrowVstdata[,seq(2, 21, 3)] # takes second replicate each
combined = combined + narrowVstdata[,seq(3, 21, 3)] # takes third replicate from each
combined = combined / 3 # average

# Filter out low variance genes
sds = rowSds(combined) > 1 # Where row STDEV > 1
narrowVstdata = narrowVstdata[sds,]



## Q3.5

# Set seed value
set.seed(42)

# Cluster genes 
k=kmeans(narrowVstdata, centers=12)$cluster

# Find ordering of samples 
ordering = order(k)

# Reorder genes
k = k[ordering]

# Plot heatmap of expressions and clusters but with colors instead of lines
heatmap(narrowVstdata[ordering,], Rowv=NA, Colv=NA, RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])

png("heatmap.jpg")



## Q3.6

# Pull out gene names from cluster 1
genes = rownames(narrowVstdata[k == 1,])

# Same gene names to a text file
write.table(genes, "cluster_genes.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)





