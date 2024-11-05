library(tidyverse)
library(broom)
library(DESeq2)

### Exercise 1

## Q1.1 

# Q1.1.1
# Load gene expression counts
counts_df <- read_delim("gtex_whole_blood_counts_downsample.txt")
  #view(counts_df)
metadata_df <- read_delim("gtex_metadata_downsample.txt")
  #view(metadata_df)

# Q1.1.2 
# Make gene names in column 1 the row names 
counts_df <- column_to_rownames(counts_df, var = "GENE_NAME")
  #view(counts_df)

# Q1.1.3
# Remove subject IDs
metadata_df <- column_to_rownames(metadata_df, var = "SUBJECT_ID")
  #view(metadata_df)

# Q1.1.4
# Look at first five rows of each df 
counts_df[1:5,]
metadata_df[1:5,]



## Q1.2

# Q1.2.1
# Make sure order is matched 
colnames(counts_df) == rownames(metadata_df)
table(colnames(counts_df) == rownames(metadata_df))


# Q1.2.2
# Make DESeq object named dds 
# create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ SEX + DTHHRDY + AGE)


## Q1.3

# Q1.3.1
# VST normalization - not used for downstream but used for PCA 
vsd <- vst(dds)

# Q1.3.2
# Plot principal components and save each
png("pca_sex.png", width = 2400, height = 1800, res = 300)
plotPCA(vsd, intgroup = "SEX") +
  labs(title = "Expression Variance by Sex") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("male" = "turquoise", "female" = "lightpink"))
dev.off()
# No difference in variance by sex 

png("pca_dthhrdy.png", width = 2400, height = 1800, res = 300)
plotPCA(vsd, intgroup = "DTHHRDY") +
  labs(title = "Expression Variance by Cause of Death") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set1")
dev.off()
# Difference exists in ventilator_case compared to fast_death_of_natural_causes

png("pca_age.png", width = 2400, height = 1800, res = 300)
plotPCA(vsd, intgroup = "AGE") +
  labs(title = "Expression Variance by Age") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_continuous("Age Group")
dev.off()
# No difference in variance by age


# Q1.3.3
  # PC1 explains 48% of the variance in the data. PC2 explains 7% of the variance in the data.
  # PC1 appears to be associated with cause of death and age group. 
  # PC1 separates cause of death into two groups, natural death and ventilator cases 
  # There is some association of PC1 with sex, but some is also associated with PC2. 



### Exercise 2 

## Q2.1 

# Extract normalized data, add to metadata df 
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()
vsd_df <- bind_cols(metadata_df, vsd_df)

# Apply multiple linear regression to WASH7P gene 
ml_WASH7P <- lm(data = vsd_df, formula = WASH7P ~ DTHHRDY + AGE + SEX) %>%
  summary() %>%
  tidy()
  view(ml_WASH7P)

# Repeat for gene SLC25A47
ml_SLC25A47 <- lm(data = vsd_df, formula = SLC25A47 ~ DTHHRDY + AGE + SEX) %>%
  summary() %>%
  tidy()
  view(ml_SLC25A47)

# Q2.1.1
  # WASH7P does not show significant evidence of sex-differential expression because p>0.05

# Q2.1.2
  # SLC25A47 does show evidence of sex-differential expression
  # There is higher expression in males because the value is positive
  # This is significant because p<0.05



## Q2.2

# Q2.2.1
# Apply DESeq2 for differential expression analysis across all genes
dds <- DESeq(dds)


## Q2.3 

# Q2.3.1
# Extract differential expression results for variable SEX
res_sex <- results(dds, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")
  view(res_sex)

# Q2.3.2
# Filter for rows with padj < 0.1 and removing where padj = NA 
#res_sex %>%
  #filter(!is.na(padj)) %>%
  #arrange(padj) %>%
  #filter(padj < 0.1) %>%
  #nrow()

# There are 262 genes with significant differential expression between males and females at 10% FDR


# Q2.3.3 
# Load mappings dataframe
loci_df <- read_delim("gene_locations.txt")
  view(loci_df)

# Merge locations_df with res_sex
merged_df <- left_join(res_sex, loci_df, by = "GENE_NAME") %>%
  arrange(padj)
  view(merged_df)

# Chromosome Y encodes the genes that are most strongly upregulated in males versus females
# There are more male-upregulated genes near the top of the list 
# This is because only males have a Y chromosome so any genes there are signifcantly differentially expressed

# Q2.3.4
WASH7P <- merged_df %>%
  filter(GENE_NAME == "WASH7P")
print(WASH7P)

SLC25A47 <- merged_df %>%
  filter(GENE_NAME == "SLC25A47")
print(SLC25A47)

# Yes, these results are broadly consistend with the linear regression findings in 2.1
# As in 2.1, no significance with WASH7P and significance with SLC25A47


## Q2.4 

# Q2.4.1
# Repeat as in 2.2/2.3 for death classification 
res_death <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes") %>%
  as_tibble(rownames = "GENE_NAME")
  view(res_death)

# Filter for rows with padj < 0.1 and removing where padj = NA 
res_death %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  filter(padj < 0.1) %>%
  nrow()

# There are 16069 genes with significant differential expression between death classifications at 10% FDR

# Q2.4.2 
# Yes, this makes sense
# In the PCA plot of cause of death, there is a clear division in the data along PC1
# This differential expression is not visible along PC1 in the pca plot by sex



### Exercise 3 

## Q3.1

volcano_plot <- ggplot(data = res_sex, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = (abs(log2FoldChange) > 1 & padj < 0.1))) +
  geom_text(data = res_sex %>% filter(abs(log2FoldChange) > 1 & padj < 0.1),
            aes(x = log2FoldChange, y = -log10(padj) + 5, label = GENE_NAME), size = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  scale_color_manual(values = c("darkgray", "coral")) +
  labs(title = "Differential Expression by Sex", y = expression(-log[10]("padj")), x = expression(log[2]("fold change")))
ggsave("volcano_plot.png", volcano_plot)

