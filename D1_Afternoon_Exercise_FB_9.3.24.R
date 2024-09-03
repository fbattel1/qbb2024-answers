library(tidyverse)
library(ggthemes)

df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
view(df)

#Q1
read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")


#Q2
df
glimpse(df)


#Q3
RNA_seq <- df %>%
  filter(SMGEBTCHT == "TruSeq.v1")


#Q4

ggplot(data = RNA_seq) +
  geom_bar(mapping = aes(x = SMTSD)) +
  xlab("Tissue Type") +
  ylab("Number of samples") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Q5

ggplot(data=RNA_seq) +
  geom_histogram(bins = 30, mapping = aes(x = SMRIN)) +
  xlab("RNA Integrity Number") +
  ylab("Frequency")

  # A histogram is best for visualizing a single continuous distribution.
  # This distribution is bimodal with one peak around RIN 7 and another at RIN 10.

#Q6

ggplot(data=RNA_seq) +
  geom_violin(mapping = aes(x = SMTSD, y = SMRIN)) +
  xlab("Tissue Type") +
  ylab("RNA Integrity Number") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # Cells, including cultured fibroblasts, EBV-transformed lymphocytes, and leukemia cell line, have higher RNA integrity number.
  # This is likely because it is easier to get intact and high quality RNA from cell culture: the conditions are highly controlled and the RNA extraction protocols are simpler, compared to extracting from tissues.
  # This leads to the conclusion, supported by the same hypothesis, that kidney medulla could be considered a tissue outlier. One would expect RNA with lower integrity extracted from an autopsy.


#Q7
ggplot(data=RNA_seq) +
  geom_violin(mapping = aes(x = SMTSD, y = SMGNSDTC)) +
  xlab("Tissue Type") +
  ylab("Genes per sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # Compared to other tissues, the testes have a greater overall number of genes per sample.
  # This is likely because the testes express more genes than other organs (Xia et al., Cell, 2021)


#Q8

ggplot(data=RNA_seq, mapping = aes(x = SMTSISCH, y = SMRIN )) +
  geom_point(size = 0.5, alpha = 0.5) +
  xlab("Ischemic Time") +
  ylab("RNA Integrity Number") +
  facet_wrap("SMTSD") +
  geom_smooth(method = "lm")

  # For most tissue types, the RNA Integrity Number does not change with ischemic time.
  # For some tissue types, the RNA Integrity Number decreases with ischemic time (ex. adrenal gland, cervix, fallopian tube, bladder, heart)
  # The relationship does appear to depend on tissue type. 


#Q9

ggplot(data=RNA_seq, mapping = aes(x = SMTSISCH, y = SMRIN )) +
  geom_point(aes(color = SMATSSCR), size = 0.5, alpha = 0.5) +
  xlab("Ischemic Time") +
  ylab("RNA Integrity Number") +
  facet_wrap("SMTSD") +
  geom_smooth(method = "lm")

  # The longer the ischemic time, the lower the RIN and the higher the autolysis score. 
  # This is observable in the bladder, lung, and adrenal gland, for example.  



