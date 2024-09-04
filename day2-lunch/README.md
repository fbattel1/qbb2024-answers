#day2-lunch answers

#BioMart

##Q1 
- (base) cmdb@QuantBio-20 day2-morning % cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c
    - There are 19618 protein coding genes 
    - I would like to learn more about miRNA because they may play a role in Ebola virus infection.


##Q2
- (base) cmdb@QuantBio-20 day2-morning % cut -f 1 hg38-gene-metadata-go.tsv | uniq -c | sort
    - ensembl_gene_id ENSG00000168036 has the most go_ids, 273

- (base) cmdb@QuantBio-20 day2-morning % grep "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -k 3 > gene_ENSG00000168036.tsv
- (base) cmdb@QuantBio-20 day2-morning % less -S gene_ENSG00000168036.tsv
    - This gene is likely involved in early development since it control and regulates the differentiation of many different cell types.


-----------
#GENCODE

##Q1
- (base) cmdb@QuantBio-20 day2-morning % grep "IG_._gene" genes.gtf | cut -f 1 | uniq -c | sort
    - The numbers below correspond to the number of IG genes on each chromosome listed.
        -  1 chr21
        -   6 chr16
        -  16 chr15
        -  48 chr22
        -  52 chr2
        -  91 chr14

- (base) cmdb@QuantBio-20 day2-morning % grep -e "IG_._pseudogene" -e "IG_pseudogene" genes.gtf | cut -f 1 | uniq -c | sort 
    - The numbers below correspond to the number of IG pseudogenes on each chromosome listed.
        -  1 chr1
        -  1 chr10
        -   1 chr18
        -  1 chr8
        - 5 chr9
        -   6 chr15
        -  8 chr16
        - 45 chr2
        - 48 chr22
        - 84 chr14

- There appear to be IG pseudogenes on more chromosomes than IG genes. Chromosomes 2, 22, and 14 contaiin the most IG genes and IG pseudogenes. 






