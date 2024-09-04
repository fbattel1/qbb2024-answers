#day2-lunch answers

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








