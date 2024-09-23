
# QB Week Exercise 


# Q1.1 
# Download and unzip SNP files
tar xzf chr1_snps.tar.gz


# Q1.2
# Downloaded 


# Q1.3
# Downloaded


# Q1.4
# Create file without overlapping ranges 
# bedtools sort and then bedtools merge for each file 
# bedtools sorts by chromosome then coordinate by default 

sortBed -i genes.bed > sorted_genes.bed
mergeBed -i sorted_genes.bed > genes_chr1.bed 

sortBed -i exons.bed > sorted_exons.bed
mergeBed -i sorted_exons.bed > exons_chr1.bed 

sortBed -i cCREs.bed > sorted_cCREs.bed
mergeBed -i sorted_cCREs.bed > cCREs_chr1.bed 


# Q1.5
# bedtools subtract to subtract exons from full genes to get introns 

subtractBed -a genes_chr1.bed -b exons_chr1.bed > introns_chr1.bed 


# Q1.6
# subtract exons, introns, and cCRE from geneome file 

subtractBed -a genome_chr1.bed -b exons_chr1.bed introns_chr1.bed cCREs_chr1.bed > other_chr1.bed







