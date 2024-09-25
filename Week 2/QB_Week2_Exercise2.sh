#!/usr/bin/env 

# In command line, chmod +x file name, then "bash ____" followed by file name 

# Script should loop (like nested for) through MAF file (provided in zip) and feature file (created in Exercise 1)
# Feature files are introns, exons, cCREs, other 
# Higher level for loop, loop through higher level SNP files
# Nested within, loop through 4 feature files 

touch snp_counts.txt # Create file with given name (written later)

# Create two variabls with array of strings, strings being file names 
SNP_files=("chr1_snps_0.1.bed" "chr1_snps_0.2.bed" "chr1_snps_0.3.bed" "chr1_snps_0.4.bed" "chr1_snps_0.5.bed")
feature_files=("exons_chr1.bed" "introns_chr1.bed" "cCREs_chr1.bed" "other_chr1.bed")
genome_file=("genome_chr1.bed")

for SNP in ${SNP_files[@]} # "Loop through each possible MAF value"
do
    #echo "MAF file: ${SNP}"
    
    bedtools coverage -a $genome_file -b $SNP > SNP_coverage_genome.txt # "Find SNP coverage whole chromosome"
    
    num_genome_SNPs=$(awk '{s+=$5}END{print s}' SNP_coverage_genome.txt) # "Sum SNP from coverage"
    #echo ${num_genome_SNPs}
   
    SNP_bp=$(awk '{s+=$6}END{print s}' SNP_coverage_genome.txt) # "Sum total bases from coverage"
    #echo ${SNP_bp}

    background=$(echo "$num_genome_SNPs / $SNP_bp" | bc -l)  #OG line didn't work: background=(bc -l -e "${num_genome_SNPs} / ${SNP_bp}")    # "Calculate the background"; a ratio, SNP per base
   
    #echo ${background}
 

    for feature in ${feature_files[@]} # "Loop through each feature name"
    do 
        #echo ${feature}
        
        bedtools coverage -a $feature -b $SNP > SNP_coverage.txt # "Find the SNP coverage of the current feature"
        
        num_SNPs=$(awk '{s+=$5}END{print s}' SNP_coverage.txt) # "Sum SNP from coverage"
        #echo ${num_SNPs}
        
        size_features=$(awk '{s+=$6}END{print s}' SNP_coverage.txt) # "Sum total bases from coverage"
        #echo ${size_features}
        
        ratio=$(echo "$num_SNPs / $size_features" | bc -l) # Calculate first ratio 
        #echo ${ratio}

        enrichment=$(echo "$ratio / $background" | bc -l) # "Calculate enrichment"
        #echo ${enrichment}
        
        echo -e "${SNP}\t${feature}\t${enrichment}" >> snp_counts.txt
    done 
done 
