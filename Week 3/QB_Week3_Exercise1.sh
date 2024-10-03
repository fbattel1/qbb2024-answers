#!/usr/bin/env bash

# Working with Illumina short-read seq
# Study: lab strain of S. cervisciae crossed with wine strain to map associations between genotype and phenotype
# Diploid offspring sporulated, colony from haploid spore sequenced 
# Data: yeast sergegants, or haploid genomes that should be mosaic of OG strains
# Mosaic = having long tracks of ancestry from one strain instead of another 

# Load seq data into library and unzip 

# Exercise 1


# Q1.1 
awk 'NR % 4 == 2 {print length($0)}' A01_09.fastq # This will calculate and print the length of every 2nd line, which corresponds to the read
    # Awk replaces the need for a for loop, it will go line by line to perform the function, here being ^ 

# The sequencing reads are 76 base pairs long 
read_length=76


# Q1.2

num_reads=$(awk 'END {print NR / 4}' A01_09.fastq) # Similar to above, uses awk to go line by line, count the number of lines, and then divide by 4 since there are 4 lines per sequence in fastq file

# There are 669548 reads in the file 



# Q1.3

# To find the length of the yeast genome
genome_length=$(grep -v '>' sacCer3.fa | tr -d '\n' | wc -c)

#grep -v is print line that don't have pattern that starts with >, as the chr headers do in the fasta file
#tr -d '\n' deletes newline characters between chromosomes
#wc -c is word count of characters 

# Genome length = 12157105 bp 

# Coverage depth = (number of reads * read length) / genome length 

echo "${num_reads} * ${read_length} / ${genome_length}" | bc -l 

# Coverage = 4.19
# The expected depth of coveage is approximately 4x



# Q1.4
du -h A01*.fastq # Finds size of all files beginning with A01 and ending with .fastq. -h for readable format

# 122M	A01_09.fastq
# 119M	A01_11.fastq
# 129M	A01_23.fastq
# 145M	A01_24.fastq
# 110M	A01_27.fastq
# 111M	A01_31.fastq
# 146M	A01_35.fastq
# 130M	A01_39.fastq
# 149M	A01_62.fastq
# 113M	A01_63.fastq

# A01_62.fastq is the largest file with 146M
# A01_72.fastq is the smallest file with 110M


# Q1.5
fastqc A01_09.fastq

# The median is approximately 30
# A higher sequence quality number means that the base is likely not an error. With 30 being high, it means that most bp are not errors.
# There is not much variation in quality along the read.