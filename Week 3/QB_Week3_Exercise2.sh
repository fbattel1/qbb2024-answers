#!/usr/bin/env bash

# Working with Illumina short-read seq
# Study: lab strain of S. cervisciae crossed with wine strain to map associations between genotype and phenotype
# Diploid offspring sporulated, colony from haploid spore sequenced 
# Data: yeast sergegants, or haploid genomes that should be mosaic of OG strains
# Mosaic = having long tracks of ancestry from one strain instead of another 

# Load seq data into library and unzip 

# Exercise 1
# Quality control. Only argument it takes is fastq. Produces HTML as output
# Open HTML in finder and will open webpage 

# Q1.1 
#wc -l A01_09.fastq


# Exercise 2

# Q2.1
#bwa index sacCer3.fa

#grep -c chr sacCer3.fa

 # There are 17 chromosomes in the yeast genome according to grep but google says 16
 # When grep "chrXVII" get error, when grep "chrXVI" get result

# Q2.2 & Q2.4

# for my_sample in *.fastq # Reads every file with this suffix 
# do
#     #echo ${my_sample}
#     sample_basename=`basename ${my_sample} .fastq` # Makes name of each file without suffix
#     #echo ${sample_basename}
#     rgroup="@RG\tID:${sample_basename}\tSM:${sample_basename}" # Creates a variable to add a header to each file
#     bwa mem -R ${rgroup} sacCer3.fa ${sample_basename}.fastq > ${sample_basename}.sam # Aligns reads of each file, given rgroup header, to ref genome and generates a sam file
#     samtools sort -@ 4 -O bam -o ${sample_basename}.bam ${sample_basename}.sam # Sorts by position in genome and converts to bam file
#     samtools index ${sample_basename}.bam # Indexes bam files for upload to IGV 
# done

# #Q2.3

# grep -c "HWI*" A01_09.sam  # Gets count of lines that start with HWI, which will exclude headers

# There are 669548 read alignments in the SAM file

# grep -c chrIII A01_09.sam

# 18196 alignments are to loci on chromosome III 


# Q2.4 within 2.2 loop

# Q2.5

    #2.4

    #2.5

# There are 3 SNPs in this window. 
# Uncertain about the third SNP becuase it only has 2x coverage while the others have 4x or 5x coverage

    #2.6

# The position is chrIV:825,834. This SNP does not fall within a gene.
