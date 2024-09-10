#!/usr/bin/env python3

import sys #import packages

import numpy #import packages 

#Q1 - Load gene-tissue pairs from gene_tissue.tsv

filename = sys.argv[1] # get gene_tissue.tsv file name, the file that is the first (1, not 0) argument of the command line
fs = open(filename, mode='r') # open file in read mode 
relevant_samples = {} # creat dictionary to hold samples for gene-tissue pairs 


for line in fs: # step through file
    fields = line.rstrip("\n").split("\t") #strip newline character and split file on tabs
    key = fields[0] #gene_ID is field 0 and tissue type is field 2 in gene_tissue.tsv = create key from gene and tissue
    relevant_samples[key] = fields[2] # initialize dict from key with list to hold samples


fs.close() #close the file

#print(relevant_samples)

#Q2 - Which tissue corresponds to which sample IDs 

filename = sys.argv[2] # get metadata file name, the 2nd file typed into command line 
fs = open(filename, mode='r') # open file in read mode 
fs.readline() # to skip one line

tissue_samples = {} # creat dictionary to hold samples for tissue name 

for line in fs: # step through file
    fields = line.rstrip("\n").split("\t") #strip newline character and split on tab 
    key = fields[6] # create key for dictionary 
    value = fields[0] # create value for dictionary
    tissue_samples.setdefault(key, []) # set empty values in dictionary 
    tissue_samples[key].append(value) # initialize dictionary from key with list to hold samples, appending 

fs.close() # close the file

#print(tissue_samples) # getting one key, the tissue type, with tons and tons of values, which are the sample IDs


#Q3

filename = sys.argv[3] # get file name, third file typed into command line
fs = open(filename, mode='r') # open file in read mode 
fs.readline() # to skip one line
fs.readline() # to skip another line
header = fs.readline().rstrip("\n").split("\t") # To just read header line, and remove newline character and split at tabs  
header = header[2:] # to get header from 2nd positon onwards


#Q4

tissue_columns = {} # create a dictionary for tissue columns; to hold samples associated with tissues 

for tissue, samples in tissue_samples.items(): # Retrieve keys and values from tissue_samples dictionary 
    tissue_columns.setdefault(tissue, []) # Make new entry in the dictionary for new tissues
    for sample in samples: # Iterate through relevant samples from tissue_samples dictionary 
        if sample in header: # If sample present in gene expression data, keep track of column index 
            position = header.index(sample) #index will tell you the position # in the list of the thing you're looking for
            tissue_columns[tissue].append(position)


fs.close() # close the file

#print(tissue_columns)


#Q5

maxValue = 0
maxValueKey = ""
for tissue, samples in tissue_columns.items(): # This loop will pull out the tissue with the maximum number of samples 
    if len(samples) > maxValue: #
        maxValue = len(samples)
        maxValueKey = tissue

#print(maxValueKey) # Muscle - Skeletal


minValue = 20000000
minValueKey = ""
for tissue, samples in tissue_columns.items(): # This loop will pull out the tissue with the minimum number of samples 
    if len(samples) < minValue:
        minValue = len(samples)
        minValueKey = tissue

#print(minValueKey) # Cells - Leukemia cell line (CML)

fs.close()

#Q6 (OLD, incomplete)

f = open("test_data.gct", "r")
for l in f:
    l = l.strip().split("\t")
    
    geneName = l[0] # name of the gene
    
    if geneName in relevant_samples.keys(): # If gene is in list of relevant samples
        myTissue = relevant_samples[geneName] # Find the tissue your gene is expressing 
        #print(tissue_columns[myTissue]) # Keys of tissue_columns are tissues 


#Q7 - see R file 













