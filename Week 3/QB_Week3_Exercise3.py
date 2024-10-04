#!/usr/bin/env python3

import numpy
import re

file = open('biallelic.vcf', 'r') # open the file

frequency = [] # Opens up a list to store frequences we get from the loop 

AF = 'AF.txt' # Variable for the allele freq file the loop will print to
outfileAF = open(AF, 'w') # This opens up the AF file in 'write' mode so the print statement can write to the file

DP = 'DP.txt' # Variable for the read depth distritbution file the loop will print to
outfileDP = open(DP, 'w') # Opens DP file in 'write' mode for loop to put output in 

for line in file: #  Goes line by line in the file opened above
    if line.startswith('#'): # For each line starting with #, continue onto rest of loop
        continue
    fields = line.rstrip().split('\t') # From each line, remove all spaces & newline characters and split by tabs
    info = fields[7].split(';') # Info column is index  7 (column 8) from file. Split everything here by ;
   
    for row in info: # For each row in the newly selected ^, see if it begins with AF=
        if row.startswith('AF='):
            allele_freq = row.split('=')[1] # If it begins with AF=, split after = and output the first index, which is the number, as allele_freq
            print(allele_freq, file = outfileAF) # Prints allele_freq and saves to the AF.txt file 

    for row in info: # For each row in info (as above)
        if row.startswith('DP='): # If the row has 'DP='
            read_depth = row.split('=')[1] # Split after = and make read_depth = the first index, which is value
            print(read_depth, file = outfileDP) # Print and safe to DP.txt file



file.close()
