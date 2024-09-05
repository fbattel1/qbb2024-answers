#!/usr/bin/env python3

#Q1

import sys

import numpy

fs = open(sys.argv[1], mode='r') #to open the file we place at position 1 and read it, mode 'r' for read 

fs.readline() #will read the first line of the code and disregard (skip one line)
fs.readline() #will read the second line of the code and disregard(skip second line)

line = fs.readline() 
fields = line.strip("\n").split("\t") #will strip newline code at end and split column header by tabs
tissues = fields[2:] #will get everything including & after the 2nd position, or the 3rd column of that first line; this is a list

gene_names = [] #a way to hold gene names
gene_IDs = [] #a way to hold gene IDs
expression = [] #a way to hold expression values

for line in fs: #for loop, for each line
    fields = line.strip("\n").split("\t") #reused from above
    gene_IDs.append(fields[0]) #gene IDs is a list, need to add to list using append function
    gene_names.append(fields[1]) #same as above for gene names
    expression.append(fields[2:]) #must use nested list to save 2+ into expression values 

fs.close()



#Q2

#converting into numpy array
tissues = numpy.array(tissues)
gene_IDs = numpy.array(gene_IDs) 
gene_names = numpy.array(gene_names)
expression = numpy.array(expression, dtype=float) #changing data type from string

print(tissues)
print(gene_IDs)
print(gene_names)
print(expression)

# Must tell numpy what type of data expression data is and not the others because we would like to run statistics on the expression data and can't do so if it's a string.



#Q3 - skip, nested for loop 



#Q4
#calculating mean expression values for first 10 genes

mean_expression = numpy.mean(expression[:10], axis=1)
print(mean_expression)



#Q5

dataset_mean = numpy.mean(expression)
dataset_median = numpy.median(expression)
print(dataset_mean)
print(dataset_median)

# As mean = 16.557814350910945 and median = 0.027107, we can infer that the data distribution is positively skewed.



#Q6

expression = expression + 1
transformed_expression = numpy.log2(expression)
transformed_mean = numpy.mean(transformed_expression)
transformed_median = numpy.median(transformed_expression)

print(transformed_mean) # This = 1.1150342022364093
print(transformed_median) # This = 0.03858718613570538 

# The transformed mean and median are much closer in value to each other than the non-transformed values are. The transformed mean is much lower than the non-transformed mean.



#Q7

transformed_expression_copy = numpy.copy(transformed_expression)
sorted_transformed_expression = numpy.sort(transformed_expression_copy, axis=1)

diff_array = sorted_transformed_expression[:,-1]-sorted_transformed_expression[:,-2]

print(diff_array) 


#Q8 

print(numpy.sum(diff_array >= 10))

