#!/usr/bin/env python3

import numpy
import scipy

#Q1.1 
# From README file, number of reads = 30,000 

genomesize = 1000000
readlength = 100
coverage = 3

num_reads = ((genomesize*coverage)/readlength)

print(num_reads)


#Q1.2 
## use an array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(1000000, int)

print(genome_coverage)

numberReads = int(num_reads)

for i in range(numberReads):
  startpos = numpy.random.randint(0, genomesize - readlength + 1) # Defining any possible index that can be a start position 
  endpos = startpos + readlength 
  genome_coverage[startpos:endpos] += 1


# Save file 
numpy.savetxt('coverages.txt', genome_coverage, delimiter = ',', fmt = '%d', header = "Coverage")


#Q1.5 - 10X COVERAGE 

genomesize = 1000000
readlength = 100
coverage = 10

num_reads = ((genomesize*coverage)/readlength)

genome_coverage = numpy.zeros(1000000, int)

numberReads = int(num_reads)

for i in range(numberReads):
  startpos = numpy.random.randint(0, genomesize - readlength + 1) 
  endpos = startpos + readlength 
  genome_coverage[startpos:endpos] += 1

numpy.savetxt('coverages10.txt', genome_coverage, delimiter = ',', fmt = '%d', header = "Coverage")


#Q1.6 - 30X COVERAGE 

genomesize = 1000000
readlength = 100
coverage = 30

num_reads = ((genomesize*coverage)/readlength)

genome_coverage = numpy.zeros(1000000, int)

numberReads = int(num_reads)

for i in range(numberReads):
  startpos = numpy.random.randint(0, genomesize - readlength + 1) 
  endpos = startpos + readlength 
  genome_coverage[startpos:endpos] += 1

numpy.savetxt('coverages30.txt', genome_coverage, delimiter = ',', fmt = '%d', header = "Coverage")


#Q1 rest - Proceed to R file 

#Q2.1
reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']
edges = {}
k = 3

for read in reads:
  for i in range(len(read) - k): # Through each read, first three letters, track the overlaps (edges)
    kmer1 = read[i: i+k] 
    kmer2 = read[i+1: i+1+k]
    edges.setdefault((kmer1,kmer2), 0)
    edges[(kmer1, kmer2)] += 1

#Q2.2 done in class

#Q2.3
print("digraph{") #digraph command within graphviz 

for left, right in edges:
  print(f"{left} -> {right}")
print("}")


#Q2.4 - 2.6 see README 




