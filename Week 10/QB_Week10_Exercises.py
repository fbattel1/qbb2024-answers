#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt # Plotting package for python 
import imageio
import plotly.express as px 
import plotly

na = np.newaxis

### ---- Exercise 1: Loading the Image Data ----

# Load directory with images and define genes, fields, and channel
image_dir = "/Users/cmdb/Documents/Quant bio/qbb2024-answers/Week 10/Images/"
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "PCNA", "nascentRNA"]

# Open an empty dictionary 
images = {} 

for gene in genes: # loops through each gene
    for field in fields: # for each gene, loops through field 
        image = np.zeros((520, 616, 3), np.uint16) # creates array with those dimensions 
        for i, channel in enumerate(channels):
            image_name = image_dir + gene + '_' + field + '_' + channel + '.tif' 
            image[:, :, i] = imageio.v3.imread(image_name)
        
        image = image.astype(np.float32) #resize
        for i in range(3):
            image[:, :, i] -= np.amin(image[:, :, i])
            image[:, :, i] /= np.amax(image[:, :, i])
        #print(np.amax(image.reshape(-1, 3), axis=0))
    
        image = (np.minimum(255, np.floor(image * 256))).astype(np.uint8) #resize
        images[f"{gene}_{field}"] = image
        
        ### ---- Exercise 2.1: Creating binary mask from DAPI channel ----
        dapi_mean = np.mean(image[:, :, 0]) #gives mean dapi channel for all 8 images
        mask = image[:, :, 0] >= dapi_mean



# Test that keys are all 8 images 
#print(images.keys())

# Test that can see different images 
#plt.imshow(images["POLR2B_field1"]) # Images still dictionary, need to call key to get image
#plt.show()


### ---- Exercise 2.2: Find labels for each image based on DAPI mask ----
# From live code

def find_labels(mask):
    # Set initial label
    l = 0
    # Create array to hold labels
    labels = np.zeros(mask.shape, np.int32)
    # Create list to keep track of label associations
    equivalence = [0]
    # Check upper-left corner
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    # For each non-zero column in row 0, check back pixel label
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If back pixel has a label, use same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    # For each non-zero row
    for x in range(1, mask.shape[0]):
        # Check left-most column, up  and up-right pixels
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If up pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If up-right pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        # For each non-zero column except last in nonzero rows, check up, up-right, up-right, up-left, left pixels
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If up pixel has label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If not up but up-right pixel has label, need to update equivalence table
                    if mask[x - 1, y - 1]:
                        # If up-left pixel has label, relabel up-right equivalence, up-left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If left pixel has label, relabel up-right equivalence, left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither up-left or left pixels are labeled, use up-right equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If not up, or up-right pixels have labels but up-left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If not up, up-right, or up-left pixels have labels but left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, add new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        # Check last pixel in row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # if up pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # if not up but up-left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # if not up or up-left but left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = np.array(equivalence)
    # Go backwards through all labels
    for i in range(1, len(equivalence))[::-1]:
        # Convert labels to the lowest value in the set associated with a single object
        labels[np.where(labels == i)] = equivalence[i]
    # Get set of unique labels
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[np.where(labels == j)] = i
    return labels

labels = find_labels(mask)

# Test to see if working 
#plt.imshow(labels)
#plt.show()


### ---- Exercise 2.3: Filter out labeled outliers based on size -----
# From live code 
sizes = np.bincount(labels.ravel())

for i in range(1, np.amax(labels)+1):
    where = np.where(labels == i)
    if sizes[i] < 100:
        labels[where] = 0


def filter_by_size(labels, minsize, maxsize):
    # Find label sizes
    sizes = np.bincount(labels.ravel())
    # Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        # If the number of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find all pixels for label
            where = np.where(labels == i)
            labels[where] = 0
    # Get set of unique labels
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[np.where(labels == j)] = i
    return labels

#labels = filter_by_size(labels, 500, 10000000)

#plt.imshow(labels == 25)
#plt.show()


### ---- Exercise 3: Find mean signal for each nucleus from PCNA and nascentRNA ----

results = [] # Opens a list to store our results

for gene in genes: # loops through genes
    for field in fields: # lopos through fields 
        image = images[f"{gene}_{field}"] # gets image name 
        dapi_mask = mask  # Use DAPI mask from ex. 2
        labels = find_labels(dapi_mask)  # Pull labels from ex. 2
        labels = filter_by_size(labels, minsize=100, maxsize=10000)  # Filters from ex. 2
        
        # Open lists for mean signal of nascentRNA and PCNA and the ratio 
        mean_nascentRNA = []
        mean_pcna = []
        ratios = []

        # Assuming each nucleus is labeled, loop through each label 
        for i in range(1, np.amax(labels) + 1):  # Skips background, which is labeled as 0
            where = np.where(labels == i) # each label, ith position 
            mean_nascentRNA.append(np.mean(image[where[0], where[1], 1]))  # Appends pixel of nucleus for nascentRNA to list created above
            mean_pcna.append(np.mean(image[where[0], where[1], 2]))  # Appends pixel of nucleus for PCNA to list created above

        for nRNA, pcna in zip(mean_nascentRNA, mean_pcna): # Loops through the lists from above
            if pcna != 0:  # Avoids dividing by 0 
                ratios.append(np.log2(nRNA / pcna)) # Calculates the ratio of nRNA to pcna mean intensity 
            else:
                ratios.append(0)  # Assigns ratio as 0 if PCNA is 0 

        for mnRNA, mPCNA, ratio in zip(mean_nascentRNA, mean_pcna, ratios): # For above
            results.append([gene, mnRNA, mPCNA, ratio]) # Append these lists to the results list opened outside the loop 


output_filename = "gene_expression_results.txt" # Makes a variable equal to the file name I want to populate 
with open(output_filename, 'w') as file: # Opens this output file in write mode 
    file.write("Gene\tMean_NascentRNA\tMean_PCNA\tRatio\n") # Writes headers to the file 
    
    for result in results: # For each part of the results list above 
        gene = result[0]  # Populates gene with the result for position 0, corresponding to gene name 
        mean_nascentRNA = result[1] # Populates spot 2 with the mean for that 
        mean_pcna = result[2] # As above for spot 3
        ratio = result[3] # As above for spot 4
        file.write(f"{gene_fied}\t{mean_nascentRNA:.4f}\t{mean_pcna:.4f}\t{ratio:.4f}\n") # Writes to file separating by newline 