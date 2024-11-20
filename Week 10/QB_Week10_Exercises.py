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


# Test that keys are all 8 images 
#print(images.keys())

# Test that can see different images 
#plt.imshow(images["POLR2B_field1"]) # Images still dictionary, need to call key to get image
#plt.show()

# Print mean of a given image as a test
#print(np.mean(images["POLR2B_field1"]))



### ---- Exercise 2: Identifying Individual Cells ----
















