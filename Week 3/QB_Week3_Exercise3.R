library(tidyverse)
library(ggthemes)
library(ggplot2)

# Q3.2 Allele Frequency 
AF <- read.table("~/Documents/Quant bio/qbb2024-answers/Week 3/AF.txt") # Loads file into data frame
#view(AF)


ggplot(data = AF, mapping = aes(x = V1)) + # Plots data loaded as AF and places frequency values along x axis
  geom_histogram(bins= 11, color = "black", fill = "lavender") + # Histogram places frequency of frequencies on Y, colors bars
  labs(color = "", x = "Allele Frequency per Variant", y = "Frequency") + # Label axes
  theme(panel.border = element_rect("black", fill = NA, size = 0.5)) + # Outline plot
  ggtitle("Allele Frequency Distribution") + # Adds title
  ggsave(filename = "AF.jpg") # Saves

  # This is a normal distribution. Given the wide range in allele frequencies observed in the AF.txt file, the histogram looks as expected.


# Q3.3 Read Depth Distribution 
DP <- read.table("~/Documents/Quant bio/qbb2024-answers/Week 3/DP.txt", header = TRUE) # Loads file into data frame
#view(DP)

DP2 <- as.numeric(as.character(DP$V1)) # Trying everything i can to understand why histogram not plotting

sum(is.na(DP$V1)) # This is 0, so there are no "NA" values interfering 

unique(DP$V1)

ggplot(data = DP, aes(x = "V1")) + 
  geom_histogram(bins= 21, color = "black", fill = "lightblue") + 
  xlim(0, 20) +
  labs(color = "", x = "Depth Read per Variant", y = "Frequency") +
  theme(panel.border = element_rect("black", fill = NA, size = 0.5)) +
  ggtitle("Depth Read Distribution")  
 # ggsave(filename = "DP.jpg") # Saves
