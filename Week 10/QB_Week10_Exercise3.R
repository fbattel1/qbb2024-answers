library(ggplot2)
library(tidyr)
library(dplyr)


# Read in file 
data <- read.table("gene_expression_results.txt", header = TRUE, sep = "\t")
#view(data) 

# Convert format for easier plotting
data_long <- data %>%
  gather(key = "Measure", value = "Value", Mean_NascentRNA, Mean_PCNA, Ratio)
#view(data_long)


# Plot nascentRNA
nascentRNA <- ggplot(subset(data_long, Measure == "Mean_NascentRNA"), # To specificallyt plot the mean_nascentRNA column
                     aes(x = Gene, y = Value, fill = Gene)) +
  geom_violin(trim = FALSE, show.legend = FALSE) + # Keeps violin plot tails # Hides legend since genes labeled on x axis 
  labs(title = "Nascent RNA Signal by Gene",
       x = "Gene",
       y = "Nascent RNA Signal") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotates the text 
plot(nascentRNA)
ggsave("nascentRNA.png", plot = nascentRNA)

# Plot PCNA
pcna <- ggplot(subset(data_long, Measure == "Mean_PCNA"), aes(x = Gene, y = Value, fill = Gene)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  labs(title = "PCNA Signal by Gene",
       x = "Gene",
       y = "PCNA Signal") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(pcna)
ggsave("pcna.png", plot = pcna)


# Plot log2 ratio 
log2ratio <- ggplot(subset(data_long, Measure == "Ratio"), aes(x = Gene, y = Value, fill = Gene)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  labs(title = "Log2 Ratio of Nascent RNA to PCNA by Gene",
       x = "Gene",
       y = "Log2 Ratio (Nascent RNA / PCNA)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(log2ratio)
ggsave("log2_ratio.png", plot = log2ratio)
