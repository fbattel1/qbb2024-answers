library(tidyverse)
library(ggthemes)

# COVERAGE SIMULATOR - 3x
# Load file
coverages <- readr::read_tsv("~/Documents/Quant bio/qbb2024-answers/Week 1/coverages.txt")

# Group coverages and count frequency of each
coverages <- coverages %>%
  group_by(`# Coverage`) %>%
  summarize(Frequency = n()) # Code to count 0 coverage occurrences (see table)

view(coverages)

# Generate poisson distribution 
poisson_pmf <- dpois(0:14, 3) * 1000000


# Generate normal distribution   
normal_pdf <- dnorm(0:14, mean = 3, sd = 1.73) * 1000000


# Plot to show distribution, PMF, and PDF 
ggplot() +
  geom_bar(data = coverages, 
           mapping = aes(x = `# Coverage`, y = Frequency, color = "Coverage across genome"), stat = "identity", fill = "pink") +
  geom_line(mapping = aes(x = 0:14, y = poisson_pmf, color = "Poisson Distribution")) +
  geom_line(mapping = aes(x = 0:14, y = normal_pdf, color = "Normal Distribution")) +
  labs(color = "Legend", x = "Coverage", y = "Frequency") 
  ggsave(filename = "ex1_3x_cov.png")


# COVERAGE SIMULATOR - 10x

coverages10 <- readr::read_tsv("~/Documents/Quant bio/qbb2024-answers/Week 1/coverages10.txt")

coverages10 <- coverages10 %>%
  group_by(`# Coverage`) %>%
  summarize(Frequency = n()) # Code to count 0 coverage occurrences (see table)

view(coverages10)

poisson_pmf_10 <- dpois(0:26, 10) * 1000000

normal_pdf_10 <- dnorm(0:26, mean = 10, sd = 3.16) * 1000000

ggplot() +
  geom_bar(data = coverages10, 
           mapping = aes(x = `# Coverage`, y = Frequency, color = "Coverage across genome"), stat = "identity", fill = "pink") +
  geom_line(mapping = aes(x = 0:26, y = poisson_pmf_10, color = "Poisson Distribution")) +
  geom_line(mapping = aes(x = 0:26, y = normal_pdf_10, color = "Normal Distribution")) +
  labs(color = "Legend", x = "Coverage", y = "Frequency") + 
  ggsave(filename = "ex1_10x_cov.png")



# COVERAGE SIMULATOR - 30x
          
coverages30 <- readr::read_tsv("~/Documents/Quant bio/qbb2024-answers/Week 1/coverages30.txt")

coverages30 <- coverages30 %>%
  group_by(`# Coverage`) %>%
  summarize(Frequency = n()) # Code to count 0 coverage occurrences (see table)

view(coverages30)

poisson_pmf_30 <- dpois(0:55, 30) * 1000000

normal_pdf_30 <- dnorm(0:55, mean = 30, sd = 5.47) * 1000000

ggplot() +
  geom_bar(data = coverages30, 
           mapping = aes(x = `# Coverage`, y = Frequency, color = "Coverage across genome"), stat = "identity", fill = "pink") +
  geom_line(mapping = aes(x = 0:55, y = poisson_pmf_30, color = "Poisson Distribution")) +
  geom_line(mapping = aes(x = 0:55, y = normal_pdf_30, color = "Normal Distribution")) +
  labs(color = "Legend", x = "Coverage", y = "Frequency") +
  ggsave(filename = "ex1_30x_cov.png")

