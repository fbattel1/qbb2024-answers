#Q2

library("tidyverse")

#Q3
df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
names(df)

print(df)

df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+"), .before=1 )

#Q4
df %>%
  group_by(SUBJECT) %>%
  summarize(num_of_samples=n()) %>%
  arrange(desc(num_of_samples))

df %>%
  group_by(SUBJECT) %>%
  summarize(num_of_samples=n()) %>%
  arrange(num_of_samples)

  # SUBJECTS K-562 and GTEX-NPJ8 have the most samples, 217
  # SUBJECTS GTEX-1JMI6 and GTEX-1PAR6 have the fewest samples, 1


#Q5
df %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples=n()) %>%
  arrange(num_of_samples)

df %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples=n()) %>%
  arrange(desc(num_of_samples))


  # Whole blood and skeletal muscle are the tissue types with the most samples.
        # This is likely because they are abundant compared to other tissue types and easy to collect.
  # Kidney medulla and cervix ectocervix are the tissue types with the least samples.
        # This is likely because the kidney medulla is difficult to access, and the ectocervix is extremely small.


#Q6
subset(df, SUBJECT=="GTEX-NPJ8")
df_npj8 <- subset(df, SUBJECT=="GTEX-NPJ8")

df_npj8 %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples=n()) %>%
  arrange(desc(num_of_samples))

view(df_npj8)

  # For subject GTEX-NPJ8, whole blood has the most samples, 9.
  # For whole blood, the samples are processed using different sequencing techniques. 

#Q7
meanSMAT <- df %>%
  filter(!is.na(SMATSSCR)) %>%
  group_by(SUBJECT) %>%
  summarize(mean(SMATSSCR))

view(meanSMAT) # Sorted descending and looked through table to get a sense of distribution values. 

sum(meanSMAT == 0)


  # 15 subjects have a mean SMATSSCR score of 0. 
  # A majority of subjects have a mean SMATSSCR below 1 and none have a mean SMATSSCR of 3.
  # In a report, this information could be presented in a histogram or a bar graph.



