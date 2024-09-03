library("tidyverse")

df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
names(df)

df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+"), .before=1 )

df %>%
  group_by(SUBJECT) %>%
  summarize(num_of_samples=n()) %>%
  arrange(desc(num_of_samples))

df %>%
  group_by(SUBJECT) %>%
  summarize(num_of_samples=n()) %>%
  arrange(num_of_samples)

#Q4: 
  # SUBJECT K-562 has the most samples, 217
  # SUBJECTS GTEX-1JMI6 and GTEX-1PAR6 have the fewest samples, 1

df %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples=n()) %>%
  arrange(num_of_samples)

df %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples=n()) %>%
  arrange(desc(num_of_samples))

#Q5
  # Whole blood and skeletal muscle are the tissue types with the most samples.
        # This is likely because they are abundant compared to other tissue types and easy to collect.
  # Kidney medulla and cervix ectocervix are the tissue types with the least samples.
        # This is likely because the kidney medulla is difficult to access, and the ectocervix is extremely small.

subset(df, SUBJECT=="GTEX-NPJ8")
df_npj8 <- subset(df, SUBJECT=="GTEX-NPJ8")

df_npj8 %>%
  group_by(SMTSD) %>%
  summarize(num_of_samples=n()) %>%
  arrange(desc(num_of_samples))

view(df_npj8)

#Q6
  # For subject GTEX-NPJ8, whole blood has the most samples, 9.
  # For whole blood, the samples are processed using different sequencing techniques. 


df %>%
  group_by(SMATSSCR) %>%
  summarize(score=n()) %>%
  filter(!is.na(SMATSSCR)) %>%
  hist(data = df, 




