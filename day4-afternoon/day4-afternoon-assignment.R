library(tidyverse)
library(broom)
library(ggthemes)

#Q1
#1.1

  #load and read data files
dnm <- read_csv(file = "/Users/cmdb/qbb2024-answers/day4-afternoon/aau1043_dnm.csv") 
ages <- read_csv(file = "/Users/cmdb/qbb2024-answers/day4-afternoon/aau1043_parental_age.csv")


#1.2

dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm=sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm=sum(Phase_combined == "mother", na.rm = TRUE))

dnm_summary

#1.3 
#see above under 1.1

#1.4

dnm_by_parental_age <- left_join(dnm_summary, ages, by = "Proband_id")

view(dnm_by_parental_age)

#Q2
  
#2.1

#Maternal
ggplot(data = dnm_by_parental_age, 
       mapping = aes(x = Mother_age, 
                     y = n_maternal_dnm)) +
  geom_point() +
  xlab("Mother age") +
  ylab("Number of maternal DNM")

#Paternal
ggplot(data = dnm_by_parental_age, 
       mapping = aes(x = Father_age, 
                     y = n_paternal_dnm)) +
  geom_point() +
  xlab("Father age") +
  ylab("Number of paternam DNM")


#2.2 - Maternal

maternal_dnm_model <- lm(data = dnm_by_parental_age,
   formula = n_maternal_dnm ~ 1 + Mother_age)

ggplot(data = dnm_by_parental_age, 
       mapping = aes(x = Mother_age, 
                     y = n_maternal_dnm)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("Mother age") +
  ylab("Number of maternal DNM")

summary(maternal_dnm_model)

  #2.2.1 
    #Mother age is positively correlated with number of maternal DMN 
    #The number of mutations increases with maternal age.
    #This matches the plot observed in 2.1

  #2.2.2
    #This relationship is significant because the p value is less than 2e-16.
    #This means that the probability of these results if the null hypothesis is true is incredibly small, less than 2e-16.


#2.3 - Paternal

paternal_dnm_model <- lm(data = dnm_by_parental_age,
                         formula = n_paternal_dnm ~ 1 + Father_age)

ggplot(data = dnm_by_parental_age, 
       mapping = aes(x = Father_age, 
                     y = n_paternal_dnm)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("Father age") +
  ylab("Number of paternal DNM")

summary(paternal_dnm_model)

  #2.3.1 
    #Father age is also positively correlated with number of paternal DMN 
    #The number of mutations increases with paternal age.
    #This matches the plot from 2.1. 

  #2.3.2
    #Like with maternal age and DNM, the relationship is significant because the p value is less than 2e-16.
    #Again, this means that the probability of these results if the null hypothesis is true is incredibly small, less than 2e-16.


#2.4

new_proband <- tibble(
  Father_age = 50.5)


print(predict(paternal_dnm_model, newdata = new_proband))

  # There would be 78.7 paternal DNMs for a proband with a father 50.5 years old at the time of birth.



#2.5

ggplot(data = dnm_by_parental_age) +
  geom_histogram(mapping = aes(x = n_maternal_dnm), binwidth = 1, fill = "pink", alpha = 0.5) +
  geom_histogram(mapping = aes(x = n_paternal_dnm), binwidth = 1, fill = "blue", alpha = 0.5) +
  xlab("Number of DNM") + 
  ylab("Frequency")



#2.6

t.test(dnm_by_parental_age$n_maternal_dnm, dnm_by_parental_age$n_paternal_dnm, paired = TRUE, alternative = "two.sided")

  #2.6.1
    # A paired t-test was used because we are comparing two different items. 
    # The result is statistically significant as the p-value was <2.2e-16. On average, the number of maternal DNMs is 39.2 fewer than paternal DNMs. 
    # This means that the probability of observing this difference if the null hypothesis is true is incredibly small, <2.2e-16.

