##ARR glmm analysis
require(tidyverse)
require(here)
require(readxl)

data<-read_xlsx(here('Data','acc_data_reformatted_sorted.xlsx'))
str(data)
data$acclimation_temperature_1<-as.numeric(data$acclimation_temperature_1)
data$acclimation_temperature_2<-as.numeric(data$acclimation_temperature_2)
data$thermal_limit_1<- as.numeric(data$thermal_limit_1)
data$thermal_limit_2<-as.numeric(data$thermal_limit_2)

#making ARR a column
data<-data %>% 
  select(-c(primary, checked_by, needs_review, envtl_data_reported, notes)) %>% 
  mutate(ARR= ((thermal_limit_2-thermal_limit_1)/(acclimation_temperature_2-acclimation_temperature_1)))

#making object with the names of the studies with more than one acclimation temp
studies<- data%>% 
  group_by(study, acclimation_temperature_1) %>% 
  summarise(mean=mean(thermal_limit_1)) %>% 
  group_by(study) %>% 
  summarise(n_acc=n()) %>% 
  mutate(n_acc2=n_acc+1) %>% 
  filter(n_acc2>2) %>% 
  select(study) 

data$acclimation_temperature_1<-factor(data$acclimation_temperature_1)
data$acclimation_temperature_2<-factor(data$acclimation_temperature_2)

test<-data %>% 
  semi_join(studies) %>% 
  group_by(study, source_population, acclimation_temperature_1, acclimation_temperature_2) %>% 
  sample_n(1) %>% 
  print(n=100)
