##ARR glmm analysis
require(tidyverse)
require(here)
require(readxl)
require(glmmTMB)

data<-read_xlsx(here('Data','acc_data_reformatted_sorted.xlsx'))
str(data)
data$acclimation_temperature_1<-as.numeric(data$acclimation_temperature_1)
data$acclimation_temperature_2<-as.numeric(data$acclimation_temperature_2)
data$thermal_limit_1<- as.numeric(data$thermal_limit_1)
data$thermal_limit_2<-as.numeric(data$thermal_limit_2)

#making ARR a column
data<-data %>% 
  select(-c(...1, primary, checked_by, needs_review, envtl_data_reported, notes)) %>% #cleaning data a little
  mutate(ARR= ((thermal_limit_2-thermal_limit_1)/(acclimation_temperature_2-acclimation_temperature_1)))

#are there NAs in the ARR column?
sum(is.na(data$ARR)) #2
which(is.na(data$ARR)) #One Armstrong and one Villeneuve data point
data<-data %>% 
  filter(ARR!='NA') #getting rid of those data points, they were missing one 
which(is.na(data$ARR)) #checking

#Why in the hell I did this originally I do not know. delete later
#making object with the names of the studies with more than one acclimation temp
#studies<- data%>% 
  #group_by(study, acclimation_temperature_1) %>% 
  #summarise(mean=mean(thermal_limit_1)) %>% 
  #group_by(study) %>% 
  #summarise(n_acc=n()) %>% 
  #mutate(n_acc2=n_acc+1) %>% 
  #filter(n_acc2>=2) %>% 
  #select(study) 
#studies %>% 
  #distinct(study)
#data %>% 
  #distinct(study)


data$acclimation_temperature_1<-factor(data$acclimation_temperature_1)
data$acclimation_temperature_2<-factor(data$acclimation_temperature_2)

data$acclimation_temperature_1<-as.numeric(as.character(data$acclimation_temperature_1))
data$acclimation_temperature_2<-as.numeric(as.character(data$acclimation_temperature_2))
#quick EDA

data %>% 
  ggplot(aes(x=temp_range, y=ARR))+
  geom_point()+
  geom_smooth(mmethod = 'loess')

data %>% 
  ggplot(aes(x=max_temp, y=ARR))+
  geom_point()+
  geom_smooth(mmethod = 'loess')

data %>% 
  ggplot(aes(x=mean_temp, y=ARR))+
  geom_point()+
  geom_smooth(mmethod = 'loess')

data %>% 
  filter(thermal_limit_type=='CTmin') %>% 
  ggplot(aes(x=thermal_limit_1, y=ARR))+
  geom_point()+
  geom_smooth(method='loess')

data %>% 
  filter(thermal_limit_type=='CTmax') %>% 
  ggplot(aes(x=thermal_limit_1, y=ARR))+
  geom_point()+
  geom_smooth(method = 'loess')

data %>% 
  ggplot(aes(x=ecosystem, y=ARR))+
  geom_boxplot()

data %>% 
  mutate(temp_diff=acclimation_temperature_2-acclimation_temperature_1) %>% 
  ggplot(aes(x=temp_diff, y=ARR))+
  geom_point(alpha=0.5)+
  geom_smooth(method='loess')
#this is good, temperature difference is already calculated into ARR and is standardized

data %>% 
  ggplot(aes(ARR))+
  geom_histogram()

#modeling
#gaussian errors for now because ARR is continuous and looks normal
#Independent variables: temp_range, thermal_limit_1, ecosystem
#response: ARR
#random effect: study and phylum

#null
glmmTMB(ARR~1, data=data)

  
  
  
  
  
  
  
  

