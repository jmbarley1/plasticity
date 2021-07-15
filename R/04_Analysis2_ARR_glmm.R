##ARR glmm analysis
require(tidyverse)
require(here)
require(readxl)
require(glmmTMB)
require(AICcmodavg)
require(MuMIn)

data<-read.csv(here('Data','acc_data_ARR_analysis.csv'))



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


data$acclimation_temperature_1<-factor(data$acclimation_temperature_1)
data$acclimation_temperature_2<-factor(data$acclimation_temperature_2)

data$acclimation_temperature_1<-as.numeric(as.character(data$acclimation_temperature_1))
data$acclimation_temperature_2<-as.numeric(as.character(data$acclimation_temperature_2))
#quick EDA
data<-data %>% 
  mutate(eco_2= case_when(
    ecosystem== 'marine' ~ 'marine',
    ecosystem== 'ocean' ~ 'marine',
    ecosystem== 'intertidal' ~ 'marine',
    ecosystem== 'terrestrial' ~ 'terrestrial',
    ecosystem== 'freshwater' ~ 'freshwater'))

#only using upper thermal limits
data<- data %>% 
  filter(upper_lower=='upper')

data %>% 
  ggplot(aes(x=thermal_limit_1, y=ARR))+
  geom_point()+
  geom_smooth(method = 'loess')

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
  ggplot(aes(x=eco_2, y=ARR))+
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

#models

modlist<- list(
  null= glmmTMB(ARR~ 1 + (1|study) + (1|phylum), data=data, family = 'gaussian'),
  mod1= glmmTMB(ARR~ thermal_limit_1 + (1|study) + (1|phylum),  data=data, family = 'gaussian'),
  mod2= glmmTMB(ARR ~ temp_range + (1|study) + (1|phylum),  data=data, family = 'gaussian'),
  mod3= glmmTMB(ARR ~ eco_2 + (1|study) + (1|phylum),  data=data, family = 'gaussian'),
  mod4= glmmTMB(ARR ~ thermal_limit_1 + temp_range + (1|study) + (1|phylum),  data=data, family = 'gaussian'),
  mod5= glmmTMB(ARR ~ thermal_limit_1 + eco_2 + (1|study) + (1|phylum),  data=data, family = 'gaussian'),  
  mod6= glmmTMB(ARR ~ eco_2 + temp_range + (1|study) + (1|phylum),  data=data, family = 'gaussian'),  
  mod7= glmmTMB(ARR ~ thermal_limit_1 + temp_range + eco_2 + (1|study) + (1|phylum),  data=data, family = 'gaussian'))
aictab(modlist)
AICc(modlist)

mod5= glmmTMB(ARR ~ thermal_limit_1 + eco_2 + (1|study) + (1|phylum),  data=data, family = 'gaussian')
summary(mod5) 

mod7= glmmTMB(ARR ~ thermal_limit_1 + temp_range + eco_2 + (1|study) + (1|phylum),  data=data, family = 'gaussian')
summary(mod7) 

mod1= glmmTMB(ARR~ thermal_limit_1 + (1|study) + (1|phylum),  data=data, family = 'gaussian')
summary(mod1)  
 
mod4= glmmTMB(ARR ~ thermal_limit_1 + temp_range + (1|study) + (1|phylum),  data=data, family = 'gaussian')
summary(mod4)  
 

modavg<-model.avg(mod1, mod5, mod4, mod7, rank = AICc)
summary(modavg)
