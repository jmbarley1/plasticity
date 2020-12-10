#Data setup and cleaning

require(here)
require(tidyverse)
require(metafor)

acc<-read.csv(here('Data','ARR_data.csv'))
limit<-read.csv(here('Data','limit_data2_updated.csv'))
arr<-read.csv(here('Data','ARR_data_simple.csv'))
dat<- limit %>% 
  group_by(study, acclimation_temperature) %>% 
  summarize(mean=mean(thermal_limit)) %>% 
  group_by(study) %>% 
  tally() %>% 
  filter(n>1)

acc.limit<- semi_join(limit, dat, by='study')
str(acc.limit)
acc.limit$study<-droplevels(acc.limit$study) #getting rid of study factor levels not in plasticity story


#need to format data with mean per group (pop, acc temp, study)

#I am testing the mean difference (Hedge's G) analysis method on all data that has an sd, or the data has n>1

#calculate the mean thermal limit for each group with sd and n, get rid of data without sd for now
data<-acc.limit %>% 
  group_by(study, source_population, acclimation_temperature, thermal_limit_type) %>% 
  summarise(thermal_limit_mean=mean(thermal_limit), sd=sd(thermal_limit), n=n()) %>% 
  filter(sd!='NA') %>% 
  print(n=91) 
#now try to get rid of studies with more than two acclimation temps
data %>% 
  group_by(study, acclimation_temperature) %>% 
  summarise(mean=mean(thermal_limit_mean)) %>% 
  filter(study!='Diamond_et_al_2018'|study!='Underwood_et_al_2012') #tells me which studies have more than 2 temps

#get rid of these studies for now
data<- data %>% 
  filter(study!='Diamond_et_al_2018'|study!='Underwood_et_al_2012')

#Now I can use pivot_wider() to make dataframe into wat escalc wants

#data %>% 
  #pivot_wider(names_glue = c('temp1','temp2'),

#I haven't been able to do this yet with this dataset

#let's try with Matt's dataset
str(acc)
acc<-acc %>% 
  group_by(study, acc_temp_1, acc_temp_2) %>% 
  summarise(m1i=mean(limit_1), m2i=mean(limit_2), sd1i=sd(limit_1), sd2i=sd(limit_2), n1i=n())
#ok this seems like what I want

m1i<-acc$m1i
n<-acc$n1i
sd1i<-acc$sd1i
m2i<- acc$m2i
sd2i<- acc$sd2i

acc_es<- escalc(measure='SMD', m1i= m1i, n1i=n, sd1i=sd1i, m2i=m2i, n2i=n, sd2i=sd2i, data=acc)

#run model- Is there evidence of phenotypic plasticity in our data?

mod<- rma.mv(yi, vi, random=~1|study, data = acc_es)
summary(mod)
#surprise! yes, there is evidence for phenotypic plasticity

#lets do the above analysis on ARR data that only has two acc temps per study
str(arr)
arr<-arr %>% 
  group_by(study, acc_temp_1, acc_temp_2) %>% 
  summarise(m1i=mean(limit_1), m2i=mean(limit_2), sd1i=sd(limit_1), sd2i=sd(limit_2), n1i=n())

m1i<-arr$m1i
n<-arr$n1i
sd1i<-arr$sd1i
m2i<- arr$m2i
sd2i<- arr$sd2i

arr_es<- escalc(measure='SMD', m1i= m1i, n1i=n, sd1i=sd1i, m2i=m2i, n2i=n, sd2i=sd2i, data=arr)
mod<- rma.mv(yi, vi, random=~1|study, data = arr_es)
summary(mod)









