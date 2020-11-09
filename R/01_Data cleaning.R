##Initial data look and exploration

require(here)
require(readxl)
require(tidyverse)
require(plotrix)

filter=dplyr::filter

#me trying to do this first####
#read in data
data<-read_xlsx(here('Data','LA_meta_data.xlsx'))

str(data)

data<- data %>% 
  select(study:n) %>% 
  mutate_each(funs(factor), study, taxon, phylum, common_name, 
              ecosystem, dispersal_mode, life_history_stage, sex,
              measurement_level, thermal_limit_type, thermal_limit_error_type) %>% 
  mutate_each(funs(as.numeric), latitude, longitude, acclimation_temperature, thermal_limit, 
              thermal_limit_error) 

data %>% 
  ggplot(aes(x=study, y=thermal_limit, fill=acclimation_temperature))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme(axis.text.x = element_text(angle = 90))

#acc<- data %>% 
  select(study, acclimation_temperature, thermal_limit) %>% 
  #mutate_each(funs(factor),acclimation_temperature) %>% 
  group_by(study, acclimation_temperature) %>% 
  summarise(mean=mean(thermal_limit)) %>% 
  tally() %>% 
  filter(n>1)
acc$study<- as.character(acc$study)
str(acc)

test<-data %>% 
  mutate_each(funs(as.character), study) 

test1<-nest_join(data, acc, copy=FALSE)

#Matt's code- identifystudies with more than one acclimation temp####
require(ggpubr)
source(file = "R/la_meta_functions.R")
limit_data = read.csv("Data/limit_data2.csv")
limit_data = as.data.frame(limit_data)
limit_data = limit_data[,-1]

acc = limit_data[limit_data$acclimation_time %in% c("1","2+"),]

a = aggregate(acc$dev_group, by = list(acc$study, acc$taxon, acc$sex_1), unique)

for(i in 1:length(a[,1])){
  a$y[i] = length(a$x[[i]])
}

acc_multi = acc[acc$study %in% b$Group.1,]
acc_multi$dev_group = as.numeric(acc_multi$dev_group)
#acc_multi = acc_multi[-which(is.na(acc_multi$dev_group)),]
acc_multi$acc_level = NA

acc_study_list = unique(acc_multi$study)

for(i in 1:length(acc_study_list)){
  max = max(acc_multi$dev_group[acc_multi$study == acc_study_list[i]])
  min = min(acc_multi$dev_group[acc_multi$study == acc_study_list[i]])
  
  acc_multi$acc_level[acc_multi$study == acc_study_list[i] & acc_multi$dev_group == max] = "high"
  acc_multi$acc_level[acc_multi$study == acc_study_list[i] & acc_multi$dev_group == min] = "low"
}

acc_multi = acc_multi[-which(is.na(acc_multi$acc_level)),]

#For now, just read in the csv that Matt made with ARR####

acc<-read.csv(here('Data','ARR_data.csv'))
str(acc)

#How many studies are from each ecosystem?
acc %>% 
  group_by(study, ecosystem) %>% 
  summarise(mean=mean(limit_1)) %>% 
  group_by(ecosystem) %>% 
  summarise(n_studies=n()) %>% 
  ggplot(aes(x=ecosystem, y=n_studies, fill=ecosystem))+
  geom_bar(stat='identity')+
  theme_classic()+
  xlab('Ecosystem')+
  ylab('Number of studies')+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))
  

#mean ARR across ecosystem
acc %>% 
  group_by(ecosystem) %>% 
  summarise(mean=mean(ARR), sd=std.error(ARR)) %>% 
  ggplot(aes(x=ecosystem, y=mean, fill=ecosystem))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
  ylab('Mean Acclimation Response Ratio')+
  xlab('Ecosystem')+
  theme_classic()+
  theme(legend.position='none')

acc %>% 
  ggplot(aes(x=ecosystem, y=ARR, fill=ecosystem))+
  geom_boxplot()+
  theme_classic()

#distribution of ARRs for intertidal, freshwater, and ocean ecosystems
acc %>% 
  dplyr::filter(ecosystem!='terrestrial') %>%  
  ggplot(aes(x=ARR))+
  geom_histogram(bins=40)+
  facet_wrap(~ecosystem)+
  theme_classic()

##distribution of ARRs for terrestrial ecosystems
acc %>% 
  dplyr::filter(ecosystem=='terrestrial') %>% 
  ggplot(aes(x=ARR))+
  geom_histogram(bins = 40)+
  theme_classic()

#separate high and low thermal limits
acc %>% 
  filter(limit_type== 'CTmax'|limit_type=='LD50_high') %>% 
  ggplot(aes(x=ARR))+
  geom_histogram(bins=40)+
  facet_wrap(~ecosystem)+
  theme_classic()

acc %>% 
  filter(limit_type== 'CTmin'|limit_type=='LD50_low') %>% 
  ggplot(aes(x=ARR))+
  geom_histogram(bins=40)+
  facet_wrap(~ecosystem)+
  theme_classic()

#ARR across ecosystem and latitude####
acc %>% 
  ggplot(aes(x=latitude, y=ARR, color=ecosystem))+
  geom_point(position = position_jitter(width=1), alpha=0.5)


#separate out CTmax ARR from CTmin ARR
#upper limits
levels(acc$limit_type)
str(acc)

acc %>% 
  filter(limit_type== 'CTmax'|limit_type=='LD50_high') %>%
  ggplot(aes(x=latitude, y=ARR, color=ecosystem))+
  geom_point()

#lower limits
acc %>% 
  filter(limit_type== 'CTmin'|limit_type=='LD50_low') %>%
  ggplot(aes(x=latitude, y=ARR, color=ecosystem))+
  geom_point()

#ARR vs thermal limit####
#upper limits
acc %>% 
  filter(limit_1>20) %>% 
  ggplot(aes(x=limit_1, y=ARR, color=ecosystem))+
  geom_point()+
  geom_smooth(method = lm)+
  theme_classic()+
  xlab('Thermal Limit (\u00B0C)')+
  ggtitle('Upper Limits')

#lower limits
acc %>% filter(limit_1<20) %>% 
  ggplot(aes(x=limit_1, y=ARR, color=ecosystem))+
  geom_point()+
  geom_smooth(method = lm)+
  theme_classic()+
  xlab('Thermal Limit (\u00B0C)')+
  ggtitle('Lower Limits')

#facet by ecosystem
acc %>% 
  filter(limit_1>20) %>% 
  ggplot(aes(x=limit_1, y=ARR))+
  geom_point()+
  facet_wrap(~ecosystem)+
  geom_smooth(method = lm)+
  theme_classic()+
  xlab('Thermal Limit (\u00B0C)')+
  ggtitle('Upper Limits')

acc %>% 
  filter(limit_1<20) %>% 
  ggplot(aes(x=limit_1, y=ARR))+
  geom_point()+
  facet_wrap(~ecosystem)+
  geom_smooth(method = lm)+
  theme_classic()+
  xlab('Thermal Limit (\u00B0C)')+
  ggtitle('Upper Limits')






