##Initial data look and exploration

require(here)
require(readxl)
require(tidyverse)
require(plotrix)

filter=dplyr::filter

#me trying to do this first####
#read in data
acc<-read.csv(here('Data','ARR_data.csv'))
head(acc)
str(acc)
levels(acc$study)

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

#How many studies for each phylum?
acc %>% 
  group_by(study, phylum) %>% 
  summarise(mean=mean(limit_1)) %>% 
  group_by(phylum) %>% 
  summarise(n_studies=n()) %>% 
  ggplot(aes(x=phylum, y=n_studies, fill=phylum))+
  geom_bar(stat='identity')+
  theme_classic()+
  xlab('Phylum')+
  ylab('Number of studies')+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))
  
#how many populations per study?
acc %>% 
  group_by(study, population) %>% 
  summarise(mean=mean(limit_1)) %>% 
  group_by(study) %>% 
  summarise(n_pop=n()) %>% 
  ggplot(aes(x=study, y=n_pop))+
  geom_bar(stat='identity')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  theme(text=element_text(size=23))+
  ylab('Number of populations')+
  xlab('Study')

#How many populations per taxon?
acc %>% 
  group_by(taxon, population) %>% 
  summarise(mean=mean(limit_1)) %>% 
  group_by(taxon) %>% 
  summarise(n_pop=n()) %>% 
  ggplot(aes(x=taxon, y=n_pop))+
  geom_bar(stat='identity')+
  ylab('Number of populations')+
  xlab('Species')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  theme(text=element_text(size=23))

#what studies do we have in this dataset?
str(acc)
levels(acc$study)
str(limit)
limit_study<-levels(limit$study)

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

#mean ARR across phylum
acc %>% 
  group_by(phylum) %>% 
  summarise(mean=mean(ARR), sd=std.error(ARR)) %>% 
  ggplot(aes(x=phylum, y=mean, fill=phylum))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
  ylab('Mean Acclimation Response Ratio')+
  xlab('Phylum')+
  theme_classic()+
  theme(legend.position='none')

acc %>% 
  ggplot(aes(x=phylum, y=ARR, fill=phylum))+
  geom_boxplot()+
  theme_classic()

#overall distribution of ARR
acc %>% 
  ggplot(aes(x=ARR))+
  geom_histogram(bins=40)
#pretty normal

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

#distribution of ARR for each phylum
acc %>% 
  ggplot(aes(x=ARR))+
  geom_histogram(bins=40)+
  facet_wrap(~phylum)+
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


  





