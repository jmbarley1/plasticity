#Hedge's g analysis- metafor package

require(here)
require(readxl)
require(tidyverse)
require(metafor)
require(plotrix)
require(MuMIn)

#read in data
source(here('R', '00_Data_setup.R'))
str(acc)

acc<- acc %>% 
  filter(n_1!=1 | n_2!=1) #remove n = 1 studies bc no variation
acc$study<-droplevels(acc$study)

#first, error measurements are not all standard deviation (which metafor needs)
acc<-acc %>% 
  filter(thermal_limit_error_1!='NA') %>% #getting rid of 2 studies that do not have error estimate with their thermal limit means
  mutate(sd1i= case_when(thermal_limit_error_type=='CI' ~ (thermal_limit_error_1*sqrt(n_1))/1.96, 
                         thermal_limit_error_type=='std_err' ~ thermal_limit_error_1*sqrt(n_1),
                         thermal_limit_error_type=='std_dev' ~ thermal_limit_error_1), #converting error estimate for thermal_limit_error_1 to standard deviation
         sd2i= case_when(thermal_limit_error_type=='CI' ~ (thermal_limit_error_2*sqrt(n_2))/1.96,
                         thermal_limit_error_type=='std_err' ~ thermal_limit_error_2*sqrt(n_2),
                         thermal_limit_error_type=='std_dev' ~ thermal_limit_error_2)) #converting error estimate for thermal_limit_error_2 to standard deviation

#hist(acc$sd1i)

#now, lets calculate hedge's g

acc_es<- escalc(measure='SMD', m1i= thermal_limit_2, n1i=n_2, sd1i=sd2i, m2i=thermal_limit_1, n2i=n_1, sd2i=sd1i, data=acc)
#Thermal limit mean #2 is mean 1 in escalc to make a positive hedges g mean more plasticity
acc_es<- as.data.frame(acc_es)
str(acc_es)
which(is.na(acc_es$yi))
acc_es[336:342,] #thermal limits = 0, 0 variation
acc_es[385:391,] #negative thermal limits, 0 variation
acc_es[406,]     #negative thermal limits, 0 variation
#getting rid of rows with NA in yi for now
acc_es<- acc_es %>%   #remove NAs from 0 variation
  filter(yi!='NA') 
acc$study<-factor(acc$study)

#trying to figure out data structure and NAs
#acc_es$yi<-as.numeric(as.character(acc_es$yi))
#acc_es$vi<-as.numeric(as.character(acc_es$vi))
#acc %>% 
  #group_by(study) %>% 
  #distinct(n_1) %>% 
  #print(n=56) %>% 
  #filter(n_1==1)
  
acc_es$ni = unlist(lapply(split(acc_es, acc_es$study), function(x) rep(sum(x$n_1) + x$n_2[1], each=nrow(x)))) #this assumes common control, which is not what we are doing right now

#now calculate variance/covariance matrix
calc.v <- function(x) {
  v <- matrix(1/x$n_2[1] + outer(x$yi, x$yi, "*")/(2*x$ni[1]), nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$vi
  v 
}

#making an ni colum (sum of sample size for an entire study)
acc_es<-acc_es %>% 
  group_by(study) %>% 
  mutate(ni= sum(n_1)+ sum(n_2)) 
acc_es$ni<-as.numeric(as.character(acc_es$ni))
str(acc_es)
acc_es$study<-factor(acc_es$study)
str(acc_es)

#write csv and read back in because tidy weirdness?
write.csv(acc_es,file="acc_es.csv")
acc_es_test<-read.csv(file="acc_es.csv",stringsAsFactors = T)
str(acc_es_test)
V<- bldiag(lapply(split(acc_es, acc_es$study), calc.v)) 
all.equal(acc_es$vi, diag(V)) #check diagonals
plot(acc_es$vi~diag(V))       

#trial run with smaller subset of data (acc_es)
acc_es2<-acc_es[c(1:343,351:538),]
write.csv(acc_es2,file="acc_es2.csv")
acc_small<-read.csv("acc_es2.csv", stringsAsFactors = T)
V<- bldiag(lapply(split(acc_small, acc_small$study), calc.v)) 
all.equal(acc_small$vi, diag(V)) #check diagonals
plot(acc_small$vi~diag(V))       
str(acc_small)

#just jensen
acc_es2<-acc_es[344:350,]
write.csv(acc_es2,file="acc_es2.csv")
acc_small<-read.csv("acc_es2.csv", stringsAsFactors = T)
V<- bldiag(lapply(split(acc_small, acc_small$study), calc.v)) 
all.equal(acc_small$vi, diag(V)) #check diagonals
plot(acc_small$vi~diag(V))       
str(acc_small)


#cluster robust approach
library(clubSandwich)
V2<-impute_covariance_matrix(vi = acc_es$vi, cluster = acc_es$study, r = 0.7)
V3<-impute_covariance_matrix(vi = acc_es$vi, cluster = acc_es$study, r = 0.3)

all.equal(acc_es$vi, diag(V2)) #check diagonals
length(acc_es$vi)
length(diag(V2))
plot(acc_es$vi~diag(V2))       



#create variable for difference in acclimation temperature
acc_es<- acc_es %>% 
  mutate(temp_diff= acclimation_temperature_2-acclimation_temperature_1)

m1<-rma.mv(yi, V2, mods = ~temp_diff, data=acc_es)
m2<-rma.mv(yi, V2, mods = ~temp_diff*temp_range, data=acc_es)
m3<-rma.mv(yi, V2, mods = ~temp_diff*upper_lower, data=acc_es)
m4<-rma.mv(yi, V2, mods = ~temp_diff+upper_lower, data=acc_es)

coef_test(m1, vcov = "CR2", cluster = acc_es$study)
coef_test(m2, vcov = "CR2", cluster = acc_es$study)
coef_test(m3, vcov = "CR2", cluster = acc_es$study)
coef_test(m4, vcov = "CR2", cluster = acc_es$study)

n1<-rma.mv(yi, V3, mods = ~temp_diff, data=acc_es)
n2<-rma.mv(yi, V3, mods = ~temp_diff*temp_range, data=acc_es)
n3<-rma.mv(yi, V3, mods = ~temp_diff*upper_lower, data=acc_es)
n4<-rma.mv(yi, V3, mods = ~temp_diff+upper_lower, data=acc_es)

coef_test(n1, vcov = "CR2", cluster = acc_es$study)
coef_test(n2, vcov = "CR2", cluster = acc_es$study)
coef_test(n3, vcov = "CR2", cluster = acc_es$study)
coef_test(n4, vcov = "CR2", cluster = acc_es$study)

#run models
full_mod<- rma.mv(yi, V, mods= ~temp_diff + ~thermal_limit_2 + ~temp_range + ~factor(phylum) + ~factor(ecosystem), 
                  slab = paste(study, sep = ""),
                  random = (~1|study), 
                  data = acc_es)

mods<- dredge(full_mod, trace=2) #takes a look at all combinations of models that have a relative AIC value less than 2
importance(mods) #not sure if other moderators are in this?

forest(
  acc_es$yi, 
  acc_es$vi, 
  annotate = FALSE, 
  slab = full_mod$slab, 
  pch = 15
)



asdsafasdf








#acc<-read.csv(here('Data','ARR_data.csv')) old stuff,probably going to delete ####
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









