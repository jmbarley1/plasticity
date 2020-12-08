##Gunderson and Stillman-type analysis lm/glm

require(here)
require(readxl)
require(tidyverse)
require(plotrix)
require(AICcmodavg)

filter=dplyr::filter

#me trying to do this first####
#read in data
acc<-read.csv(here('Data','ARR_data.csv'))

#Hypothesis 1: Species/populations from more variable environments will have more plasticity####
#from EDA, data look pretty normal, so I will use gaussian errors to start with

mlist<-list(
null=glm(ARR~1, data=acc),
mod1=glm(ARR~temp_range, data = acc),
mod2=glm(ARR~temp_range*latitude, data=acc),
mod3=glm(ARR~temp_range*ecosystem*latitude, data=acc),
mod4=glm(ARR~temp_range*ecosystem*latitude*phylum, data=acc),
mod5=glm(ARR~temp_range+latitude, data=acc),
mod6=glm(ARR~temp_range+ecosystem+latitude, data=acc),
mod7=glm(ARR~temp_range+ecosystem+latitude+phylum, data=acc)
)
aictab(mlist)

#residuals 
mod4<-glm(ARR~temp_range*ecosystem*latitude*phylum, data=acc)
summary(mod4)
mod4.resid<-summary(mod4$residuals)
plot(mod4)
hist(mod4.resid, breaks = 10)

resid_mod<-resid(mod4)
fitted_mod<-fitted(mod4)
plot(resid_mod~fitted_mod)
#definitely seems to have some patterns

#Hypothesis 2: More thermally tolerant species/populations will have less plasticity####

mlist2<-list(
  null=glm(ARR~1, data=acc), 
  mod1=glm(ARR~limit_1, data = acc),
  mod2=glm(ARR~limit_1*phylum, data = acc),
  mod3=glm(ARR~limit_1*phylum*ecosystem, data = acc),
  mod4=glm(ARR~limit_1*phylum*ecosystem*latitude, data = acc),
  mod5=glm(ARR~limit_1+phylum, data = acc),
  mod6=glm(ARR~limit_1+phylum+ecosystem, data = acc),
  mod7=glm(ARR~limit_1+phylum+ecosystem+latitude, data = acc)
)
aictab(mlist2)

mod4b<-glm(ARR~limit_1*phylum*ecosystem*latitude, data = acc)
summary(mod4b)
mod4b.resid<-resid(mod4b)
mod4b.fitted<-fitted(mod4b)
hist(mod4b.resid)
plot(mod4b.resid~mod4b.fitted)





