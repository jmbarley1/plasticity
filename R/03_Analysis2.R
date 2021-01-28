#Hedge's g analysis- metafor package

require(here)
require(readxl)
require(tidyverse)
require(metafor)
require(plotrix)
require(MuMIn)
require(clubSandwich)
require(metaviz)
require(AICcmodavg)
require(bbmle)
require(rJava)
require(glmulti)

#read in data
source(here('R', '00_Data_setup.R'))
str(acc)

acc<- acc %>% 
  filter(n_1!=1 | n_2!=1) #remove n = 1 studies bc no variation
acc$study<-factor(acc$study) #also makes sure that levels that are relics are removed

#first, error measurements are not all standard deviation (which escalc needs)
acc<-acc %>% 
  filter(thermal_limit_error_1!='NA') %>% #getting rid of 2 studies that do not have error estimate with their thermal limit means
  mutate(sd1i= case_when(thermal_limit_error_type=='CI' ~ (thermal_limit_error_1*sqrt(n_1))/1.96, 
                         thermal_limit_error_type=='std_err' ~ thermal_limit_error_1*sqrt(n_1),
                         thermal_limit_error_type=='std_dev' ~ thermal_limit_error_1), #converting error estimate for thermal_limit_error_1 to standard deviation
         sd2i= case_when(thermal_limit_error_type=='CI' ~ (thermal_limit_error_2*sqrt(n_2))/1.96,
                         thermal_limit_error_type=='std_err' ~ thermal_limit_error_2*sqrt(n_2),
                         thermal_limit_error_type=='std_dev' ~ thermal_limit_error_2)) #converting error estimate for thermal_limit_error_2 to standard deviation

hist(acc$sd1i) #making sure that sd looks right

#In order to make Hedge's g more interpretable for our study, we made m1i, n1i, and sd1i thermal_limit_2
#becuase this way a positive Hedge's g mean more plasticity


#now, lets calculate hedge's g

acc_es<- escalc(measure='SMD', m1i= thermal_limit_2, n1i=n_2, sd1i=sd2i, m2i=thermal_limit_1, n2i=n_1, sd2i=sd1i, data=acc)



str(acc_es)
acc_es<- as.data.frame(acc_es)
str(acc_es)
which(is.na(acc_es$yi))#there are some rows with NA in the yi column because of an sd of 0 (study DID estimate error, but just didn't have any?)
#sutdies with partial NAs in yo: Fangue et al. (4 rows) and Yu et al. (4 rows)

#getting rid of rows with NA in yi for now
acc_es<- acc_es %>%   #remove NAs from 0 variation
  filter(yi!='NA') 
which(is.na(acc_es$yi)) #making sure this worked

#calculating ni, which is the total sample size for the whole study
acc_es<-acc_es %>% 
  group_by(study) %>% 
  mutate(ni= sum(n_1)+ sum(n_2)) 
  
acc_es$ni<-as.numeric(as.character(acc_es$ni))
acc_es$study<-factor(acc_es$study) #again, making sure that study is a factor with levels dropped
str(acc_es)

#now calculate variance/covariance matrix
#making function to calculate variance/covariance matrix
calc.v <- function(x) {
  v <- matrix(1/x$n_1[1] + outer(x$yi, x$yi, "*")/(2*x$ni[1]), nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$vi
  v 
}

#calculating matrix
V<- bldiag(lapply(split(acc_es, acc_es$study), calc.v)) 
all.equal(acc_es$vi, diag(V)) #check diagonals, making sure that diagonals are equal to variance
plot(acc_es$vi~diag(V)) #visualizing above and double checking  


#create variable for difference in acclimation temperature
acc_es<- acc_es %>% 
  mutate(temp_diff= acclimation_temperature_2-acclimation_temperature_1)

#starting with just the studies that used two acclimation temperatures

studies<- acc_es %>% #this extracts the studies that used 2 acclimation temps
  group_by(study, acclimation_temperature_1) %>% 
  summarise(mean=mean(thermal_limit_1)) %>% 
  group_by(study) %>% 
  summarise(n_temps=n()) %>% 
  filter(n_temps==1) %>% 
  dplyr::select(study)

str(studies)
str(acc)

dat<-acc_es %>% #this creates a new dataframe that filters only the studies that used 2 temps
  semi_join(studies, by='study') 

#have to recalculate V matrix for dat
dat$study<-factor(dat$study)
V<- bldiag(lapply(split(dat, dat$study), calc.v)) 

full_mod<- rma.mv(yi, V, mods= ~temp_diff + ~thermal_limit_1 + ~factor(phylum) + ~factor(ecosystem), 
                  slab = paste(study, sep = ""),
                  random = (~1|study), 
                  data = dat)

mods<- dredge(full_mod, trace=2) #takes a look at all combinations of models that have a relative AIC value less than 2
importance(mods) #not sure if other moderators are in this?
coef_test(full_mod, vcov = "CR2", cluster = dat$study)

#making new ecosystem variable lumping ocean and intertidal into 'marine'
dat<-dat %>% 
  mutate(eco_2= case_when(ecosystem== 'ocean' ~ 'marine',
                   ecosystem== 'intertidal' ~ 'marine',
                   ecosystem== 'terrestrial' ~ 'terrestrial',
                   ecosystem== 'freshwater' ~ 'freshwater'))
dat$eco_2<-factor(dat$eco_2)


#simple model
simple<-rma.mv(yi, V, mods=~temp_diff, 
               random = (~1|study),
               data=dat)
range(dat$temp_diff)
hist(dat$temp_diff)

simple2<- rma.mv(yi, V, mods=~upper_lower, 
                 random = (~1|study),
                 data=dat)
summary(simple2)

simple3<- rma.mv(yi, V, mods=~upper_lower + ~temp_diff, 
                 random = (~1|study),
                 data=dat)
summary(simple3)

dat %>% 
  ggplot(aes(x=temp_diff, y=yi))+
  geom_point()
dat %>% 
  ggplot(aes(x=temp_diff, y=ARR))+
  geom_point()
dat %>% 
  ggplot(aes(x=upper_lower, y=yi))+
  geom_boxplot()+
  geom_point()
summary(dat$upper_lower)

#separate CTmax
upper<- dat %>% 
  filter(upper_lower=='upper')
upper$study<-factor(upper$study)


V<- bldiag(lapply(split(upper, upper$study), calc.v)) 
str(upper)
full_mod<- rma.mv(yi, V, mods= ~ temp_diff * thermal_limit_1, 
                  random = (~1|study), 
                  data = upper,
                  method="ML")
summary(full_mod)
mods<- dredge(full_mod)
mods

ggplot(data=upper, aes(x=temp_diff,y=yi, color = thermal_limit_1))+geom_point()
ggplot(data=upper, aes(x=thermal_limit_1,y=yi, color = temp_diff))+geom_point()
ggplot(data=upper, aes(x=thermal_limit_1,y=ARR))+geom_point()
which(upper$ARR>1)
upper[37:40,]
#create3 models long hand to see what is going on
m1<- rma.mv(yi, V, mods= ~temp_diff, 
            slab = paste(study, sep = ""),
            method = 'ML',
            random = (~1|study), 
            data = upper)
m2<- rma.mv(yi, V, mods= ~thermal_limit_1, 
            slab = paste(study, sep = ""),
            method = 'ML',
            random = (~1|study), 
            data = upper)
m3<- rma.mv(yi, V, mods= ~temp_diff * thermal_limit_1, 
            slab = paste(study, sep = ""),
            method = 'ML',
            random = (~1|study), 
            data = upper)
AICctab(m1, m2, m3, weights=TRUE)



#glmulti


viz_forest(
  x=full_mod,
  group = dat[1:14, c('study','population')],
  study_labels = dat[1:14, 'study'],
  xlab = 'Hedges g'
) #doesnt work

#forrest plot
forest(
  dat$yi, 
  dat$vi, 
  annotate = FALSE, 
  slab = full_mod$slab, 
  ilab=cbind(dat$population),
  pch = 15
) #works, but populations have their own row

#ok lets try ggplot
#first, lets set up a variable that ranks sourve_population by ascending latitude
test<-dat %>% 
  group_by(source_population) %>% 
  mutate(rank= order(order(latitude, decreasing = FALSE)))

plot<-dat %>% 
  ggplot(aes(x=study, y=yi, ymin=yi-vi, ymax=yi+vi, col=source_population, fill=source_population))+
  geom_linerange(size=5, position=position_dodge(width=0.5))+
  geom_hline(yintercept = 0, lty=2)+
  geom_point(size=3, shape=21, color='white', stroke=0.5, position = position_dodge(width=0.5))+
  scale_x_discrete(name="Study") +
  scale_y_continuous(name="Hedges g", limits = c(-2,20)) +
  coord_flip()+
  theme_classic()
plot

asdfhasdfkj


#before 1/21/21####
full_mod<- rma.mv(yi, V, mods= ~temp_diff + ~thermal_limit_2 + ~temp_range + ~factor(phylum) + ~factor(ecosystem), 
                  slab = paste(study, sep = ""),
                  random = (~1|study), 
                  data = acc_es)

mods<- dredge(full_mod, trace=2) #takes a look at all combinations of models that have a relative AIC value less than 2
importance(mods) #not sure if other moderators are in this?

#Brians code
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












