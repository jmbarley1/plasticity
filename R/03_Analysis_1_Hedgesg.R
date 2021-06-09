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
require(knitr)
require(viridis)
require(car)
source(here('R', '00_Data_setup.R'))

#things to delete maybe####
require(readxl)
acc<-read_xlsx(here('Data','acc_data_reformatted_sorted.xlsx'))
acc %>% 
  group_by(study) %>% 
  distinct(study) #30 studies

acc %>% 
  filter(thermal_limit_type=='LD50_high'|thermal_limit_type=='LD50_low') %>% 
  distinct(study)
acc %>% 
  filter(thermal_limit_type=='CTmin') %>% 
  distinct(study)

acc$acclimation_temperature_1<-factor(acc$acclimation_temperature_1)
acc$acclimation_temperature_2<-factor(acc$acclimation_temperature_2)

acc<-acc %>% 
  filter(n_1!=1|n_2!=1, upper_lower=='upper', thermal_limit_type=='CTmax')


acc %>% 
  group_by(study, acclimation_temperature_1, acclimation_temperature_2) %>% 
  distinct(acclimation_temperature_1, acclimation_temperature_2) %>% 
  print(n=92)

acc %>% 
  ggplot(aes(x=mean_temp, y=acclimation_temperature_1))+
  geom_point()+
  geom_abline(intercept=0, slope=1)

acc$acclimation_temperature_1<-as.numeric(as.character(acc$acclimation_temperature_1))
acc$acclimation_temperature_2<-as.numeric(as.character(acc$acclimation_temperature_2))

acc<- acc %>% 
  mutate(acc_anom=acclimation_temperature_1-mean_temp)

acc %>% 
  ggplot(aes(x=study, y=acc_anom, color=source_population))+
  geom_point()+
  coord_flip()

#looking at studies with two acclimation temps
studies<- acc %>% #this extracts the studies that used 2 acclimation temps
  group_by(study, acclimation_temperature_1) %>% 
  summarise(mean=mean(thermal_limit_1)) %>% 
  group_by(study) %>% 
  summarise(n_temps=n()) %>% 
  filter(n_temps==1) %>% 
  dplyr::select(study)

dat<-acc %>% #this creates a new dataframe that filters only the studies that used 2 temps
  semi_join(studies, by='study') 



#common control####

#data exploration
comm %>% 
  group_by(study, ecosystem, phylum, common_name) %>% 
  summarise(mean=mean(thermal_limit_1))

#average populations per study?
average<-comm%>% 
  group_by(study, taxon, ecosystem) %>% 
  summarise(n_pops=mean(number_of_populations))
sum(average$n_pops)

mean(average$n_pops) #average pop
sd(average$n_pops) #sd of average number of pop

#getting error estimates for all studies to be in the right format
comm<- comm %>% 
  filter(thermal_limit_error_1!='NA', n_1!=1|n_2!=1, upper_lower=='upper', thermal_limit_type=='CTmax') %>% #getting rid of 2 studies that do not have error estimate with their thermal limit means
  mutate(sd1i= case_when(thermal_limit_error_type=='CI' ~ (thermal_limit_error_1*sqrt(n_1))/1.96, 
                         thermal_limit_error_type=='std_err' ~ thermal_limit_error_1*sqrt(n_1),
                         thermal_limit_error_type=='std_dev' ~ thermal_limit_error_1), #converting error estimate for thermal_limit_error_1 to standard deviation
         sd2i= case_when(thermal_limit_error_type=='CI' ~ (thermal_limit_error_2*sqrt(n_2))/1.96,
                         thermal_limit_error_type=='std_err' ~ thermal_limit_error_2*sqrt(n_2),
                         thermal_limit_error_type=='std_dev' ~ thermal_limit_error_2))  #converting error estimate for thermal_limit_error_2 to standard deviation
#acc_anom calculated with acc_temp_2, because acc_temp_1 is the common control for studies that have more than 2 temps

#calculation of effect size- SMD (Hedges g)
comm_es<- escalc(measure='SMD', m1i= thermal_limit_2, n1i=n_2, sd1i=sd2i, m2i=thermal_limit_1, n2i=n_1, sd2i=sd1i, data=comm)
#this comes up with NAs in some rows because there are some rows that report 0 error (ie cold tolerance of a fish is 0C and there is no variation in this)

#changing ecosystem factor levels, combining ocean and intertidal into on 'marine' category
comm_es<-comm_es %>% 
  mutate(eco_2= case_when(ecosystem== 'ocean' ~ 'marine',
                        ecosystem== 'intertidal' ~ 'marine',
                        ecosystem== 'terrestrial' ~ 'terrestrial',
                        ecosystem== 'freshwater' ~ 'freshwater'))
comm_es$eco_2<-factor(comm_es$eco_2)

#chaning one vertabrata to chordata
comm_es$phylum<-recode(comm_es$phylum,"Vertebrata" ="Chordata")

#making column for difference in acclimation temp for each study
comm_es<- comm_es %>% 
  mutate(temp_diff= acclimation_temperature_2-acclimation_temperature_1)


comm_es$study<-factor(comm_es$study)

#rosenthall fail safe number
fsn(yi, vi, data=comm_es, type="Rosenthal") #very high number, BUT this is probably inflated because it is technically saying how many studies would 
                                            #have to show no plasticity for our results to change, NOT null result of difference in plasticity between populations


comm_es %>% 
  distinct(study)

#ARR graphs
comm_es %>% 
  ggplot(aes(x=temp_range, y=ARR, color=eco_2))+
  geom_point()+
  labs(x = 'Annual temperature range (C)', y = 'ARR', color = 'Ecosystem')+
  theme_classic()+
  geom_smooth(method='lm')+
  theme(text=element_text(size=24)) #decided against using these
  
  
#making function for var/cov matrix, from the Gleser & Olkin example on the Metafor website
calc.v <- function(x) {
  v <- matrix(1/x$n_1[1] + outer(x$yi, x$yi, "*")/(2*x$ni[1]), nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$vi
  v 
}


V <- bldiag(lapply(split(comm_es, comm_es$study), calc.v)) #creating var/cov matrix
lapply(split(comm_es, comm_es$study), calc.v) #chekcing
V #checking
all.equal(comm_es$vi, diag(V)) #check diagonals, making sure that diagonals are equal to variance
plot(comm_es$vi~diag(V)) #visualizing above and double checking  

#full model, with all moderators
full_mod<- rma.mv(yi, V, mods= ~temp_diff + temp_range + thermal_limit_1 + eco_2, 
                     slab = paste(study, sep = ""),
                     random = list(~1|study, ~1|phylum), 
                     data = comm_es)

summary(full_mod)
plot(comm_es$yi~comm_es$thermal_limit_1)
plot(comm_es$yi~comm_es$temp_range)

range(comm_es$thermal_limit_1)
range(comm_es$temp_range)
#need to center and scale continuous predictors because ranges are so different, estimates are hard to interpret

#center and scale thermal_limit_1 and temp_range and temp_diff
comm_es<- comm_es %>% 
  mutate(temp_diff_std=scale(temp_diff, center = TRUE, scale=TRUE), 
         temp_range_std=scale(temp_range, center = TRUE, scale=TRUE),
         limit_1_std=scale(thermal_limit_1, center = TRUE, scale=TRUE))
mean(comm_es$temp_diff_std)
sd(comm_es$temp_diff_std)

comm_es %>% 
  distinct(study)

#rerun full model with centered and scaled variables
full_mod_std<- rma.mv(yi, V, mods= ~temp_diff_std + temp_range_std + limit_1_std + factor(eco_2), 
                  slab = paste(study, sep = ""),
                  method='ML',
                  random = list(~1|phylum, ~1|study), 
                  data = comm_es)
summary(full_mod_std)

#collinearity
vif.rma(full_mod_std)

#MuMIn- for model comparison and averaging
eval(metafor:::.MuMIn)

comm_mods<-dredge(full_mod_std, trace = 2) #looking at all the possible models, weighting them by AICc
comm_mods
subset(comm_mods, delta<=2, recalc.weights=FALSE) #subsets full model list with all models that are within 2 delta AIC units away from top model

sw(comm_mods) #sum of model weights for each moderator

#log likelihood profiles for the full model
par(mfrow=c(1,1))
profile.rma.mv(full_mod_std)
profile(full_mod_std, sigma2=2)

#model averaging top models
avg<-model.avg(comm_mods, revised.var=FALSE)
summary(avg)

#extracting model averaged estimates
est<-coef(avg, complete=TRUE)
est<-as.data.frame(est)
est$se<-rbind(0.61036, 0.26733, 0.05864, 0.04781, 1.08772, 1.01927) #cannot figure out how to extract sd
est <- cbind(rownames(est), est)
rownames(est) <- NULL
colnames(est) <- c("mods","estimate","se")

#Figure 3: estimate plot
jpeg(file= here('Output','figure_3_estimate_plot.jpg'), width = 1500, height = 1128)
est %>% 
  mutate(mods=fct_recode(mods,'Ann. Temp. Range'='temp_range_std', 
         'Acc. Temp. Diff.'='temp_diff_std',
         'Mean Thermal Limit'='limit_1_std',
         'Freshwater'='intrcpt',
         'Terrestrial'='factor(eco_2)terrestrial',
         'Marine'='factor(eco_2)marine')) %>% 
  mutate(se=se*1.96) %>% 
  ggplot(aes(x=mods, y=estimate))+
  geom_point(size=6)+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=0.3, size=1.5)+
  theme_classic()+
  theme(text=element_text(size=24))+
  geom_hline(yintercept=0, linetype="dashed", size=2)+
  coord_flip()+
  xlab('Predictor')+
  ylab('Parameter Estimate')+
  theme(text=element_text(size=38), 
      axis.text.x = element_text(size=36), 
      axis.line = element_line(size = 1.5))
dev.off()

#estimates for full mod
summary(full_mod_std)
estimates<-coef(full_mod_std, complete=TRUE) #extracts estimates
estimates<-as.data.frame(estimates) #makes data frame so se can be added
estimates$se<-summary(full_mod_std)$se #adding se 

#diagnostics

funnel(full_mod_std)
funnel(full_mod_std, yaxis="vi", main="Sampling Variance")
funnel(full_mod_std, yaxis="seinv", main="Inverse Standard Error")
funnel(full_mod_std, yaxis="vinv", main="Inverse Sampling Variance")

comm_es$study[8]

comm_es$resid<-resid(full_mod_std)
comm_es$fitted <-fitted.rma(full_mod_std)
comm_es$rowid<-row.names(comm_es)

comm_es %>% 
  ggplot(aes(x=fitted, y=resid))+
  geom_text(aes(label=rowid))

hist(comm_es$temp_diff)

comm_es %>% 
  ggplot(aes(x=study, y=temp_diff))+
  geom_point()+
  coord_flip()

#testing to see if variation is still present in residual when you filter data by small temp_diffs
#comm_temp<-comm_es %>% 
 # filter(temp_diff<10)
#comm_temp$study<-factor(comm_temp$study)

#V
#V_temp <- bldiag(lapply(split(comm_temp, comm_temp$study), calc.v))
#lapply(split(comm_temp, comm_temp$study), calc.v)
#full_mod_temp<- rma.mv(yi, V_temp, mods= ~temp_diff_std + temp_range_std + limit_1_std + factor(eco_2), 
                      #slab = paste(study, sep = ""),
                      #random = list(~1|study, ~1|phylum), 
                      #data = comm_temp)
#summary(full_mod_temp)
#funnel(full_mod_temp)
#funnel(full_mod_temp, yaxis="vi", main="Sampling Variance")
#funnel(full_mod_temp, yaxis="seinv", main="Inverse Standard Error")
#funnel(full_mod_temp, yaxis="vinv", main="Inverse Sampling Variance")

#modeling the data that has smaller differences between acclimation temps (which seems to be increasing variation)
#does not seem to affect the models all that much. Above, we modeled only the data with less than 10 degrees difference 
#between acclimation temps and the models were relatively the same, similar residuals and estimates.
#Therefore, we are going to model all of the data

#going back to full_mod_std

#extracting estimates
summary(full_mod_std)
estimates<-coef(full_mod_std, complete=TRUE) #extracts estimates
estimates<-as.data.frame(estimates) #makes data frame so se can be added
estimates$se<-summary(full_mod_std)$se #adding se 


#making sure column names are right
estimates <- cbind(rownames(estimates), estimates)
rownames(estimates) <- NULL
colnames(estimates) <- c("mods","estimate","se")


#plotting estimates with standard error for full_mod_std

estimates %>% 
  filter(mods!='intrcpt') %>% 
  ggplot(aes(x=mods, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=0.2)
  
#predictions- full_mod_std
X <- cbind(-0.275,-0.02,seq(from=min(comm_es$limit_1_std),to=max(comm_es$limit_1_std),length.out = 100),0,1)
X          
pred<-predict.rma(full_mod_std, newmods=X)
pred<- as.data.frame(pred)
pred$limit_1_std<- X[,3]
plot(pred$pred~pred$limit_1_std)

ggplot(data=comm_es, aes(x=limit_1_std, y=yi, color=eco_2))+
  geom_point()
tapply(comm_es$temp_diff_std, comm_es$eco_2, mean)
tapply(comm_es$temp_range_std, comm_es$eco_2, mean)

ggplot(data=comm_es, aes(x=limit_1_std, y=yi))+
  geom_point(aes(color=eco_2))+
  geom_line(data= pred, aes(x=limit_1_std, y=pred))+
  geom_line(data=pred, aes(x=limit_1_std, y=ci.lb), lty=2)+
  geom_line(data=pred, aes(x=limit_1_std, y=ci.ub), lty=2)
#predictions based on thermal_limit_1 and terrestrial ecosystem, temp_diff and temp_range set to 0 (mean because they are centered and scaled)

#plot of influence diagnostics
inf<- cooks.distance.rma.mv(full_mod_std)
plot(inf)

#there are two that have high cook's distances
inf<- as.data.frame(inf)
which(inf$inf>20)
#Dong et al. ?? because ni is huge

#taking dong et al. out- do we get the same results?
test<- comm_es %>% 
  filter(study!="Dong_et_al_2015")
test$study<-factor(test$study)
V_test <- bldiag(lapply(split(test, test$study), calc.v))

test_mod<- rma.mv(yi, V_test, mods= ~temp_diff_std + temp_range_std + limit_1_std + factor(eco_2), 
                      slab = paste(study, sep = ""),
                      method='ML',
                      random = list(~1|study, ~1|phylum), 
                      data = test)
summary(test_mod) #same results
inf_test<-cooks.distance.rma.mv(test_mod)
plot(inf_test) #still some data points that have relatively high influence

which(inf_test>0.5)
#van_heerwarden et al. shows up as influential-- looks like this study has small yi and vi values

#taking van_heerwaarden et al.- same results?
test2<- test %>% 
  filter(study!='van_heerwaarden_et_al_2017')
test2$study<-factor(test2$study)
V_test2 <- bldiag(lapply(split(test2, test2$study), calc.v))

test_mod2<- rma.mv(yi, V_test2, mods= ~temp_diff_std + temp_range_std + limit_1_std + factor(eco_2), 
                  slab = paste(study, sep = ""),
                  method='ML',
                  random = list(~1|study, ~1|phylum), 
                  data = test2)
summary(test_mod2) #relatively the same results-- temp range become significant, but the effect is still pretty small
inf_test2<-cooks.distance.rma.mv(test_mod2)
plot(inf_test2) #one point from Diamond et al. has higher influence

profile.rma.mv(test_mod2)

profile.rma.mv(full_mod_std)
#taking out Dong et al. and van Heerwaarden et al. does not change the results or the LL profile. Therefore, we elected to leave those studies in

#top model- diagnostic plots
require(kableExtra)
kable(subset(comm_mods, delta<=4, recalc.weights=FALSE))
top<- rma.mv(yi, V, mods= ~temp_diff_std + limit_1_std, 
              slab = paste(study, sep = ""),
              random = list(~1|study, ~1|phylum), 
              data = comm_es)
summary(top)

##funnel plots
funnel(top)
funnel(top, yaxis="vi", main="Sampling Variance")
funnel(top, yaxis="seinv", main="Inverse Standard Error")
funnel(top, yaxis="vinv", main="Inverse Sampling Variance")

##LL profile
profile.rma.mv(top) #sigma2=2 plot shows log likelihood at 0 because there is not a phylum effect

##cook's distance
inf_top<- cooks.distance.rma.mv(top)
plot(inf_top, ylab='Cooks distance')
which(inf_top>0.5) #same points as in the full_mod_std

## figure 1: forest plot
jpeg(file= here('Output','forrest_plot.jpg'), width = 1500, height = 1128)
comm_es %>% 
  group_by(study) %>% 
  mutate(mean_yi= mean(yi)) %>%
  mutate(study=fct_recode(study,'Dong et al. 2015'='Dong_et_al_2015', 
                         'Kellerman et al. 2017'='van_heerwaarden_et_al_2017',
                         'Chen et al. 2001'='Chen_et_al_2001',
                         'Healy et al. 2019'='Healy_et_al_2019',
                         'Philips et al. 2015'='Philips_et_al_2015',
                         'Jensen et al. 2019'='Jensen_et_al_2019',
                         'Weldon et al. 2018'='Weldon_et_al_2018',
                         'Diamond et al. 2018'='Diamond_et_al_2018',
                         'Kelley et al. 2011'='Kelley_et_al_2011',
                         'Tepolt and Somero 2014'='Tepolt_et_al_2014',
                         'Barria et al. 2017'='Barria_et_al_2017',
                         'Underwood et al. 2012'='Underwood_et_al_2012',
                         'Yu et al. 2018'='Yu_et_al_2018',
                         'Bugg et al. 2020'='Bugg_et_al_2020',
                         'Manis and Claussen 1986'='Manis_and_Claussen_1986',
                         'Darveau et al. 2012'='Darveau_et_al_2012',
                         'Fernando et al. 2016'='Fernando_et_al_2016',
                         'Fangue et al. 2006'='Fangue_et_al_2006')) %>% 
  mutate(eco_2=fct_recode(eco_2, 'Freshwater'='freshwater',
                          'Marine'='marine',
                          'Terrestrial'='terrestrial')) %>% 
  ggplot(aes(x=reorder(study, -mean_yi), y=yi, color=eco_2))+
  #geom_point(position = position_dodge(width = 0.5))+
  geom_pointrange(aes(ymin=yi-vi, ymax=yi+vi), size=1.5, position = position_jitter(width = 0.2))+
  scale_color_viridis(begin=0, end=0.85, discrete = TRUE)+
  coord_flip()+
  geom_hline(size=1, yintercept=0, lty=2)+
  theme_classic()+
  theme(text=element_text(size=36))+
  labs(x = 'Study', y = 'Effect size (Hedges g)', color = 'Ecosystem') 
dev.off()

##estimate plot
summary(top)
estimates<-coef(top, complete=TRUE) #extracts estimates
estimates<-as.data.frame(estimates) #makes data frame so se can be added
estimates$se<-summary(top)$se #adding se 


#making sure column names are right
estimates <- cbind(rownames(estimates), estimates)
rownames(estimates) <- NULL
colnames(estimates) <- c("mods","estimate","se")

estimates %>% 
  filter(mods!='intrcpt') %>% 
  ggplot(aes(x=mods, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=0.2)

##making predictions with top model- thermal limit
X <- cbind(-0.275,seq(from=min(comm_es$limit_1_std),to=max(comm_es$limit_1_std),length.out = 500))
X          
pred<-predict.rma(top, newmods=X)
pred<- as.data.frame(pred)
pred$limit_1_std<- X[,2]
plot(pred$pred~pred$limit_1_std)
#predictions are made with the thermal limit with temp_diff held at 0 (which is the mean because it is centered and scaled)

ggplot(data=comm_es, aes(x=limit_1_std, y=yi, color=eco_2))+
  geom_point()
tapply(comm_es$temp_diff_std, comm_es$eco_2, mean)
tapply(comm_es$temp_range_std, comm_es$eco_2, mean)

#figure 2a
jpeg(file= here('Output','figure_2a.jpg'), width = 1500, height = 1128)
comm_es %>% 
  mutate(eco_2=fct_recode(eco_2, 'Freshwater'='freshwater',
                        'Marine'='marine',
                        'Terrestrial'='terrestrial')) %>% 
  ggplot(aes(x=limit_1_std, y=yi))+
  geom_point(aes(color=eco_2), size=6)+
  geom_line(data= pred, aes(x=limit_1_std, y=pred), size=1)+
  geom_line(data=pred, aes(x=limit_1_std, y=ci.lb), lty=2, size=.7)+
  geom_line(data=pred, aes(x=limit_1_std, y=ci.ub), lty=2, size=.7)+
  labs(x='Standardized Mean Thermal Limit', y='Plasticity (Hedges g)', color='Ecosystem')+
  scale_color_viridis(begin=0, end=0.85, discrete = TRUE)+
  theme_classic()+
  theme(text=element_text(size=38), 
        axis.text.x = element_text(size=36), 
        axis.line = element_line(size = 1.5))
dev.off()

#figure 2b
#predictions- temp-diff
Z<- cbind(seq(from=min(comm_es$temp_diff_std), to=max(comm_es$temp_diff_std), length.out = 500), 0)
Z
pred.z<-predict.rma(top, newmods = Z)
pred.z<- as.data.frame(pred.z)
pred.z$temp_diff_std<- Z[,1]
plot(pred.z$pred~pred.z$temp_diff_std)

jpeg(file= here('Output','figure_2b.jpg'), width = 1500, height = 1128)
comm_es %>% 
  mutate(eco_2=fct_recode(eco_2, 'Freshwater'='freshwater',
                          'Marine'='marine',
                          'Terrestrial'='terrestrial')) %>%
  ggplot(aes(x=temp_diff_std, y=yi))+
  geom_point(aes(color=eco_2), size=6)+
  geom_line(data= pred.z, aes(x=temp_diff_std, y=pred), size=1)+
  geom_line(data=pred.z, aes(x=temp_diff_std, y=ci.lb), lty=2, size=0.7)+
  geom_line(data=pred.z, aes(x=temp_diff_std, y=ci.ub), lty=2, size=0.7)+
  labs(x='Standardized Acc. Temperature Difference', y='Plasticity (Hedges g)', color='Ecosystem')+
  scale_color_viridis(begin=0, end=0.85, discrete = TRUE)+
  scale_y_continuous(limits=c(-2,20))+
  theme_classic()+
  theme(text=element_text(size=38), 
        axis.text.x = element_text(size=36), 
        axis.line = element_line(size = 1.5))
dev.off()

#figure 2c
jpeg(file= here('Output','figure_2c.jpg'), width = 1500, height = 1128)
comm_es %>% 
  mutate(eco_2=fct_recode(eco_2, 'Freshwater'='freshwater',
                          'Marine'='marine',
                          'Terrestrial'='terrestrial')) %>%
  ggplot(aes(x=temp_range_std, y=yi))+
  geom_point(aes(color=eco_2), size=6)+
  labs(x='Standardized Ann. Temperature Range', y='Plasticity (Hedges g)', color='Ecosystem')+
  scale_color_viridis(begin=0, end=0.85, discrete = TRUE)+
  scale_y_continuous(limits=c(-2,20))+
  theme_classic()+
  theme(text=element_text(size=38), 
        axis.text.x = element_text(size=36), 
        axis.line = element_line(size = 1.5))
dev.off()

#hedge's g against temp_range

comm_es %>% 
  ggplot(aes(x=temp_range_std, y=yi, color=eco_2))+
  geom_point(size=3)+
  labs(x = 'Temperature range (standardized)', y = 'Plasticity (Hedges g)', color = 'Ecosystem')+
  theme_classic()+
  theme(text=element_text(size=24))
  
##figure 2: making predictions with the model average
summary(avg)
#predictions not working because averaging models from a model selection object is not what predict() wants
X <- cbind(-0.275,-0.02,seq(from=min(comm_es$limit_1_std),to=max(comm_es$limit_1_std),length.out = 100),0,1)
X          
pred<-predict(avg, newdata=X)
pred<- as.data.frame(pred)
pred$limit_1_std<- X[,3]
