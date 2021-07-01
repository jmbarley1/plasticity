##Data clean and setup
require(here)
require(raster)
require(readxl)
require(tidyverse)

comm<-read.csv(here('Data','acc_data_commoncontrol.csv'), stringsAsFactors = TRUE)
str(comm)
#study:                       study data was extracted from
#taxon:                       latin name of species
#phylum:                      phylum species bolongs to
#common_name:                 common name of species
#number_of_populations:       number of populations for each species
#ecosystem:                   ecosystem species inhabits (ocean, intertidal, terrestrial, freshwater)  
#dispersal_mode:              mode of dispersal for larval phase
#acclimation_time:            time researchers acclimated individuals in GENERATION TIME (ie. less than one generation, one generation, etc)
#source_population:           population that individual was obtained from, as defined by the researcher
#latitude:                    latitude of the source population
#longitude:                   longitude of the source population
#life_history_stage:          life history stage of the individual used in experimentation
#sex:                         sex of the individual used in experimentation
#acclimation_temperature_1:   lower temperature in pairwaise contrast between acclimation temperature used 
#acclimation_temperature_2:   upper temperature in pairwaise contrast between acclimation temperature used 
#measurement_level:           level of the thermal tolerance measurement (ie. CTmax is an individual level measurement, LT50 is a population level measurement)
#thermal_limit_type:          type of thermal limit measuremed (CTmax, CTmin, LD50_high, LD50_lower)
#thermal_limit_1:             thermal limit for acclimation temperature #1
#thermal_limit_2:             thermal limit for acclimation temperature #2
#thermal_limit_error_type:    type of error association with the mean thermal limit (standard error, standard deviation, confidence intervals)
#thermal_limit_error_1:       error measurement for thermal_limit_1, given by the researchers
#thermal_limit_error_2:       error measurement for thermal_limit_2, given by the researchers
#n_1:                         sample size for a given population 1
#n_2:                         sample size for a given population 2
#ni:                          total sample sized used for each population 
#max_temp:                    maximum temp. measured at the source population lat/long 
#min_temp:                    minimum temp. measured at the source population lat/long 
#mean_temp:                   mean temp. measured at the source population lat/long 
#temp_range:                  difference between max_temp and min_temp
#upper_lower:                 factor distinguishing between upper and lower thermal limits
#ARR:                         acclimation response ratio for a given comparison between two acclimation temperatures; computed by 
                              #higher ARR means more plasticity

require(here)
require(raster)
require(readxl)
require(tidyverse)

acc<- read_xlsx(here('Data','acc_data_reformatted_sorted.xlsx'))
acc$study<-factor(acc$study)
acc$ecosystem<- factor(acc$ecosystem)
acc$phylum<- factor(acc$phylum)
acc$acclimation_temperature_1<-as.numeric(as.character(acc$acclimation_temperature_1))
acc$acclimation_temperature_2<-as.numeric(as.character(acc$acclimation_temperature_2))
acc$thermal_limit_1<-as.numeric(acc$thermal_limit_1)
acc$thermal_limit_2<-as.numeric(acc$thermal_limit_2)
acc$thermal_limit_error_type<-factor(acc$thermal_limit_error_type)
acc$thermal_limit_error_1<-as.numeric(as.character(acc$thermal_limit_error_1))
acc$thermal_limit_error_2<-as.numeric(as.character(acc$thermal_limit_error_2))
acc$upper_lower<-factor(acc$upper_lower)

acc<-acc %>% 
  dplyr::select(study,                          #study data was extracted from
         taxon,                          #latin name of species
         phylum,                         #phylum species bolongs to
         common_name,                    #common name of species
         number_of_populations,          #number of populations for each species 
         ecosystem,                      #ecosystem species inhabits (ocean, intertidal, terrestrial, freshwater)  
         dispersal_mode,                 #mode of dispersal for larval phase
         acclimation_time,               #time researchers acclimated individuals in GENERATION TIME (ie. less than one generation, one generation, etc)
         source_population,              #population that individual was obtained from, as defined by the researcher
         latitude,                       #latitude of the source population
         longitude,                      #longitude of the source population
         life_history_stage,             #life history stage of the individual used in experimentation
         sex,                            #sex of the individual used in experimentation
         acclimation_temperature_1,      #lower temperature in pairwaise contrast between acclimation temperature used 
         acclimation_temperature_2,      #upper temperature in pairwaise contrast between acclimation temperature used 
         measurement_level,              #level of the thermal tolerance measurement (ie. CTmax is an individual level measurement, LT50 is a population level measurement)
         thermal_limit_type,             #type of thermal limit measuremed (CTmax, CTmin, LD50_high, LD50_lower)
         thermal_limit_1,                #thermal limit for acclimation temperature #1
         thermal_limit_2,                #thermal limit for acclimation temperature #2
         thermal_limit_error_type,       #type of error association with the mean thermal limit (standard error, standard deviation, confidence intervals)
         thermal_limit_error_1,          #error measurement for thermal_limit_1, given by the researchers
         thermal_limit_error_2,          #error measurement for thermal_limit_2, given by the researchers
         n_1,                            #sample size for a given population 1
         n_2,                            #sample size for a given population
         max_temp,                       #maximum temp. measured at the source population lat/long 
         min_temp,                       #minimum temp. measured at the source population lat/long 
         mean_temp,                      #mean temp. measured at the source population lat/long 
         temp_range,                     #difference between max_temp and min_temp
         upper_lower                     #factor distinguishing between upper and lower thermal limits
         ) %>% 
  mutate(ARR=(thermal_limit_2-thermal_limit_1)/(acclimation_temperature_2-acclimation_temperature_1)) #adding column for Acclimation Response Ratio (ARR), measurement of plasticity
  #higher ARR means more plasticity


#Other things to note: 

#checked Manis and Claussen 1986- Thermal limit values are correct 
#comm data is common control, meaning that lowest acc_temp is always acclimation_temp_1, and in cases where there are more than two 
  #acclimation temperatures used in one study, the lowest temp is used as the common control
#got rid of Bible et al. 2020- second acclimation temp was a heat shock and not a true acclimation 
#kept CTmin in this dataset
#calculated ni by hand in excel


