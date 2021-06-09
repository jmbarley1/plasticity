#Adding temps to data set from lon/lat
require(here)
require(readxl)
require(tidyverse)
require(metafor)
require(plotrix)

#read in data
source(here('R', '00_Data_setup.R'))


#Loads sea surface temperature data from Bio-ORACLE (Files downloaded from the BioOracle download manager: https://www.bio-oracle.org/downloads-to-email.php)
#files too large to host on github
#sst_max = raster(here('Data','BioOracle',"BO_sstmax_lonlat.tif"))
#sst_min = raster(here('Data','BioOracle','BO_sstmin_lonlat.tif'))
#sst_mean = raster(here('Data','BioOracle','BO_sstmean_lonlat.tif'))


#these tif files are too big to be hosted on the github server
#downloaded from https://chelsa-climate.org/downloads/
#Loads terrestrial data from CHELSA
#terr_mean = raster(here('Data','CHELSA','CHELSA_bio10_01.tif'))
#terr_max = raster(here('Data','CHELSA','CHELSA_bio10_05.tif'))
#terr_min = raster(here('Data','CHELSA','CHELSA_bio10_06.tif'))

#Downloads CHELSA data directly from http://chelsa-climate.org/downloads/ (Files are too large to host in GitHub Repo)
#terr_mean = raster("ftp://envicloud.wsl.ch/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_01.tif")
#terr_max = raster("ftp://envicloud.wsl.ch/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_05.tif")
#terr_min = raster("ftp://envicloud.wsl.ch/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_06.tif")


#Makes empty columns for environmental data
acc$max_temp = NA
acc$min_temp = NA
acc$mean_temp = NA
acc$temp_range = NA

#Identifies the difference ecosystems in the data set
eco_types = unique(acc$ecosystem)

#For each of the ecosystems... 
for(h in 1:length(eco_types)){
  hab = eco_types[h]
  #...record which rows in the full data set are from that ecosystem
  row_num = which(acc$ecosystem == hab)
  
  #Pull out the data for just that ecosystem
  hab_data = acc[acc$ecosystem == hab,]
  
  #Convert coordinate columns into spatial points
  cds = cbind(hab_data$longitude, hab_data$latitude)
  coords = SpatialPoints(cds)
  
  #Extract data from different environmental data sets depending on the ecosystem
  if(hab == "ocean"){
    #Using raster and bio-oracle data to extract environmental data; Buffering set to 1 km
    tmax = raster::extract(x = sst_max, y = coords, buffer = 1000, fun = mean)
    tmin = raster::extract(x = sst_min, y = coords, buffer = 1000, fun = mean)
    tmean = raster::extract(x = sst_mean, y = coords, buffer = 1000, fun = mean)
    
    #Identify which sites did not have environmental data
    missing = which(is.na(tmax))
    #If sites are missing data...
    if(length(missing) > 0){
      #...increase the buffer distance for those sites to 2 km
      tmax_missing = raster::extract(x = sst_max, y = coords[missing], buffer = 2000, fun = mean)
      
      tmin_missing = raster::extract(x = sst_min, y = coords[missing], buffer = 2000, fun = mean)
      
      tmean_missing = raster::extract(x = sst_mean, y = coords[missing], buffer = 2000, fun = mean)
      
      #Add the new data
      tmax[missing] = tmax_missing
      tmin[missing] = tmin_missing
      tmean[missing] = tmean_missing
    }
    
    #Calculate the range (difference between max and min temperatures)
    trange = as.numeric(tmax) - as.numeric(tmin)
    
  } else {
    #Using raster and CHELSA data to extract environmental data; Buffering set to 1 km
    tmax = raster::extract(x = terr_max, y = coords, buffer = 1000, fun = mean)
    tmin = raster::extract(x = terr_min, y = coords, buffer = 1000, fun = mean)
    tmean = raster::extract(x = terr_mean, y = coords, buffer = 1000, fun = mean)
    
    #If no temperature data was extracted, expand the buffering zone to 2 km
    missing = which(is.na(tmax))
    if(length(missing) > 0){
      
      tmax_missing = raster::extract(x = terr_max, y = coords[missing], buffer = 2000, fun = mean)
      tmin_missing = raster::extract(x = terr_min, y = coords[missing], buffer = 2000, fun = mean)
      tmean_missing = raster::extract(x = terr_mean, y = coords[missing], buffer = 2000, fun = mean)
      
      tmax[missing] = tmax_missing
      tmin[missing] = tmin_missing
      tmean[missing] = tmean_missing
    }
    
    #Correct the temperatures and calculate the range
    tmax = tmax / 10
    tmin = tmin / 10 
    tmean = tmean / 10
    
    trange = as.numeric(tmax) - as.numeric(tmin)
    
  }
  
  #Make sure the environmental data is numeric 
  acc$max_temp[row_num] = as.numeric(unlist(tmax))
  acc$min_temp[row_num] = as.numeric(unlist(tmin))
  acc$mean_temp[row_num] = as.numeric(unlist(tmean))
  acc$temp_range[row_num] = as.numeric(unlist(trange))
  
  still_missing = which(is.na(tmax))
  
  #print useful information as the entries are processed
  if(length(still_missing) == 0){
    print(paste("Able to recover environmental data for all", hab, "sites"))
  }else{
    print(paste("Unable to recover environmental data for", length(which(is.na(tmax))), hab, "sites"))
  }
}



#checking
which(acc$temp_range=='NA')

#write the file
#write.csv(x=acc, file = 'acc_data.csv')




