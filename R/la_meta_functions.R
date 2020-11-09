require(readxl)
require(purrr)
require(rlist)
require(sdmpredictors)
require(raster)
require(ncdf4)
require(ggplot2)
require(maps)
require(viridis)  
require(ggthemes)
require(paletteer)
require(vegan)
require(metafor)
require(here)
require(tidyverse)
require(ggExtra)
require(ggpubr)
require(broom)

#This function takes two sets of lat-long coordinates and determines the distance between them
coord_to_dist = function(lat1, long1, lat2, long2, units = "km") {
  loadNamespace("purrr")
  loadNamespace("geosphere")
  loadNamespace("rlist")
  longlat1 = purrr::map2(long1, lat1, function(x,y) c(x,y))
  longlat2 = purrr::map2(long2, lat2, function(x,y) c(x,y))
  distance_list = purrr::map2(longlat1, longlat2, function(x,y) geosphere::distHaversine(x, y))
  distance_m = rlist::list.extract(distance_list, position = 1)
  if (units == "km") {
    distance = distance_m / 1000.0;
  }
  else if (units == "miles") {
    distance = distance_m / 1609.344
  }
  else {
    distance = distance_m
    # This will return in meters 
  }
  distance
}

#This function can be used to calculate standard error via aggregrate when there's missing data
st.err = function(x, na.rm=FALSE) {
  if(na.rm==TRUE) x = na.omit(x)
  sd(x)/sqrt(length(x))
}

#This function is used to calculate varcov matrices for use with standardized mean difference
calc.v <- function(x) {
  v <- matrix(1/x$n2[1] + outer(x$yi, x$yi, "*")/(2*x$Ni[1]), nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$vi
  v 
}