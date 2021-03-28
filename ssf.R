
# set working directory and load required packages
setwd("C:/Users/Tess Bronnvik/Desktop/Br-nnvik_honey_buzzard_ssf")
ssf_packs <- c("lubridate", "amt", "sf", "move", "dplyr", "purrr", "ggpubr", "ggplot2", "RNCEP", "mapview")
lapply(ssf_packs, require, character.only = TRUE)

# import required functions
NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}
source("C:/Users/Tess Bronnvik/Desktop/wind_support_Kami.R")

## set criteria for tracks
stepNumber <- 9 # random steps
stepLength <- 120 # minutes between bursts
toleranceLength <- 15 # tolerance in minutes
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")


### Prepare data for Movebank
#############################################################################################
full_data <- read.csv("full_buzz_set.csv",stringsAsFactors = F, header = T)
full_data$timestamp <- as.POSIXct(strptime(full_data$dt, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")


# get only autumn migrations
all_autumns <- full_data[grep("autumn", full_data$phase),]

# find and remove duplicate observations
ind2 <- all_autumns %>% select(long, lat, name, timestamp) %>%
  duplicated
sum(ind2) # 59 duplicates
all_autumns <- all_autumns[ind2!=TRUE,]

## order all the data by timestamp
all_autumns <- all_autumns %>% #group_by(name) %>% 
  arrange(timestamp)#, .by_group = T) %>%
  #ungroup()

first_autumns <- all_autumns[all_autumns$phase == "first.autumn",]
autumns <- all_autumns[all_autumns$phase == "autumn",]


## all the juveniles
Autumnsl <- first_autumns %>% filter(name == "Aida" | name == "Anni" | name == "Edit" | name == "Ella" | 
                                name == "Emma" | name == "Gilda" | name == "Hans" | name == "Heidi" | 
                                #name == "Jaana"| 
                                  name == "Julia" | name == "Kirsi"| name == "Lars"| name == "Lisa"|
                                name == "Matti"| name == "Miikka" | name == "Mohammed" | name == "Per"| 
                                #name == "Piff"| name == "Puff" | name == "Roosa"| 
                                  #name == "Rudolf" |
                                  #name == "Senta"| 
                                name == "Sven" | #name == "Tor"|
                                  name == "Ulla" | name == "Valentin"| name == "Venus" | name == "Viljo")

Autumnsl <- first_autumns %>% filter(name == "Jaana"| name == "Senta")

# do simple plots of lat/long to check for outliers
ggplot(Autumns, aes(x=long, y=lat)) + geom_point()+ facet_wrap(~name, scales="free")
ggplot(Autumns, aes(x=long, y=lat, color=as.factor(name))) + geom_point()

# create a move object for plotting
mv <- move(x = Autumns$long, y = Autumns$lat, time = Autumns$timestamp, data = Autumns, animal = Autumns$name, proj = wgs)
maps::map("world", xlim = c(-20,58), ylim = c(-25,70))
points(mv, cex = 0.3, pch = 16) #adds points to the previous plot

# create a spatial object for plotting
data_sp <- Autumns
coordinates(data_sp) <- ~ long + lat
proj4string(data_sp) <- wgs
mapView(data_sp, zcol = "cy")

# make the data of class track
#trk <- mk_track(Autumns, .x=long, .y=lat, .t=timestamp, id = name, crs = CRS("+init=epsg:4326"))

#transform to geographic coordinates
#trk <- transform_coords(trk, CRS("+init=epsg:3857")) # pseudo-Mercator in meters (https://epsg.io/3857)


## prepare data for the SSF 
set.seed(8675309)
autumn_track <- lapply(split(Autumnsl, Autumnsl$name), function(x){
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs) %>%
    #transform_coords(CRS("+init=epsg:3857")) %>% 
    track_resample(rate = minutes(stepLength), tolerance = minutes(toleranceLength))
  
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  
  
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
  #rnd_stps<-rnd_stps %>% mutate(id = x$name)
  
}) %>% reduce(rbind)
unique(autumn_track$id)
autumn_track <- as_tibble(autumn_track)
head(autumn_track)

all(complete.cases(autumn_track$case_))

# plot some comparisons of random and matched points
Anni <- autumn_track[autumn_track$id == "Anni",]
Edit <- autumn_track[autumn_track$id == "Edit",]
Emma <- autumn_track[autumn_track$id == "Emma",]
Gilda <- autumn_track[autumn_track$id == "Gilda",]
Anni_plot <- ggplot(Anni, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id, scales="free")
Edit_plot <- ggplot(Edit, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id, scales="free")
Emma_plot <- ggplot(Emma, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id, scales="free")
Gilda_plot <- ggplot(Gilda, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id, scales="free")
ggarrange(Anni_plot, Edit_plot, Emma_plot,  Gilda_plot, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

# plot the step lengths for observed and random steps
true_steps <- autumn_track[autumn_track$case_ == TRUE,]
false_steps <- autumn_track[autumn_track$case_ == FALSE,]
summary(false_steps$sl_)

obs <- true_steps %>% select(id, sl_) %>% 
  ggplot(aes(sl_, fill = factor(id))) + ggtitle("Observed") + geom_density(alpha = 0.4)
rand <- false_steps %>% select(id, sl_) %>% 
  ggplot(aes(sl_, fill = factor(id))) + ggtitle("Alternative") + geom_density(alpha = 0.4)

ggarrange(obs, rand, ncol=2, nrow=1, legend = "none")

# create a new gamma distribution
#library(tidyverse)
# Calculate summary statistics
#stats <- true_steps %>% summarise(Mean=mean(true_steps$sl_), Variance=var(true_steps$sl_))
#stats
# Derive parameter estimates
#scale <- stats$Variance / stats$Mean
#shape <-  stats$Mean / scale
# Derive fitted PDF and compare with empirical PDF
#fitted <- tibble(x=seq(from = 0,to = 53.26001, by = 0.25), y=dgamma(x, shape=shape, scale=scale))
#true_steps %>% ggplot + geom_histogram(aes(x=sl_, y=..density..), bins=120) + geom_line(data=fitted, aes(x=x, y=y), colour="blue")
#summary(fitted$x)


autumn_track %>%
  mutate(used = as.numeric(case_)) %>%
  filter(used == 0) %>%
  group_by(step_id_) %>%
  summarise(n=n()) %>%
  filter(n != 9)


#ggplot(autumn_track, aes(x2_, y2_, color=case_))+geom_point()+facet_wrap(~id, scales="free")


## export for Movebank annotation

#relabel the columns to clarify utms
#autumn_track$utm.easting<-autumn_track$x2_
#autumn_track$utm.northing<-autumn_track$y2_

autumn_track2 <- SpatialPointsDataFrame(autumn_track[,c("x2_","y2_")], autumn_track, proj4string = CRS("+init=epsg:3857"))  
ssf.df <- autumn_track#data.frame(spTransform(autumn_track2, CRS("+proj=longlat +datum=WGS84"))) 
names(ssf.df)[c(8,2,4)] <-c("individual.local.identifier", "location-long", "location-lat")
ssf.df$timestamp<-ssf.df$t1_
ssf.df$timestamp <- paste0(ssf.df$timestamp,".000" )
#write.csv(ssf.df, "12.juv.1A.60.15.csv", row.names = F)

## make a track for the subsequent years

# find and remove duplicate observations
ind2 <- autumns %>% select(long, lat, name, timestamp) %>%
  duplicated
sum(ind2) # 59 duplicates
autumns <- autumns[ind2!=TRUE,]

## order all the data by timestamp
autumns <- autumns %>% 
  arrange(timestamp)

# remove unwanted individuals 
autumn2 <- autumns %>% filter( 
    name == "Jaana"|
    name == "Lars" |
    name == "Mohammed" |
    name == "Senta"| 
    name == "Valentin")

## set criteria for tracks
stepNumber <- 29 # random steps
stepLength <- 120 # minutes between bursts
toleranceLength <- 15 # tolerance in minutes
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")


set.seed(8675309)
autumn2_track <- lapply(split(autumn2, autumn2$name), function(x){
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs) %>%
    #transform_coords(CRS("+init=epsg:3857")) %>% 
    track_resample(rate = minutes(stepLength), tolerance = minutes(toleranceLength))
  
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  
  
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
  #rnd_stps<-rnd_stps %>% mutate(id = x$name)
  
}) %>% reduce(rbind)
unique(autumn2_track$id)
autumn2_track <- as_tibble(autumn2_track)
head(autumn2_track)

#############################################################################################

## The for loop approach
#############################################################################################
# create the tracks, re-sample them, filter bursts, convert to steps, and create alternatives
ssfdat<-NULL
temptrk<-with(trk, track(x=x_, y=y_, t=t_, id=id))
uid<-unique(trk$id) # individual identifiers
luid<-length(uid) # number of unique individuals
for(i in 1:luid){
  # Subset individuals & regularize track
  temp<-temptrk%>% filter(id==uid[i]) %>% 
    track_resample(rate=minutes(stepLength), tolerance=minutes(toleranceLength))
  
  # Get rid of any bursts without at least 2 points
  temp<-filter_min_n_burst(temp, 3)
  
  # burst steps
  stepstemp<-steps_by_burst(temp)
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- stepstemp %>%  random_steps(n = stepNumber)
  
  # append id
  rnd_stps<-rnd_stps%>%mutate(id=uid[i])
  
  # append new data to data from other individuals
  ssfdat<-rbind(rnd_stps, ssfdat)
}
ssfdat<-as_tibble(ssfdat)
ssfdat

#############################################################################################


### Prepare data for models
#############################################################################################
ssfdata <- read.csv("juvenile_autumns-895583402217231081.csv", stringsAsFactors = F)
ssfdata$timestamp <- as.POSIXct(strptime(ssfdata$timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
ssfdata$year <- format(as.POSIXct(ssfdata$timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y")
ssfdata$case_ <- as.numeric(ssfdata$case_)


# get the years for each individual
ssfdata$a <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Aida", "1", NA)
ssfdata$b <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Anni", "1", NA)
ssfdata$c <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Edit", "1", NA)
ssfdata$d <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Ella", "1", NA)
ssfdata$e <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Emma", "1", NA)
ssfdata$f <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Gilda", "1", NA)
ssfdata$g <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Hans", "1", NA)
ssfdata$h <- ifelse(ssfdata$year == "2011" & ssfdata$id == "Jaana", "1", NA)
ssfdata$i <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Jaana", "2", NA)
ssfdata$j <- ifelse(ssfdata$year == "2014" & ssfdata$id == "Jaana", "3", NA)
ssfdata$k <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Julia", "1", NA)
ssfdata$l <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Kirsi", "1", NA)
ssfdata$m <- ifelse(ssfdata$year == "2011" & ssfdata$id == "Lars", "1", NA)
ssfdata$n <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Lars", "2", NA)
ssfdata$o <- ifelse(ssfdata$year == "2014" & ssfdata$id == "Lars", "3", NA)   
ssfdata$p <- ifelse(ssfdata$year == "2015" & ssfdata$id == "Lars", "4", NA)
ssfdata$q <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Lisa", "1", NA)
ssfdata$r <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Matti", "1", NA)
ssfdata$s <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Miikka", "1", NA)
ssfdata$t <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Mohammed", "1", NA)
ssfdata$u <- ifelse(ssfdata$year == "2018" & ssfdata$id == "Mohammed", "2", NA)
ssfdata$v <- ifelse(ssfdata$year == "2019" & ssfdata$id == "Mohammed", "3", NA)   
ssfdata$w <- ifelse(ssfdata$year == "2020" & ssfdata$id == "Mohammed", "4", NA)
ssfdata$x <- ifelse(ssfdata$year == "2016" & ssfdata$id == "Per", "1", NA)

ssfdata$order <- NA
ssfdata$order[!is.na(ssfdata$a)] = ssfdata$a[!is.na(ssfdata$a)]
ssfdata$order[!is.na(ssfdata$b)] = ssfdata$b[!is.na(ssfdata$b)]
ssfdata$order[!is.na(ssfdata$c)] = ssfdata$c[!is.na(ssfdata$c)]    
ssfdata$order[!is.na(ssfdata$d)] = ssfdata$d[!is.na(ssfdata$d)]    
ssfdata$order[!is.na(ssfdata$e)] = ssfdata$e[!is.na(ssfdata$e)]    
ssfdata$order[!is.na(ssfdata$f)] = ssfdata$f[!is.na(ssfdata$f)]    
ssfdata$order[!is.na(ssfdata$g)] = ssfdata$g[!is.na(ssfdata$g)]    
ssfdata$order[!is.na(ssfdata$h)] = ssfdata$h[!is.na(ssfdata$h)]    
ssfdata$order[!is.na(ssfdata$i)] = ssfdata$i[!is.na(ssfdata$i)]    
ssfdata$order[!is.na(ssfdata$j)] = ssfdata$j[!is.na(ssfdata$j)]    
ssfdata$order[!is.na(ssfdata$k)] = ssfdata$k[!is.na(ssfdata$k)]    
ssfdata$order[!is.na(ssfdata$l)] = ssfdata$l[!is.na(ssfdata$l)]    
ssfdata$order[!is.na(ssfdata$m)] = ssfdata$m[!is.na(ssfdata$m)]    
ssfdata$order[!is.na(ssfdata$n)] = ssfdata$n[!is.na(ssfdata$n)]    
ssfdata$order[!is.na(ssfdata$o)] = ssfdata$o[!is.na(ssfdata$o)]    
ssfdata$order[!is.na(ssfdata$p)] = ssfdata$p[!is.na(ssfdata$p)]    
ssfdata$order[!is.na(ssfdata$q)] = ssfdata$q[!is.na(ssfdata$q)]
ssfdata$order[!is.na(ssfdata$r)] = ssfdata$r[!is.na(ssfdata$r)]
ssfdata$order[!is.na(ssfdata$s)] = ssfdata$s[!is.na(ssfdata$s)]
ssfdata$order[!is.na(ssfdata$t)] = ssfdata$t[!is.na(ssfdata$t)]
ssfdata$order[!is.na(ssfdata$u)] = ssfdata$u[!is.na(ssfdata$u)]
ssfdata$order[!is.na(ssfdata$v)] = ssfdata$v[!is.na(ssfdata$v)]
ssfdata$order[!is.na(ssfdata$w)] = ssfdata$w[!is.na(ssfdata$w)]
ssfdata$order[!is.na(ssfdata$x)] = ssfdata$x[!is.na(ssfdata$x)]


ssfdata$a <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Rudolf", "1", NA)
ssfdata$b <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Senta", "1", NA)   
ssfdata$c <- ifelse(ssfdata$year == "2015" & ssfdata$id == "Senta", "2", NA)
ssfdata$d <- ifelse(ssfdata$year == "2016" & ssfdata$id == "Senta", "3", NA)
ssfdata$e <- ifelse(ssfdata$year == "2017" & ssfdata$id == "Senta", "4", NA)
ssfdata$f <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Sven", "1", NA)
ssfdata$g <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Tor", "1", NA)
ssfdata$h <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Ulla", "1", NA)
ssfdata$i <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Valentin", "1", NA)   
ssfdata$j <- ifelse(ssfdata$year == "2015" & ssfdata$id == "Valentin", "2", NA)
ssfdata$k <- ifelse(ssfdata$year == "2012" & ssfdata$id == "Venus", "1", NA)
ssfdata$l <- ifelse(ssfdata$year == "2014" & ssfdata$id == "Viljo", "1", NA)

ssfdata$order[!is.na(ssfdata$a)] = ssfdata$a[!is.na(ssfdata$a)]
ssfdata$order[!is.na(ssfdata$b)] = ssfdata$b[!is.na(ssfdata$b)]
ssfdata$order[!is.na(ssfdata$c)] = ssfdata$c[!is.na(ssfdata$c)]    
ssfdata$order[!is.na(ssfdata$d)] = ssfdata$d[!is.na(ssfdata$d)]    
ssfdata$order[!is.na(ssfdata$e)] = ssfdata$e[!is.na(ssfdata$e)]    
ssfdata$order[!is.na(ssfdata$f)] = ssfdata$f[!is.na(ssfdata$f)]    
ssfdata$order[!is.na(ssfdata$g)] = ssfdata$g[!is.na(ssfdata$g)]    
ssfdata$order[!is.na(ssfdata$h)] = ssfdata$h[!is.na(ssfdata$h)]    
ssfdata$order[!is.na(ssfdata$i)] = ssfdata$i[!is.na(ssfdata$i)]    
ssfdata$order[!is.na(ssfdata$j)] = ssfdata$j[!is.na(ssfdata$j)]    
ssfdata$order[!is.na(ssfdata$k)] = ssfdata$k[!is.na(ssfdata$k)]    
ssfdata$order[!is.na(ssfdata$l)] = ssfdata$l[!is.na(ssfdata$l)]   

colnames(ssfdata)
#ssfdata[c(25:48)] <- NULL
ssfdata <- ssfdata[ssfdata$id != "Jaana",]

## derive the cross and tail wind components
colnames(ssfdata)[c(1,2,3,4,16,17)] <- c("lon1","lon2","lat1","lat2","u_wind","v_wind")


ssfdata$heading <- NCEP.loxodrome.na(lat1 = ssfdata$lat1, lat2 = ssfdata$lat2, lon1 = ssfdata$lon1, lon2 = ssfdata$lon2)
ssfdata$tail <- wind_support(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)
ssfdata$cross <- cross_wind(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)

summary(ssfdata) # 254 NAs in heading

### normalize the data
ssfdata <- ssfdata %>% 
  #group_by(id) %>% 
  mutate(scaled_cross = abs(as.numeric(scale(cross, center = T, scale = T))),
         scaled_tail = as.numeric(scale(tail, center = T, scale = T))) %>%
  ungroup()
#############################################################################################
data_sf <- ssfdata %>% st_as_sf(coords = c("lon2", "lat2"), crs = wgs)


### Build the models
#############################################################################################

# nest the data
ind_by_year <- ssfdata %>% 
  group_by(id,year) %>% 
  nest()

# build the model function
ind_model <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_), data = df)
}

# model
ind_by_year <- ind_by_year %>% 
  mutate(model = purrr::map(data, ind_model))

ind_by_year2 <- ind_by_year %>% 
  mutate(coefs = purrr::map(model, coef)) %>% 
  unnest(model) %>%
  unnest(data)

#ind_by_year2 <- ind_by_year2 %>% separate(coefs, c("scaled_tail", "scaled_cross"))

coefs <- lapply(split(ind_by_year, ind_by_year$id, ind_by_year$model), function(x){
         coef(summary(x$model))["scaled_tail",]
}) %>% reduce(rbind)

## build individual models
modelAida1 <- ssfdata %>% 
  dplyr::filter(id == "Aida" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelAnni1 <- ssfdata %>% 
  dplyr::filter(id == "Anni" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelEdit1 <- ssfdata %>% 
  dplyr::filter(id == "Edit" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelElla1 <- ssfdata %>% 
  dplyr::filter(id == "Ella" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelEmma1 <- ssfdata %>% 
  dplyr::filter(id == "Emma" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelGilda1 <- ssfdata %>% 
  dplyr::filter(id == "Gilda" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelHans1 <- ssfdata %>% 
  dplyr::filter(id == "Hans" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelJaana1 <- ssfdata %>% 
#  dplyr::filter(id == "Jaana" & order == "1") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelJaana2 <- ssfdata %>% 
  dplyr::filter(id == "Jaana" & order == "2") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelJaana3 <- ssfdata %>% 
#  dplyr::filter(id == "Jaana" & order == "3") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelJulia1 <- ssfdata %>% 
  dplyr::filter(id == "Julia" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelKirsi1 <- ssfdata %>% 
  dplyr::filter(id == "Kirsi" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelLars1 <- ssfdata %>% 
#  dplyr::filter(id == "Lars" & order == "1") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelLars2 <- ssfdata %>% 
  dplyr::filter(id == "Lars" & order == "2") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelLars3 <- ssfdata %>% 
#  dplyr::filter(id == "Lars" & order == "3") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelLars4 <- ssfdata %>% 
  dplyr::filter(id == "Lars" & order == "4") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelLisa1 <- ssfdata %>% 
  dplyr::filter(id == "Lisa" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelMatti1 <- ssfdata %>% 
  dplyr::filter(id == "Matti" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelMiikka1 <- ssfdata %>% 
  dplyr::filter(id == "Miikka" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelMohammed1 <- ssfdata %>% 
  dplyr::filter(id == "Mohammed" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelMohammed2 <- ssfdata %>% 
  dplyr::filter(id == "Mohammed" & order == "2") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelMohammed3 <- ssfdata %>% 
#  dplyr::filter(id == "Mohammed" & order == "3") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelMohammed4 <- ssfdata %>% 
#  dplyr::filter(id == "Mohammed" & order == "4") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelPer1 <- ssfdata %>% 
  dplyr::filter(id == "Per" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelRudolf1 <- ssfdata %>% 
  dplyr::filter(id == "Rudolf" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelSenta1 <- ssfdata %>% 
  dplyr::filter(id == "Senta" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelSenta2 <- ssfdata %>% 
  dplyr::filter(id == "Senta" & order == "2") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelSenta3 <- ssfdata %>% 
  dplyr::filter(id == "Senta" & order == "3") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

#modelSenta4 <- ssfdata %>% 
#  dplyr::filter(id == "Senta" & order == "4") %>%
#  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelSven1 <- ssfdata %>% 
  dplyr::filter(id == "Sven" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelTor <- ssfdata %>% 
  dplyr::filter(id == "Tor" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelUlla1 <- ssfdata %>% 
  dplyr::filter(id == "Ulla" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelValentin1 <- ssfdata %>% 
  dplyr::filter(id == "Valentin" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelValentin2 <- ssfdata %>% 
  dplyr::filter(id == "Valentin" & order == "2") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelVenus1 <- ssfdata %>% 
  dplyr::filter(id == "Venus" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))

modelViljo1 <- ssfdata %>% 
  dplyr::filter(id == "Viljo" & order == "1") %>%
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_))
#############################################################################################


### Extract model coefficients for plotting
#############################################################################################

model_ids <- c(modelAida1, modelAnni1,modelEdit1,modelElla1,modelEmma1,modelGilda1,
               modelHans1,#modelJaana1,
               #modelJaana2,modelJaana3,
               modelJulia1,modelKirsi1,# modelLars1,
               modelLars2,#modelLars3,
               modelLars4,modelLisa1,modelMatti1,
               modelMiikka1,modelMohammed1,modelMohammed2,#modelMohammed3, modelMohammed4, 
               modelPer1,modelRudolf1,modelSenta1,modelSenta2,modelSenta3,#modelSenta4,
               modelSven1,modelTor,modelUlla1,modelValentin1,modelValentin2,modelVenus1,
               modelViljo1)

coefs <- lapply(model_ids, function(x){
  coefs <- coef(summary(x))["scaled_tail",]
}) %>% reduce(rbind)



AI1t <- coef(summary(modelAida1))["scaled_tail",]
AN1t <- coef(summary(modelAnni1))["scaled_tail",]
ED1t <- coef(summary(modelEdit1))["scaled_tail",]
EL1t <- coef(summary(modelElla1))["scaled_tail",]
EM1t <- coef(summary(modelEmma1))["scaled_tail",]
GI1t <- coef(summary(modelGilda1))["scaled_tail",]
HA1t <- coef(summary(modelHans1))["scaled_tail",]
JA1t <- coef(summary(modelJaana1))["scaled_tail",]
JA2t <- coef(summary(modelJaana2))["scaled_tail",]
JA3t <- coef(summary(modelJaana3))["scaled_tail",]
JL1t <- coef(summary(modelJulia1))["scaled_tail",]
KI1t <- coef(summary(modelKirsi1))["scaled_tail",]
LA1t <- coef(summary(modelLars1))["scaled_tail",]
LA2t <- coef(summary(modelLars2))["scaled_tail",]
LA3t <- coef(summary(modelLars3))["scaled_tail",]
LA4t <- coef(summary(modelLars4))["scaled_tail",]
LI1t <- coef(summary(modelLisa1))["scaled_tail",]
MA1t <- coef(summary(modelMatti1))["scaled_tail",]
MII1t <- coef(summary(modelMiikka1))["scaled_tail",]
MO1t <- coef(summary(modelMohammed1))["scaled_tail",]
MO2t <- coef(summary(modelMohammed2))["scaled_tail",]
MO3t <- coef(summary(modelMohammed3))["scaled_tail",]
MO4t <- coef(summary(modelMohammed4))["scaled_tail",]
PE1t <- coef(summary(modelPer1))["scaled_tail",]
RU1t <- coef(summary(modelRudolf1))["scaled_tail",]
SE1t <- coef(summary(modelSenta1))["scaled_tail",]
SE2t <- coef(summary(modelSenta2))["scaled_tail",]
SE3t <- coef(summary(modelSenta3))["scaled_tail",]
SE4t <- coef(summary(modelSenta4))["scaled_tail",]
SV1t <- coef(summary(modelSven1))["scaled_tail",]
TO1t <- coef(summary(modelTor))["scaled_tail",]
UL1t <- coef(summary(modelUlla1))["scaled_tail",]
VA1t <- coef(summary(modelValentin1))["scaled_tail",]
VA2t <- coef(summary(modelValentin2))["scaled_tail",]
VE1t <- coef(summary(modelVenus1))["scaled_tail",]
VI1t <- coef(summary(modelViljo1))["scaled_tail",]

tails_stats <- rbind(AI1t,AN1t,ED1t,EL1t,EM1t,GI1t,HA1t,#JA1t,
                     #JA2t,#JA3t,
                     JL1t,KI1t,#LA1t,
                     LA2t,#LA3t,
                     LA4t,LI1t,MA1t,MII1t,MO1t,MO2t,#MO3t,MO4t,
                     PE1t,RU1t,SE1t,SE2t,SE3t,#SE4t,
                     SV1t,TO1t,UL1t,VA1t,VA2t,VE1t,VI1t)
colnames(tails_stats)[] <- c("tail_coef", "tail_exp(coef)","tail_se(coef)", "tail_z", "tail_Pr(>|z|)")


AI1c <- coef(summary(modelAida1))["scaled_cross",]
AN1c <- coef(summary(modelAnni1))["scaled_cross",]
ED1c <- coef(summary(modelEdit1))["scaled_cross",]
EL1c <- coef(summary(modelElla1))["scaled_cross",]
EM1c <- coef(summary(modelEmma1))["scaled_cross",]
GI1c <- coef(summary(modelGilda1))["scaled_cross",]
HA1c <- coef(summary(modelHans1))["scaled_cross",]
JA1c <- coef(summary(modelJaana1))["scaled_cross",]
JA2c <- coef(summary(modelJaana2))["scaled_cross",]
JA3c <- coef(summary(modelJaana3))["scaled_cross",]
JL1c <- coef(summary(modelJulia1))["scaled_cross",]
KI1c <- coef(summary(modelKirsi1))["scaled_cross",]
LA1c <- coef(summary(modelLars1))["scaled_cross",]
LA2c <- coef(summary(modelLars2))["scaled_cross",]
LA3c <- coef(summary(modelLars3))["scaled_cross",]
LA4c <- coef(summary(modelLars4))["scaled_cross",]
LI1c <- coef(summary(modelLisa1))["scaled_cross",]
MA1c <- coef(summary(modelMatti1))["scaled_cross",]
MII1c <- coef(summary(modelMiikka1))["scaled_cross",]
MO1c <- coef(summary(modelMohammed1))["scaled_cross",]
MO2c <- coef(summary(modelMohammed2))["scaled_cross",]
MO3c <- coef(summary(modelMohammed3))["scaled_cross",]
MO4c <- coef(summary(modelMohammed4))["scaled_cross",]
PE1c <- coef(summary(modelPer1))["scaled_cross",]
RU1c <- coef(summary(modelRudolf1))["scaled_cross",]
SE1c <- coef(summary(modelSenta1))["scaled_cross",]
SE2c <- coef(summary(modelSenta2))["scaled_cross",]
SE3c <- coef(summary(modelSenta3))["scaled_cross",]
SE4c <- coef(summary(modelSenta4))["scaled_cross",]
SV1c <- coef(summary(modelSven1))["scaled_cross",]
TO1c <- coef(summary(modelTor))["scaled_cross",]
UL1c <- coef(summary(modelUlla1))["scaled_cross",]
VA1c <- coef(summary(modelValentin1))["scaled_cross",]
VA2c <- coef(summary(modelValentin2))["scaled_cross",]
VE1c <- coef(summary(modelVenus1))["scaled_cross",]
VI1c <- coef(summary(modelViljo1))["scaled_cross",]

crosses_stats <- rbind(AI1c,AN1c,ED1c,EL1c,EM1c,GI1c,HA1c,#JA1c,
                     #JA2c,#JA3c,
                     JL1c,KI1c,#LA1c,
                     LA2c,#LA3c,
                     LA4c,LI1c,MA1c,MII1c,MO1c,MO2c,#MO3c,MO4c,
                     PE1c,RU1c,SE1c,SE2c, SE3c,#SE4c,
                     SV1c,TO1c,UL1c,VA1c,VA2c,VE1c,VI1c)
colnames(crosses_stats)[] <- c("cross_coef", "cross_exp(coef)", "cross_se(coef)", "cross_z", "cross_Pr(>|z|)")

coefs_df <- NULL
coefs_df$id <- c("Aida", "Anni", "Edit", "Ella", "Emma", "Gilda", "Hans",# "Jaana", 
                 "Julia", "Kirsi", "Lars", "Lars", "Lisa", "Matti", "Miikka", "Mohammed", "Mohammed",
                 "Per", "Rudolf", "Senta", "Senta", "Senta","Sven", "Tor", "Ulla","Valentin",
                 "Valentin", "Venus","Viljo")
coefs_df$autumn <- c(1,1,1,1,1,1,1,1,1,2,4,1,1,1,1,2,1,1,1,2,3,1,1,1,1,2,1,1)
coefs_df <- as.data.frame(coefs_df)

wind_stats <- cbind(coefs_df, tails_stats, crosses_stats)
#wind_stats <- wind_stats[wind_stats$autumn == 1,] # select which years to plot
#write.csv(wind_stats, "wind_stats.csv", row.names = F)
#wind_stats <- read.csv("wind_stats.csv", stringsAsFactors = F)
#############################################################################################


### Plot the wind coefficients
#############################################################################################

tailwind_plot <-  ggplot(data = wind_stats, aes(x = id, y = tail_coef, colour = id, group = id)) + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  geom_point(size = 2.5) + ylim(-3,3) +
  #geom_smooth(method=lm, se = F) # Don't add shaded confidence region
  geom_hline(yintercept = 0)


crosswind_plot <-  ggplot(data = wind_stats, aes(x = id, y = cross_coef, colour = id, group = id)) + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.position  = "none") +
        geom_point(size = 2.5) +
        ylim(-3,3)

ggarrange(tailwind_plot, crosswind_plot, ncol=2, nrow=1, legend = "none")# common.legend = TRUE, legend="right")
#############################################################################################

