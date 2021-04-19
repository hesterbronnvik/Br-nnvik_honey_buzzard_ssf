
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
stepNumber <- 99 # random steps
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
                                  name == "Rudolf" |
                                  #name == "Senta"| 
                                name == "Sven" | #name == "Tor"|
                                  name == "Ulla" | name == "Valentin"| name == "Venus" | name == "Viljo")
## all the adults
Autumnsl <- autumns %>% filter(name == "Aarne" | name == "Annika" |name == "Jouko" |name == "Kari" |name == "Mikko" |name == "Paivi")

## the multi-year juveniles
Autumnsl <- autumns %>% filter(name == "Lars" | name == "Mohammed" |name == "Valentin" |name == "Senta" |name == "Jaana")

Autumnsl$Migration <- "First"
Autumnsl <- Autumnsl %>% select(lat,long,Migration, name)

all_autumns$id_year <- paste(all_autumns$name,all_autumns$yr,sep="_")
secaut <- all_autumns[which(all_autumns$id_year == "Jaana_2013" | 
                              all_autumns$id_year == "Lars_2013" | 
                                    all_autumns$id_year == "Mohammed_2018" | 
                                    all_autumns$id_year == "Senta_2015" | 
                                    all_autumns$id_year == "Valentin_2015"),]
secaut$Migration <- "Second"
secaut <- secaut %>% select(lat,long,Migration, name)
thiraut <- all_autumns[which(all_autumns$id_year == "Jaana_2014" | 
                                    all_autumns$id_year == "Lars_2014" | 
                                    all_autumns$id_year == "Senta_2016"),]
thiraut$Migration <- "Third"
thiraut <- thiraut %>% select(lat,long,Migration, name)

fouraut <- all_autumns[which(all_autumns$id_year == "Lars_2015" |
                                    all_autumns$id_year == "Senta_2017"),]
fouraut$Migration <- "Fourth"
fouraut <- fouraut %>% select(lat,long,Migration, name)

map_autumns <- rbind(Autumnsl,secaut,thiraut,fouraut)

# do simple plots of lat/long to check for outliers
ggplot(Autumnsl, aes(x=long, y=lat)) + geom_point()+ facet_wrap(~name, scales="free")
ggplot(Autumnsl, aes(x=long, y=lat, color=as.factor(name))) + geom_point()

# create a move object for plotting
map_autumns <- Autumnsl %>% group_by(name) %>% 
  arrange(timestamp, .by_group = T) %>%
  ungroup()
mv <- move(x = map_autumns$long, y = map_autumns$lat, time = map_autumns$timestamp, data = map_autumns, animal = map_autumns$name, proj = wgs)
maps::map("world", xlim = c(-20,58), ylim = c(-25,70))
points(mv, cex = 0.3, pch = 16) #adds points to the previous plot

# create a spatial object for plotting
data_sp <- map_autumns#[map_autumns$yr == "2014",]
coordinates(data_sp) <- ~ long + lat
proj4string(data_sp) <- wgs
mapView(data_sp, zcol = "Migration", burst = T, cex = 3, color = rainbow)
#mapshot(all_track_map,file = paste0(getwd(), "/map.png"),remove_controls = c("zoomControl", "layersControl"))

library(ggmap)
require(mapproj)
A_mapDF <- as.data.frame(map_autumns)
loc_point <- c(left = -18, bottom = -15, right = 49, top = 66)#c(9.186824662880852,47.69034606201426)
myMap <- get_map(location = loc_point, source="stamen", zoom=5, crop = F, maptype = "watercolor")
ggmap(myMap)+
  geom_path(data=A_mapDF, aes(x=long, y=lat, color = name)) + 
  #scale_colour_manual(values=rainbow(37))+
  geom_point(alpha = .5)+
  scale_size(range=c(3,20)) +
  theme(legend.position = "none") +
  labs(x="Longitude", y = "Latitude")

# make the data of class track
#trk <- mk_track(Autumns, .x=long, .y=lat, .t=timestamp, id = name, crs = CRS("+init=epsg:4326"))

#transform to geographic coordinates
#trk <- transform_coords(trk, CRS("+init=epsg:3857")) # pseudo-Mercator in meters (https://epsg.io/3857)


## prepare data for the SSF 
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
names(ssf.df)[c(2,4)] <-c("location-long", "location-lat")
ssf.df$timestamp<-ssf.df$t1_
ssf.df$timestamp <- paste0(ssf.df$timestamp,".000" )
write.csv(ssf.df, "all.adults.99s.120.15.csv", row.names = F)

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
juveniles <- read.csv("all.juv.1A.99s.240&120.15-5221204997286432070.csv", stringsAsFactors = F)
juveniles$timestamp <- as.POSIXct(strptime(juveniles$timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
juveniles$year <- format(as.POSIXct(juveniles$timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y")
juveniles$id_year <- paste(juveniles$id,juveniles$year,sep="_")
juveniles$Migration <- "First"

juveniles_experienced <- read.csv("all.juv.234A.99s.120.15-5922315139798958994.csv", stringsAsFactors = F)
juveniles_experienced$timestamp <- as.POSIXct(strptime(juveniles_experienced$timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
juveniles_experienced$year <- format(as.POSIXct(juveniles_experienced$timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y")
juveniles_experienced$id_year <- paste(juveniles_experienced$id,juveniles_experienced$year,sep="_")

juveniles_experienced$Migration <- NA
juveniles_experienced$Migration[which(juveniles_experienced$id_year == "Jaana_2013" | 
                                        juveniles_experienced$id_year == "Lars_2013" | 
                                        juveniles_experienced$id_year == "Mohammed_2018" |
                                        juveniles_experienced$id_year == "Senta_2015" |
                                        juveniles_experienced$id_year == "Valentin_2015")] <- "Second"

juveniles_experienced$Migration[which(juveniles_experienced$id_year == "Jaana_2014" | 
                                        juveniles_experienced$id_year == "Lars_2014" | 
                                        juveniles_experienced$id_year == "Senta_2016")] <- "Third"

juveniles_experienced$Migration[which(juveniles_experienced$id_year == "Lars_2015" |
                                        juveniles_experienced$id_year == "Senta_2017")] <- "Fourth"

adults <- read.csv("all.adults.99s.120.15-579802980818675201.csv", stringsAsFactors = F)
adults$timestamp <- as.POSIXct(strptime(adults$timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
adults$year <- format(as.POSIXct(adults$timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y")
adults$id_year <- paste(adults$id,adults$year,sep="_")
adults$Migration <- ">4"
adults <- adults %>% filter(id_year == "Annika_2015" |
                            id_year == "Jouko_2015" |
                            id_year == "Mikko_2015" |
                            id_year == "Paivi_2017")

ssfdata <- rbind(juveniles,juveniles_experienced,adults)
ssfdata$case_ <- as.numeric(ssfdata$case_)


# remove Lisa and Sven, the outliers in coefs
#ssfdata <- ssfdata %>% filter(id != "Sven" & id != "Lisa" & id != "Kirsi")

## how many steps per day left?
juveniles$ymd <- format(as.POSIXct(juveniles$timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y-%m-%d")
perday <- ssfdata %>%  
  group_by(ymd, id) %>% 
  mutate(stepno = max(step_id_)) %>% 
  ungroup()


## derive the cross and tail wind components
colnames(ssfdata)[c(1,2,3,4,16,17,20,21)] <- c("lon1","lon2","lat1","lat2","u_wind","v_wind","orographic","thermal")


ssfdata$heading <- NCEP.loxodrome.na(lat1 = ssfdata$lat1, lat2 = ssfdata$lat2, lon1 = ssfdata$lon1, lon2 = ssfdata$lon2)
ssfdata$tail <- wind_support(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)
ssfdata$cross <- cross_wind(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)


# find and plot NA values in thermal uplift
thermalNA <- is.na(ssfdata$thermal)
sum(thermalNA) # 162 NA values
missing_thermals <- ssfdata[thermalNA == T,] # select the observations missing values
coordinates(missing_thermals) <- ~ lon1 + lat1 # prepare for a map
proj4string(missing_thermals) <- wgs
mapView(missing_thermals, zcol = "id")

# which individual has the most missing values
sum(with(missing_thermals,id == "Edit"))# Edit has 94 missing values in step 16
sum(with(missing_thermals,id == "Annika")) # 57 (48 in 177 and 9 in 171)
sum(with(missing_thermals,id == "Miikka")) # 4
sum(with(missing_thermals,id == "Mohammed")) # 2
sum(with(missing_thermals,id == "Emma")) # 1
sum(with(missing_thermals,id == "Matti")) # 1
sum(with(missing_thermals,id == "Roosa")) # 1
sum(with(missing_thermals,id == "Tor")) # 1

## drop rows with NA in thermal column and then correct for new step numbers

# remove the missing values
thermalNA <- is.na(ssfdata$thermal)
ssfdata <- ssfdata[thermalNA != T,]

# select data with missing values
Edit_ssfdata <- ssfdata[ssfdata$id == "Edit",] # burst 1 step 16 (True) & burst 3 step 380
Annika_ssfdata <- ssfdata[ssfdata$id == "Annika",] # burst 1 417 step 171 (True) & burst 419 step 177

#remove those data from the data set
ssfdata <- ssfdata %>% filter(id != "Edit" &
                              id != "Annika")

# edit each ID with missing data to have the right number of steps remaining
Edit_ssfdata <- Edit_ssfdata[which(Edit_ssfdata$step_id_ != 16),] # remove step with 94 NAs
Annika_ssfdata <- Annika_ssfdata[which(Annika_ssfdata$step_id_ != 171),] # remove step with 9 NAs

# add the IDs back to the data
ssfdata <- rbind(ssfdata,Edit_ssfdata,Annika_ssfdata)
  
# select 50 steps from each step_id_. Choosing the first 50 allows protection of True steps.
ssfdata <- ssfdata %>% 
  group_by(id_year, burst_, step_id_) %>% 
  slice(1:50) %>% 
  ungroup()


### normalize the data
ssfdata <- ssfdata %>% 
  #group_by(id) %>% 
  mutate(scaled_cross = abs(as.numeric(scale(cross, center = T, scale = T))),
         scaled_tail = as.numeric(scale(tail, center = T, scale = T)),
         scaled_thermal = scale(thermal, center = T, scale = T),
         scaled_orographic = scale(orographic, center = T, scale = T)) %>%
  ungroup()
#############################################################################################
data_sf <- ssfdata %>% st_as_sf(coords = c("lon2", "lat2"), crs = wgs)


### Build the models
#############################################################################################

# build the model function
ssf_model4 <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_orographic + strata(step_id_), data = df)
}

ssf_model3 <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_cross + scaled_orographic + strata(step_id_), data = df)
}

ssf_model2 <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_), data = df)
}


# nest the data, then apply above models to each individual for each year and extract AIC
SSF_results <- ssfdata %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_model4 = purrr::map(data, ssf_model4), # add each model as a column of lists 
         ssf_model3 = purrr::map(data, ssf_model3),
         ssf_model2 = purrr::map(data, ssf_model2),
         ssf_coefs4 = purrr::map(ssf_model4, coef), # extract the coefficients for each model
         ssf_coefs3 = purrr::map(ssf_model3, coef),
         ssf_coefs2 = purrr::map(ssf_model2, coef),
         AIC_4 = map_dbl(ssf_model4, ~AIC(.)),
         AIC_3 = map_dbl(ssf_model3, ~AIC(.)),
         AIC_2 = map_dbl(ssf_model2, ~AIC(.)))

# For model selection, build 16 models to cover all combinations of variables
ssf_modelCO <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_orographic + strata(step_id_), data = df)
}
ssf_modelCOt <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_orographic + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelCt <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelCTt <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_tail + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelTO <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_orographic + strata(step_id_), data = df)
}
ssf_modelTOt <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_orographic + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelTt <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelOt <- function(df) {
  fit_clogit(case_ ~ scaled_orographic + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelC <- function(df) {
  fit_clogit(case_ ~ scaled_cross + strata(step_id_), data = df)
}
ssf_modelT <- function(df) {
  fit_clogit(case_ ~ scaled_tail + strata(step_id_), data = df)
}
ssf_modelO <- function(df) {
  fit_clogit(case_ ~ scaled_orographic + strata(step_id_), data = df)
}
ssf_modelt <- function(df) {
  fit_clogit(case_ ~ scaled_thermal + strata(step_id_), data = df)
}

SSF_selections <- ssfdata %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_model4 = purrr::map(data, ssf_model4), # add each model as a column of lists 
         ssf_model3 = purrr::map(data, ssf_model3),
         ssf_model2 = purrr::map(data, ssf_model2),
         ssf_modelCO = purrr::map(data, ssf_modelCO),
         ssf_modelCOt = purrr::map(data, ssf_modelCOt),
         ssf_modelCt = purrr::map(data, ssf_modelCt),
         ssf_modelCTt = purrr::map(data, ssf_modelCTt),
         ssf_modelTO = purrr::map(data, ssf_modelTO),
         ssf_modelTOt = purrr::map(data, ssf_modelTOt),
         ssf_modelTt = purrr::map(data, ssf_modelTt),
         ssf_modelOt = purrr::map(data, ssf_modelOt),
         ssf_modelC = purrr::map(data, ssf_modelC),
         ssf_modelT = purrr::map(data, ssf_modelT),
         ssf_modelO = purrr::map(data, ssf_modelO),
         ssf_modelt = purrr::map(data, ssf_modelt),
         tail_cross_thermal_orographic = map_dbl(ssf_model4, ~AIC(.)), # extract AIC for each model
         tail_cross_orographic = map_dbl(ssf_model3, ~AIC(.)),
         tail_cross = map_dbl(ssf_model2, ~AIC(.)),
         cross_orographic = map_dbl(ssf_modelCO, ~AIC(.)),
         cross_thermal_orographic = map_dbl(ssf_modelCOt, ~AIC(.)),
         cross_thermal = map_dbl(ssf_modelCt, ~AIC(.)),
         tail_cross_thermal = map_dbl(ssf_modelCTt, ~AIC(.)),
         tail_orographic = map_dbl(ssf_modelTO, ~AIC(.)),
         tail_thermal_orographic = map_dbl(ssf_modelTOt, ~AIC(.)),
         tail_thermal = map_dbl(ssf_modelTt, ~AIC(.)),
         thermal_orographic = map_dbl(ssf_modelOt, ~AIC(.)),
         cross = map_dbl(ssf_modelC, ~AIC(.)),
         tail = map_dbl(ssf_modelT, ~AIC(.)),
         orographic = map_dbl(ssf_modelO, ~AIC(.)),
         thermal = map_dbl(ssf_modelt, ~AIC(.)))

# Extract deltaAICs
SSF_AICs <- colSums(SSF_selections[19:33]) # get the summed AICs from all 15 models
SSF_AICs <- as.data.frame(SSF_AICs) %>%
  rownames_to_column() %>% 
  rename(model = rowname,
         AIC = (SSF_AICs)) %>% 
  mutate(deltaAIC = AIC - min(AIC)) %>% # calculate the change in AIC from the best model
  arrange(deltaAIC) # order the df by deltaAIC size

# print those values
library(gridExtra)
pdf("AICs.pdf", height=11, width=10)
grid.table(SSF_AICs)
dev.off()

##########################################################################################

## extract model coefficients for plotting (all years)
##########################################################################################
SSF_results$id_year <- paste(SSF_results$id,SSF_results$year,sep="_")

ssf_coefs4 <- unnest(SSF_results, ssf_coefs4)

tail_coefs4 <- do.call(rbind, lapply(unique(ssf_coefs4$id_year), function(x){head(ssf_coefs4[ssf_coefs4$id_year == x,], 4)[1,]}))
cross_coefs4 <- do.call(rbind, lapply(unique(ssf_coefs4$id_year), function(x){head(ssf_coefs4[ssf_coefs4$id_year == x,], 4)[2,]}))
thermal_coefs4 <- do.call(rbind, lapply(unique(ssf_coefs4$id_year), function(x){head(ssf_coefs4[ssf_coefs4$id_year == x,], 4)[3,]}))
orographic_coefs4 <- do.call(rbind, lapply(unique(ssf_coefs4$id_year), function(x){head(ssf_coefs4[ssf_coefs4$id_year == x,], 4)[4,]}))


tail_coefs4$var <- "tailwind"
cross_coefs4$var <- "crosswind"
thermal_coefs4$var <- "thermal"
orographic_coefs4$var <- "orographic"

plot_data4 <- rbind(tail_coefs4, cross_coefs4,thermal_coefs4,orographic_coefs4)

# add number of each migration to the plot data
plot_data4$Migration <- NA # create an empty column to store values
plot_data4$Migration[which(plot_data4$id_year == "Jaana_2013" | # add year 2 info
                                    plot_data4$id_year == "Lars_2013" | 
                                    plot_data4$id_year == "Mohammed_2018" | 
                                    plot_data4$id_year == "Senta_2015" | 
                                    plot_data4$id_year == "Valentin_2015")] <- "2"
plot_data4$Migration[which(plot_data4$id_year == "Jaana_2014" | # add year 3 info
                                    plot_data4$id_year == "Lars_2014" | 
                                    plot_data4$id_year == "Senta_2016")] <- "3"
plot_data4$Migration[which(plot_data4$id_year == "Lars_2015" | # add year 4 info
                                    plot_data4$id_year == "Senta_2017")] <- "4"
plot_data4$Migration[which(plot_data4$id_year == "Annika_2015" | # add adult info
                             plot_data4$id_year == "Jouko_2015" | 
                             plot_data4$id_year == "Mikko_2015" | 
                             plot_data4$id_year == "Paivi_2017")] <- "5 +" 
plot_data4$Migration[which(is.na(plot_data4$Migration))] <- "1" # everything else is year 1


# repeat for each model
ssf_coefs3 <- unnest(SSF_results, ssf_coefs3)


tail_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id), function(x){head(ssf_coefs3[ssf_coefs3$id == x,], 3)[1,]}))
cross_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id), function(x){head(ssf_coefs3[ssf_coefs3$id == x,], 3)[2,]}))
orographic_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id), function(x){tail(ssf_coefs3[ssf_coefs3$id == x,], 3)[3,]}))


tail_coefs3$var <- "tailwind"
cross_coefs3$var <- "crosswind"
orographic_coefs3$var <- "orographic"

plot_data3 <- rbind(tail_coefs3, cross_coefs3,orographic_coefs3)

# model 2
ssf_coefs2 <- unnest(SSF_results, ssf_coefs2)


tail_coefs2 <- do.call(rbind, lapply(unique(ssf_coefs2$id), function(x){head(ssf_coefs2[ssf_coefs2$id == x,], 2)[1,]}))
cross_coefs2 <- do.call(rbind, lapply(unique(ssf_coefs2$id), function(x){head(ssf_coefs2[ssf_coefs2$id == x,], 2)[2,]}))

tail_coefs2$var <- "tailwind"
cross_coefs2$var <- "crosswind"

plot_data2 <- rbind(tail_coefs2, cross_coefs2)

#############################################################################################


### Plot the wind coefficients
#############################################################################################

## create a scatterplot showing each individual coef/year and connecting paths
library(hrbrthemes)
library(viridis)
library(patchwork)
ggplot(plot_data4, aes(x=Migration, y=ssf_coefs4, color = Migration)) + 
  geom_jitter(aes(fill=Migration), size=2, width=0.2) +
  geom_line(aes(group = id),
            color = "grey",
            arrow = arrow(type = "closed",length=unit(0.075, "inches"))) +
  labs(x="Migration", y = "Coefficients") +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(var))


## create plots of the coefficients for all individuals
violins4 <- ggplot(A2_plot_data4, aes(x=ssf_coefs4, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)#+scale_fill_brewer(palette="Blues")

violins3 <- ggplot(plot_data3, aes(x=ssf_coefs3, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)

violins2 <- ggplot(plot_data2, aes(x=ssf_coefs2, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title = "first autumn migration",x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)

ggarrange(violins4,                                                 # First row with best model
          ggarrange(violins3, violins2, ncol = 2, labels = c("B", "C")), # Second row with other models
          nrow = 2, 
          labels = "A"                                        
) 

## create plots of the coefficients for all individuals in all the years
first_autumn_violins <- ggplot(plot_data4, aes(x=ssf_coefs4, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)
second_autumn_violins <- ggplot(A2_plot_data4, aes(x=ssf_coefs4, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)#+scale_fill_brewer(palette="Blues")
third_autumn_violins <- ggplot(A3_plot_data4, aes(x=ssf_coefs4, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)#+scale_fill_brewer(palette="Blues")
fourth_autumn_violins <- ggplot(A4_plot_data4, aes(x=ssf_coefs4, y=var, color=var, shape=var)) + 
  geom_point(size = 3)+
  #geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=1.2)+
  #geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)#+scale_fill_brewer(palette="Blues")

ggarrange(first_autumn_violins, second_autumn_violins,third_autumn_violins, fourth_autumn_violins, ncol=2, nrow=2,labels = c("1", "2","3","4"), legend = "none")# common.legend = TRUE, legend="right")

## create a single plot with all the individuals for each year
ggplot(plot_data4, aes(x=ssf_coefs4, y=Migration, fill = Migration)) +
  geom_violin(trim=FALSE)+
  #scale_fill_discrete(labels= c("First (n=28)","Second (n=5)","Third (n=3)","Fourth (n=2)","Adult (n=5)"))+
  geom_boxplot(width=0.1, fill="white")+ 
  scale_y_discrete(limits=rev)+
  theme(legend.position  = "none")+
  #xlim(-1,3)+
  labs(x="Coefficients", y = "Migration")+ 
  #scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  geom_vline(xintercept = 0)+ facet_wrap(vars(var))#,rows = vars(Autumn)) 


## create a plot for each covariate
cross4 <- whole_plot_data[which(whole_plot_data$var == "crosswind"),]
cross_plot4 <- ggplot(cross4, aes(x = ssf_coefs4, y = Autumn, fill = Autumn)) +
  geom_boxplot() +
  labs(title = "Crosswind", x="Coefficients", y = "Migration")
tail4 <- whole_plot_data[which(whole_plot_data$var == "tailwind"),]
tail_plot4 <- ggplot(tail4, aes(x = ssf_coefs4, y = Autumn, fill = Autumn)) +
  geom_boxplot() +
  labs(title = "Tailwind", x="Coefficients", y = "Migration")
orog4 <- whole_plot_data[which(whole_plot_data$var == "orographic"),]
orog_plot4 <- ggplot(orog4, aes(x = ssf_coefs4, y = Autumn, fill = Autumn)) +
  geom_boxplot() +
  labs(title = "Orographic", x="Coefficients", y = "Migration")
ther4 <- whole_plot_data[which(whole_plot_data$var == "thermal"),]
ther_plot4 <- ggplot(ther4, aes(x = ssf_coefs4, y = Autumn, fill = Autumn)) +
  geom_boxplot() +
  labs(title = "Thermal", x="Coefficients", y = "Migration")

ggarrange(cross_plot4, tail_plot4, orog_plot4, ther_plot4, ncol=4, nrow=1, legend = "none")# common.legend = TRUE, legend="right")


## create plots of each coefficient for each individual
tailwind_plot <-  ggplot(data = tail_coefs, aes(x = id, y = coefs, colour = id, group = id)) + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)) +
  geom_point(size = 2.5) + ylim(-1,1) +
  #geom_smooth(method=lm, se = F) # Don't add shaded confidence region
  geom_hline(yintercept = 0)


crosswind_plot <-  ggplot(data = cross_coefs, aes(x = id, y = coefs, colour = id, group = id)) + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.position  = "none") +
        geom_point(size = 2.5) +
        ylim(-1,1) +
  geom_hline(yintercept = 0)

ggarrange(tailwind_plot, crosswind_plot, ncol=2, nrow=1, legend = "none")# common.legend = TRUE, legend="right")
#############################################################################################

