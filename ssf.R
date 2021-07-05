
# set working directory and load required packages
setwd("C:/Users/Tess Bronnvik/Desktop/Br-nnvik_honey_buzzard_ssf")
ssf_packs <- c("lubridate", "amt", "sf", "move", "dplyr", "purrr", "ggpubr", "ggplot2", "RNCEP", "mapview","tidyverse","viridis","hrbrthemes","patchwork", "groupdata2","cvms", "survival")
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
source("wind_support_Kami.R")

## set criteria for tracks
stepNumber <- 150 # random steps
stepLength <- 120 # minutes between bursts
toleranceLength <- 15 # tolerance in minutes
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")


### Prepare data for Movebank
#############################################################################################
full_data <- read.csv("./original_data/full_buzz_set.csv",stringsAsFactors = F, header = T)
full_data$timestamp <- as.POSIXct(strptime(full_data$dt, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")


# get only autumn migrations
all_autumns <- full_data[grep("autumn", full_data$phase),]

# find and remove duplicate observations
ind2 <- all_autumns %>% dplyr::select(long, lat, name, timestamp) %>%
  duplicated
sum(ind2) # 59 duplicates
all_autumns <- all_autumns[ind2!=TRUE,]

# order all the data by timestamp
all_autumns <- all_autumns %>% #group_by(name) %>% 
  arrange(timestamp)#, .by_group = T) %>%
  #ungroup()

# get the 178 locations when the birds did not move and remove them from the data
resting <- read.csv("resting.csv", stringsAsFactors = F)
all_autumns2 <- all_autumns %>% mutate(id_ts = paste(name,timestamp, sep="_")) %>% 
  filter(!id_ts %in% resting$id_ts)

## the individuals with 1 hour median sampling rate
Autumnsl <- all_autumns %>% filter(name == "Anni" | name == "Edit" | name == "Emma" | name == "Julia" |
                                     name == "Matti"| name == "Miikka" | name == "Per"| 
                                     name == "Sven" | name == "Ulla" | name == "Viljo" |
                                     name == "Jari" | name == "Johannes")
## the individuals with 2 hour median sampling rate
Autumnsl <- all_autumns %>% filter(name == "Aida" | name == "Ella" | name == "Heidi" | name == "Kirsi"| 
                                     name == "Lisa"| name == "Rudolf" | name == "Valentin" | 
                                     name == "Aarne" | name == "Annika" | name == "Venus" | name == "Gilda") # Annika, Jouko, and Venus throw errors in the spring
## the individuals with 3 hour median sampling rate
Autumnsl <- all_autumns %>% filter(name == "Mohammed")
## the individuals with 4 hour median sampling rate
Autumnsl <- all_autumns %>% filter(name == "Jaana"| name == "Kari" | name == "Lars"| name == "Piff"| 
                                     name == "Puff" | name == "Roosa"| name == "Senta"| 
                                     name == "Tor" | name == "Tiina" |name == "Jouko" |
                                     name == "Mikko" |name == "Paivi" | name == "Hans")

## for mapping by migration, add info
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
thiraut$Migration <- "fifth"
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
data_sp <- autumn_track[autumn_track$case_ == 1,]
coordinates(data_sp) <- ~ x1_ + y1_
proj4string(data_sp) <- wgs
mapView(data_sp, zcol = "id", burst = F, cex = 3, color = rainbow)
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


#trk <- mk_track(Autumns, .x=long, .y=lat, .t=timestamp, id = name, crs = CRS("+init=epsg:4326"))

#transform to geographic coordinates
#trk <- transform_coords(trk, CRS("+init=epsg:3857")) # pseudo-Mercator in meters (https://epsg.io/3857)


## make the data of class track
stepLength <- 240

Autumnsl$id_year <- paste(Autumnsl$name,Autumnsl$yr, sep="_")
unique(Autumnsl$id_year)


autumn_track <- lapply(split(Autumnsl, Autumnsl$id_year), function(x){
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs)
  
  # resample to a consistent time between steps
  trk <- track_resample(trk, rate = minutes(stepLength), tolerance = minutes(toleranceLength))
  
  # remove bursts with fewer than 3 steps to allow pythagorean headings
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
}) %>% reduce(rbind)
autumn_track <- as_tibble(autumn_track)
head(autumn_track)

# find steps of length 0
autumn_track$difflon <- autumn_track$x2_ - autumn_track$x1_
autumn_track$difflat <- autumn_track$y2_ - autumn_track$y1_
no_movement <- autumn_track[which(autumn_track$difflon == 0 & autumn_track$difflat ==0),]# & autumn_track$case_ == T),]
no_movement$step_id_ <- NA
no_movement$case_ <- NA
no_movement$id_ts <- paste(no_movement$id,no_movement$t1_,sep="_")

loc1 <- cbind(autumn_track$x1_,autumn_track$y1_)
loc2 <- cbind(autumn_track$x2_,autumn_track$y2_)
autumn_track$heading <- bearingRhumb(loc1,loc2)
autumn_track$dist <- distRhumb(loc1,loc2)

no_move <- rbind(no_move, no_movement)
write.csv(no_move, file = "resting.csv", row.names = F)


# save tracks from each resampling rate
one.hours <- autumn_track
two.hours <- autumn_track
three.hours <- autumn_track
four.hours <- autumn_track

# combine those back into one data set
autumn_track <- rbind(one.hours,two.hours,three.hours,four.hours)

autumn_track <- autumn_track %>% rowwise() %>% 
  mutate(heading = NCEP.loxodrome.na(y1_, y2_, x1_, x2_)) %>% 
  ungroup()

half1 <- autumn_track[1:765645,]
half2 <- autumn_track[765646:nrow(autumn_track),]

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

# find times when the bird did not move
autumn_track$difflon <- autumn_track$x2_ - autumn_track$x1_
autumn_track$difflat <- autumn_track$y2_ - autumn_track$y1_
no_movement <- autumn_track[which(autumn_track$difflon == 0 & autumn_track$difflat == 0),]
Hour <- as.POSIXlt(no_movement$t2_)$hour
hist(Hour, breaks=seq(0, 23), main="Time (hour)", xlim =c(0,23))

autumn_track2 <- SpatialPointsDataFrame(autumn_track[,c("x2_","y2_")], autumn_track, proj4string = CRS("+init=epsg:3857"))  
ssf.df <- half1#data.frame(spTransform(autumn_track2, CRS("+proj=longlat +datum=WGS84"))) 
names(ssf.df)[c(2,4)] <-c("location-long", "location-lat")
ssf.df$timestamp<-ssf.df$t1_
ssf.df$timestamp <- paste0(ssf.df$timestamp,".000" )
write.csv(ssf.df, "half_track1_june_150s.csv", row.names = F)

#############################################################################################

### Prepare data for models
#############################################################################################
half1 <- read.csv("half_track1_june_150s-4433398109146553161.csv",stringsAsFactors = F, header = T)
half2 <- read.csv("half_track2_june_150s-4573417134914709132.csv",stringsAsFactors = F, header = T)
ssfdata <- rbind(half1,half2) %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         year = format(as.POSIXct(timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y"),
         id_year = paste(id,year,sep="_"),
         case_ = as.numeric(case_))

# check how many true steps survived filtering
print(ssfdata %>% group_by(id_year) %>% tally(case_ == 1), n= 75)


## derive the cross and tail wind components
colnames(ssfdata)[c(1,2,3,4,17,18,20,21)] <- c("lon1","lon2","lat1","lat2","u_wind","v_wind","thermal","sea_surface")


ssfdata$tail <- wind_support(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)
ssfdata$cross <- cross_wind(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)
ssfdata_saved <- ssfdata



# find and remove sea crossings
seacrossing_points <- !is.na(ssfdata$sea_surface) # get points over the sea (not NA sea surface)
sum(seacrossing_points) # 226235 NA values
ssfdata$id_year_burst_step <- paste(ssfdata$id,ssfdata$year,
                                         ssfdata$burst_,ssfdata$step_id_,sep = "_") # create a new column of each point for a unique vector
seacrossings <- ssfdata[seacrossing_points == T,] # make a df of the points over seas
ssfdata <- ssfdata %>% # remove all steps containing seacrossings
  filter(!id_year_burst_step %in% seacrossings$id_year_burst_step)


# map the data removing seacrossings to double-check accuracy
ssfdata_test <- ssfdata %>% distinct(lon1, lat1, .keep_all = TRUE)
coordinates(ssfdata_test) <- ~lon1 + lat1
proj4string(ssfdata_test) <- wgs
mapView(ssfdata_test)

# find and plot NA values in thermal uplift
thermalNA <- is.na(ssfdata$thermal)
sum(thermalNA) # 27487 NA values
missing_thermals <- ssfdata[thermalNA == T,] # select the observations missing values

coordinates(missing_thermals) <- ~ lon1 + lat1 # prepare for a map
proj4string(missing_thermals) <- wgs
mapView(missing_thermals, zcol = "id")

# which individual has the most missing values
sum(with(missing_thermals,id == "Mohammed"))# Mohammed has 27,482 missing values in 2019 & 2020
sum(with(missing_thermals,id == "Emma")) # 1
sum(with(missing_thermals,id == "Jouko")) # 1
sum(with(missing_thermals,id == "Jaana")) # 1
sum(with(missing_thermals,id == "Lars")) # 2

## drop rows with NA in thermal column and then correct for new step numbers

# remove the missing values
thermalNA <- is.na(ssfdata$thermal)
ssfdata <- ssfdata[thermalNA != T,]

# select data with missing values
Edit_ssfdata <- ssfdata[ssfdata$id == "Edit",] # burst 1 step 31 (True) 


#remove those data from the data set
ssfdata <- ssfdata %>% filter(id != "Edit")# & id != "Annika")

# edit each ID with missing data to have the right number of steps remaining
Edit_ssfdata <- Edit_ssfdata[which(Edit_ssfdata$step_id_ != 31),] # remove step with 101 NAs

# add the IDs back to the data
ssfdata <- rbind(ssfdata,Edit_ssfdata)

  
# select 101 steps from each step_id_. Choosing the first 101 allows protection of True steps.
ssfdata <- ssfdata %>% 
  #mutate(id_year = paste(id,year,sep = "_")) %>% 
  group_by(id_year, burst_, step_id_) %>% 
  slice(1:101) %>% 
  ungroup()


### normalize the data
ssfdata <- ssfdata %>% 
  #group_by(id) %>% 
  mutate(scaled_cross = abs(as.numeric(scale(cross, center = T, scale = T))),
         scaled_tail = as.numeric(scale(tail, center = T, scale = T)),
         scaled_thermal = as.numeric(scale(thermal, center = T, scale = T)),
         #scaled_orographic = as.numeric(scale(orographic, center = T, scale = T))
         ) %>%
  ungroup()
#############################################################################################
data_sf <- ssfdata %>% st_as_sf(coords = c("lon2", "lat2"), crs = wgs)


### Build the models
#############################################################################################
# load the prepped data
load("ssfdata_prepped.RData")

# build the model function
ssf_modelTCt <- function(df) {
  clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_), data = df)
}

# check how many true steps survived filtering
print(ssfdata %>% group_by(id_year) %>% tally(case_ == 1), n= 75)

# exclude certain data
ssf_modeling_data <- ssfdata %>% filter(id_year != "Johannes_2015" & # only three steps
                                        id_year != "Jouko_2015" & # only one step
                                        id_year != "Mohammed_2019") # incomplete record


SSF_results <- ssf_modeling_data %>% 
  group_by(id_year) %>% 
  nest() %>% 
  mutate(ssf_modelTCt = purrr::map(data, ssf_modelTCt),
         ssf_coefsTCt = purrr::map(ssf_modelTCt, coef),
         AIC_TCt = map_dbl(ssf_modelTCt, ~AIC(.)))

#############################################################################################


### Model selection
#############################################################################################
# build the model function
ssf_modelTCti <- function(df) {
  clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_tail*scaled_thermal + strata(step_id_), data = df)
}

ssf_modelTC <- function(df) {
  clogit(case_ ~ scaled_tail + scaled_cross + strata(step_id_), data = df)
}


# For stepwise model selection, 
# nest the data, then apply above models to each individual for each year and extract AIC
SSF_selections <- ssf_modeling_data %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_modelTCti = purrr::map(data, ssf_modelTCti), # add each model as a column of lists 
         ssf_modelTCt = purrr::map(data, ssf_modelTCt),
         ssf_modelTC = purrr::map(data, ssf_modelTC),
         ssf_coefsTCti = purrr::map(ssf_modelTCti, coef), # extract the coefficients for each model
         ssf_coefsTCt = purrr::map(ssf_modelTCt, coef),
         ssf_coefsTC = purrr::map(ssf_modelTC, coef),
         AIC_TCti = map_dbl(ssf_modelTCti, ~AIC(.)),
         AIC_TCt = map_dbl(ssf_modelTCt, ~AIC(.)),
         AIC_TC = map_dbl(ssf_modelTC, ~AIC(.)))

# to get stepwise regression on each id_year the model has to be built as clogit
ssf_aic_data <- ssf_modeling_data %>% dplyr::select(case_,scaled_tail, scaled_cross, scaled_thermal, ta_, sl_,id,year, id_year_burst_step,step_id_) %>% 
  filter(!is.na(ta_)) %>% filter(!is.na(scaled_tail))
ssf_modelTCtts <- function(df) {clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + ta_ + sl_ + strata(id_year_burst_step), data = ssf_aic_data)}
stepAIC(ssf_modelTCtts())


# automated stepwise
library(MASS)
ssf_modeling_data <- ssf_modeling_data %>% filter(stage == "adult" & Migration == "5" | stage == "juvenile") 

first_data <- ssf_modeling_data %>% filter(Migration == "1") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
first_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_tail*scaled_thermal + strata(step_id_), data = first_data)
stepAIC(first_model)

second_data <- ssf_modeling_data %>% filter(Migration == "2") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
second_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_tail*scaled_thermal + strata(step_id_), data = second_data)
stepAIC(second_model)

third_data <- ssf_modeling_data %>% filter(Migration == "3") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
third_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_tail*scaled_thermal + strata(step_id_), data = third_data)
stepAIC(fifth_model)

fourth_data <- ssf_modeling_data %>% filter(Migration == "4") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
fourth_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_tail*scaled_thermal + strata(step_id_), data = fourth_data)
stepAIC(fourth_model)

fifth_data <- ssf_modeling_data %>% filter(Migration == "5") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_) %>% filter(!is.na(scaled_tail))
fifth_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_tail*scaled_thermal + strata(step_id_), data = fifth_data)
stepAIC(fifth_model)

select_result <- data.frame()
uid <- unique(ssf_aic_data$id_year)
for (i in 1:length(uid)) {
  temp_df <- ssf_aic_data[which(ssf_aic_data$id_year == i),]
  ssf_modelTCtts <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + ta_ + sl_ + strata(step_id_), data = temp_df)
  selections <- stepAIC(ssf_modelTCtts())
  temp_result <- rbind(select_result, selections)
}

ssf_modelTCtts <- function(df) {clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + ta_ + sl_ + strata(step_id_), data = Senta)}
  
# For all-possible-regressions model selection, build 16 models to cover all combinations of variables
ssf_modelCO <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_orographic + strata(step_id_), data = df)
}
ssf_modelCOt <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_orographic + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelCt <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelCTO <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_tail + scaled_orographic + strata(step_id_), data = df)
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

ssfdata <- ssfdata %>% filter(id_year != "Jaana_2014" & id_year != "Senta_2017")
SSF_selections <- ssfdata %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_modelTCtO = purrr::map(data, ssf_modelTCtO), # add each model as a column of lists 
         ssf_modelTCt = purrr::map(data, ssf_modelTCt),
         ssf_modelCOt = purrr::map(data, ssf_modelCOt),
         ssf_modelCTO = purrr::map(data, ssf_modelCTO),
         ssf_modelTOt = purrr::map(data, ssf_modelTOt),
         ssf_modelTC = purrr::map(data, ssf_modelTC),
         ssf_modelTO = purrr::map(data, ssf_modelTO),
         ssf_modelCt = purrr::map(data, ssf_modelCt),
         ssf_modelTt = purrr::map(data, ssf_modelTt),
         ssf_modelOt = purrr::map(data, ssf_modelOt),
         ssf_modelCO = purrr::map(data, ssf_modelCO),
         ssf_modelC = purrr::map(data, ssf_modelC),
         ssf_modelT = purrr::map(data, ssf_modelT),
         ssf_modelO = purrr::map(data, ssf_modelO),
         ssf_modelt = purrr::map(data, ssf_modelt),
         tail_cross_thermal_orographic = map_dbl(ssf_modelTCtO, ~AIC(.)), # extract AIC for each model
         tail_cross_thermal = map_dbl(ssf_modelTCt, ~AIC(.)),
         cross_thermal_orographic = map_dbl(ssf_modelCOt, ~AIC(.)),
         tail_cross_orographic = map_dbl(ssf_modelCTO, ~AIC(.)),
         tail_thermal_orographic = map_dbl(ssf_modelTOt, ~AIC(.)),
         tail_cross = map_dbl(ssf_modelTC, ~AIC(.)),
         cross_orographic = map_dbl(ssf_modelCO, ~AIC(.)),
         cross_thermal = map_dbl(ssf_modelCt, ~AIC(.)),
         tail_orographic = map_dbl(ssf_modelTO, ~AIC(.)),
         tail_thermal = map_dbl(ssf_modelTt, ~AIC(.)),
         thermal_orographic = map_dbl(ssf_modelOt, ~AIC(.)),
         cross = map_dbl(ssf_modelC, ~AIC(.)),
         tail = map_dbl(ssf_modelT, ~AIC(.)),
         orographic = map_dbl(ssf_modelO, ~AIC(.)),
         thermal = map_dbl(ssf_modelt, ~AIC(.)))

# for only the all-four model and each one removing only one variable at a time
SSF_selections_Round1 <- ssfdata %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_modelTCtO = purrr::map(data, ssf_modelTCtO), # add each model as a column of lists 
         ssf_modelTCt = purrr::map(data, ssf_modelTCt),
         ssf_modelCOt = purrr::map(data, ssf_modelCOt),
         ssf_modelCTO = purrr::map(data, ssf_modelCTO),
         ssf_modelTOt = purrr::map(data, ssf_modelTOt),
         tail_cross_thermal_orographic = map_dbl(ssf_modelTCtO, ~AIC(.)), # extract AIC for each model
         tail_cross_thermal = map_dbl(ssf_modelTCt, ~AIC(.)),
         cross_thermal_orographic = map_dbl(ssf_modelCOt, ~AIC(.)),
         tail_cross_orographic = map_dbl(ssf_modelCTO, ~AIC(.)),
         tail_thermal_orographic = map_dbl(ssf_modelTOt, ~AIC(.)))

# looking at each model with only two variables left from the model removing orographic uplift
SSF_selections_Round2 <- ssfdata %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_modelTC = purrr::map(data, ssf_modelTC),
         ssf_modelCt = purrr::map(data, ssf_modelCt),
         ssf_modelTt = purrr::map(data, ssf_modelTt),
         tail_cross = map_dbl(ssf_modelTC, ~AIC(.)),
         cross_thermal = map_dbl(ssf_modelCt, ~AIC(.)),
         tail_thermal = map_dbl(ssf_modelTt, ~AIC(.)))

## Extract deltaAICs
# assign migrations so that AICs can be taken for each
SSF_selections$id_year <- paste(SSF_selections$id,SSF_selections$year,sep="_")
SSF_selections$Migration <- NA # create an empty column to store values
SSF_selections$Migration[which(SSF_selections$id_year == "Jaana_2013" | # add year 2 info
                                 SSF_selections$id_year == "Lars_2013" | 
                                 SSF_selections$id_year == "Mohammed_2018" | 
                                 SSF_selections$id_year == "Senta_2015" | 
                                 SSF_selections$id_year == "Valentin_2015" | 
                                 SSF_selections$id_year == "Aarne_2013" |
                                 SSF_selections$id_year == "Annika_2012" | 
                                 SSF_selections$id_year == "Jari_2013" | 
                                 SSF_selections$id_year == "Johannes_2014" | 
                                 SSF_selections$id_year == "Jouko_2012" |
                                 SSF_selections$id_year == "Paivi_2014" |
                                 SSF_selections$id_year == "Mikko_2012" |
                                 SSF_selections$id_year == "Tiina_2012")] <- "2"

SSF_selections$Migration[which(SSF_selections$id_year == "Jaana_2014" | # add year 3 info
                                 SSF_selections$id_year == "Lars_2014" | 
                                 SSF_selections$id_year == "Senta_2016" |
                                 SSF_selections$id_year == "Annika_2013" | 
                                 SSF_selections$id_year == "Jari_2014" | 
                                 SSF_selections$id_year == "Johannes_2015" |
                                 SSF_selections$id_year == "Mikko_2013" | # missing
                                 SSF_selections$id_year == "Paivi_2015" |
                                 SSF_selections$id_year == "Mohammed_2019" |
                                 SSF_selections$id_year == "Tiina_2013")] <- "3"

SSF_selections$Migration[which(SSF_selections$id_year == "Lars_2015" | # add year 4 info
                                 SSF_selections$id_year == "Senta_2017" |
                                 SSF_selections$id_year == "Annika_2014" |
                                 SSF_selections$id_year == "Jouko_2014" |
                                 SSF_selections$id_year == "Mikko_2014" |
                                 SSF_selections$id_year == "Paivi_2016" |
                                 SSF_selections$id_year == "Tiina_2014")] <- "4"

SSF_selections$Migration[which(SSF_selections$id_year == "Annika_2015" | # add adult infoSSF_selections$id_year == "Annika_2015" |
                                 SSF_selections$id_year == "Jouko_2015" | 
                                 SSF_selections$id_year == "Mikko_2015" | 
                                 SSF_selections$id_year == "Paivi_2017" | 
                                 SSF_selections$id_year == "Tiina_2015")] <- "5" 
SSF_selections$Migration[which(is.na(SSF_selections$Migration))] <- "1" # everything else is year 1

SSF_selections$stage <- NA
SSF_selections$stage[which(SSF_selections$id == "Aarne" | 
                             SSF_selections$id == "Annika" |
                             SSF_selections$id == "Jari" |
                             SSF_selections$id == "Johannes" |
                             SSF_selections$id == "Jouko" |
                             SSF_selections$id == "Mikko" |
                             SSF_selections$id == "Paivi" |
                             SSF_selections$id == "Tiina" |
                             SSF_selections$id == "Kari")] <- "adult"
SSF_selections$stage[which(is.na(SSF_selections$stage))] <- "juvenile"
# get the average AIC for each migration
SSF_selections5 <- SSF_selections %>% filter(stage == "juvenile" | stage == "adult" & Migration == "5")

SSF_AICs <- lapply(split(SSF_selections5, SSF_selections5$Migration), function(x){
  colMeans(x[10:12])
}) %>% 
  reduce(rbind) %>% 
  t() 
colnames(SSF_AICs) <- c("first_migration", "second_migration", "fifth_migration", "fourth_migration","adult_migration")

# set the order to display models in
model_order <- c("tail_cross_thermal_orographic", "tail_cross_thermal", "cross_thermal_orographic",
                 "tail_cross_orographic", "tail_thermal_orographic","tail_cross","cross_orographic",
                 "cross_thermal","tail_orographic", "tail_thermal", "thermal_orographic","cross",
                 "tail","thermal","orographic")
model_order <- c("tail_cross_thermal_tailXthermal","tail_cross_thermal", "tail_cross")

# calculate the deltaAICs and order
SSF_AICs <- as.data.frame(SSF_AICs) %>%
  rownames_to_column() %>% 
  rename(variables = rowname) %>% 
  mutate(deltaAIC_1 = first_migration - min(first_migration),# calculate the change in AIC from the best model for each migration
         deltaAIC_2 = second_migration - min(second_migration),
         deltaAIC_3 = fifth_migration - min(fifth_migration),
         deltaAIC_4 = fourth_migration - min(fourth_migration),
         deltaAIC_5 = adult_migration - min(adult_migration))#,
         #variables =  factor(variables, levels = model_order)) %>%
  #arrange(variables)# order the df by variables size

library(sjPlot)
tab_df(SSF_AICs,
       alternate.rows = T, # this colors the rows
       title = "Model selection",
       file = "SSF_AICs_inter.doc")

#separate the migrations again in order to arrange each by best model
SSF_AICs_1 <- SSF_AICs %>% 
  dplyr::select(variables,first_migration,deltaAIC_1) %>% 
  arrange(deltaAIC_1)

SSF_AICs_2 <- SSF_AICs %>% 
  dplyr::select(variables,second_migration,deltaAIC_2) %>% 
  arrange(deltaAIC_2)

SSF_AICs_3 <- SSF_AICs %>% 
  dplyr::select(variables,fifth_migration,deltaAIC_3) %>% 
  arrange(deltaAIC_3)

SSF_AICs_4 <- SSF_AICs %>% 
  dplyr::select(variables,fourth_migration,deltaAIC_4) %>% 
  arrange(deltaAIC_4)

SSF_AICs_5 <- SSF_AICs %>% 
  dplyr::select(variables,adult_migration,deltaAIC_5) %>% 
  arrange(deltaAIC_5)

allaicstab <- cbind(SSF_AICs_1,SSF_AICs_2,SSF_AICs_3,SSF_AICs_4,SSF_AICs_5)

# print those values
library(gridExtra)
pdf("AICs_1.pdf", height=11, width=10,onefile = F)
grid.table(SSF_AICs_1)
dev.off()

pdf("AICs_2.pdf", height=11, width=10,onefile = F)
grid.table(SSF_AICs_2)
dev.off()

pdf("AICs_3.pdf", height=11, width=10,onefile = F)
grid.table(SSF_AICs_3)
dev.off()

pdf("AICs_4.pdf", height=11, width=10,onefile = F)
grid.table(SSF_AICs_4)
dev.off()

pdf("AICs_5.pdf", height=11, width=10,onefile = F)
grid.table(SSF_AICs_5)
dev.off()

library(pdftools)
pdf_combine(c("AICs_1.pdf", "AICs_2.pdf", "AICs_3.pdf", "AICs_4.pdf", "AICs_5.pdf"), output = "AICs.pdf")

##########################################################################################

### Cross-validation
#############################################################################################

# Partition data to list of data frames
# The first data frame will be train (80% of the data)
# The second will be test (20% of the data)

load("validation_result.RData") # a tibble from only 30 of the steps due to RAM

library(doParallel)
registerDoParallel(4)

validation_data <- ssfdata %>% select(id,case_,year,burst_,step_id_,id_year_burst_step,scaled_tail,scaled_cross,scaled_thermal)

validation_data <- validation_data %>% group_by(id_year_burst_step) %>% 
  mutate(id_grouping = as.factor(id))

data_partitioned <- partition(
  validation_data,
  p = 0.8,
  cat_col = "case_",
  list_out = T
)


# Validate a model
# Gaussian
valid_tibble <- validate(
  data_partitioned,
  formulas = "case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_)",
  family = "gaussian",
  REML = FALSE
)
save(valid_tibble, file ="validation_result.RData")


ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE)
mod_fit <- train(Class ~ Age + ForeignWorker + Property.RealEstate + Housing.Own + 
                   CreditHistory.Critical,  data=GermanCredit, method="glm", family="binomial",
                 trControl = ctrl, tuneLength = 5)
pred = predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$Class)
#############################################################################################


## Extract model coefficients for plotting (all years)
##########################################################################################
ssf_coefs3 <- unnest(SSF_results, ssf_coefsTCt)


tail_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id_year), function(x){head(ssf_coefs3[ssf_coefs3$id_year == x,], 4)[1,]}))
cross_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id_year), function(x){head(ssf_coefs3[ssf_coefs3$id_year == x,], 4)[2,]}))
thermal_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id_year), function(x){tail(ssf_coefs3[ssf_coefs3$id_year == x,], 4)[3,]}))

wind_heat_coefs <- do.call(rbind, lapply(unique(ssf_coefs3$id_year), function(x){tail(ssf_coefs3[ssf_coefs3$id_year == x,], 4)[4,]}))

ta_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id_year), function(x){tail(ssf_coefs3[ssf_coefs3$id_year == x,], 5)[4,]}))
sl_coefs3 <- do.call(rbind, lapply(unique(ssf_coefs3$id_year), function(x){tail(ssf_coefs3[ssf_coefs3$id_year == x,], 5)[5,]}))


tail_coefs3$var <- "Wind support"
cross_coefs3$var <- "Crosswind"
thermal_coefs3$var <- "Thermal uplift"

wind_heat_coefs$var <- "TailxLift"

ta_coefs3$var <- "Turn angle"
sl_coefs3$var <- "Step length"

# add stage to the data
plot_data3 <- rbind(tail_coefs3, cross_coefs3,thermal_coefs3)#, wind_heat_coefs)#,ta_coefs3,sl_coefs3)

plot_data3 <- plot_data3 %>% separate(id_year, c("id","year"), sep = "_")
plot_data3$id_year <- paste(plot_data3$id,plot_data3$year,sep="_")
plot_data3$stage <- NA
plot_data3$stage[which(plot_data3$id == "Aarne" | 
                           plot_data3$id == "Annika" |
                           plot_data3$id == "Jari" |
                           plot_data3$id == "Johannes" |
                           plot_data3$id == "Jouko" |
                           plot_data3$id == "Mikko" |
                           plot_data3$id == "Paivi" |
                           plot_data3$id == "Tiina" |
                           plot_data3$id == "Kari")] <- "adult"
plot_data3$stage[which(is.na(plot_data3$stage))] <- "juvenile"

# add number of each migration to the plot data
plot_data3$Migration <- NA # create an empty column to store values
plot_data3$Migration[which(plot_data3$id_year == "Jaana_2013" | # add year 2 info
                                 plot_data3$id_year == "Lars_2013" | 
                                 plot_data3$id_year == "Mohammed_2018" | 
                                 plot_data3$id_year == "Senta_2015" | 
                                 plot_data3$id_year == "Valentin_2015" | 
                                 plot_data3$id_year == "Aarne_2013" |
                                 plot_data3$id_year == "Annika_2012" |
                                 plot_data3$id_year == "Jari_2013" | 
                                 plot_data3$id_year == "Johannes_2014" | 
                                 plot_data3$id_year == "Jouko_2012" |
                                 plot_data3$id_year == "Kari_2012" |
                                 plot_data3$id_year == "Paivi_2014" |
                                 plot_data3$id_year == "Mikko_2012" |
                                 plot_data3$id_year == "Tiina_2012")] <- "2"

plot_data3$Migration[which(plot_data3$id_year == "Jaana_2014" | # add year 3 info
                                 plot_data3$id_year == "Lars_2014" | 
                                 plot_data3$id_year == "Senta_2016" |
                                 plot_data3$id_year == "Annika_2013" | 
                                 plot_data3$id_year == "Jari_2014" | 
                                 plot_data3$id_year == "Johannes_2015" |
                                 plot_data3$id_year == "Jouko_2013" |
                                 plot_data3$id_year == "Mikko_2013" | # missing
                                 plot_data3$id_year == "Paivi_2015" |
                                 plot_data3$id_year == "Mohammed_2019" |
                                 plot_data3$id_year == "Tiina_2013")] <- "3"

plot_data3$Migration[which(plot_data3$id_year == "Lars_2015" | # add year 4 info
                                 plot_data3$id_year == "Senta_2017" |
                                 plot_data3$id_year == "Annika_2014" |
                                 plot_data3$id_year == "Jouko_2014" |
                                 plot_data3$id_year == "Mikko_2014" |
                                 plot_data3$id_year == "Paivi_2016" |
                                 plot_data3$id_year == "Tiina_2014")] <- "4"

plot_data3$Migration[which(plot_data3$id_year == "Annika_2015" | # add adult infoplot_data3$id_year == "Annika_2015" |
                                 plot_data3$id_year == "Jouko_2015" | 
                                 plot_data3$id_year == "Mikko_2015" | 
                                 plot_data3$id_year == "Paivi_2017" | 
                                 plot_data3$id_year == "Tiina_2015")] <- "5" 
plot_data3$Migration[which(is.na(plot_data3$Migration))] <- "1" # everything else is year 1
#############################################################################################


### Plot the wind coefficients
#############################################################################################
## create a scatterplot showing each individual coef/year and connecting paths
id_colors <- c("#542344","#1478A3","#00A354","#F34213","#FFBF00")

geom_boxplot2 <- function(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge2", 
                          ..., outlier.colour = NA, outlier.color = NA, outlier.fill = NA, 
                          outlier.shape = NA, outlier.size = NA, outlier.stroke = NA, 
                          outlier.alpha = NA, notch = FALSE, notchwidth = 0.5, varwidth = FALSE, 
                          na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                          linetype = "dashed"){
  list(
    geom_boxplot(mapping = mapping, data = data, stat = stat, position = position,
                 outlier.colour = outlier.colour, outlier.color = outlier.color, 
                 outlier.fill = outlier.fill, outlier.shape = outlier.shape, 
                 outlier.size = outlier.size, outlier.stroke = outlier.stroke, 
                 outlier.alpha = outlier.alpha, notch = notch, 
                 notchwidth = notchwidth, varwidth = varwidth, na.rm = na.rm, 
                 show.legend = show.legend, inherit.aes = inherit.aes, 
                 linetype = linetype, ...),
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) ,
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) ,
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) ,
    theme_classic(), # remove panel background and gridlines
    theme(plot.title = element_text(hjust = 0.5,  # hjust = 0.5 centers the title
                                    size = 14,
                                    face = "bold"),
          panel.border = element_rect(linetype = "solid",
                                      colour = "black", fill = "NA", size = 0.5)))}

tail_plot_data <- plot_data3 %>% filter(var == "Wind support" & stage == "juvenile")
tail_plot_highlight <- tail_plot_data %>% filter(id == "Jaana" | id == "Senta" | id == "Valentin" | id == "Mohammed" | id == "Lars")
tail_plot_data_points <- tail_plot_data %>% filter(id != "Jaana" & id != "Senta" & id != "Valentin" & id != "Mohammed" & id != "Lars")
tail_plot <- ggplot(tail_plot_data, aes(x=Migration, y=ssf_coefsTCt)) + 
  geom_boxplot(outlier.shape = NA, color = "grey30", fill = "grey80") +
  geom_jitter(data = tail_plot_data_points, aes(), size=2, width=0.2, color = "grey8") +
  geom_point(data = tail_plot_highlight, aes(x=Migration, y=ssf_coefsTCt, fill = id), shape = 21, size = 3, color = "black") +
  geom_line(data = tail_plot_highlight, aes(group = id, color =id)) +
  scale_fill_manual(values = id_colors) + 
  scale_color_manual(values = id_colors) + 
  labs(title = "Wind support", x="Migration", y = "Coefficients" ) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position = "none")
cross_plot_data <- plot_data3 %>% filter(var == "Crosswind" & stage == "juvenile")
cross_plot_highlight <- cross_plot_data %>% filter(id == "Jaana" | id == "Senta" | id == "Valentin" | id == "Mohammed" | id == "Lars")
cross_plot_data_points <- cross_plot_data %>% filter(id != "Jaana" & id != "Senta" & id != "Valentin" & id != "Mohammed" & id != "Lars")
cross_plot <- ggplot(cross_plot_data, aes(x=Migration, y=ssf_coefsTCt)) + 
  geom_boxplot(outlier.shape = NA, color = "grey30", fill = "grey80") +
  geom_jitter(data = cross_plot_data_points, aes(), size=2, width=0.2, color = "grey8") +
  geom_point(data = cross_plot_highlight, aes(x=Migration, y=ssf_coefsTCt, fill = id), shape = 21, size = 3, color = "black")  + 
  geom_line(data = cross_plot_highlight, aes(group = id, color =id)) + 
  scale_fill_manual(values = id_colors) + 
  scale_color_manual(values = id_colors) + 
  labs(title = "Crosswind", x="Migration", y = "Coefficients" ) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position = "none")
heat_plot_data <- plot_data3 %>% filter(var == "Thermal uplift" & stage == "juvenile")
heat_plot_highlight <- heat_plot_data %>% filter(id == "Jaana" | id == "Senta" | id == "Valentin" | id == "Mohammed" | id == "Lars")
heat_plot_data_points <- heat_plot_data %>% filter(id != "Jaana" & id != "Senta" & id != "Valentin" & id != "Mohammed" & id != "Lars")
heat_plot <- ggplot(heat_plot_data, aes(x=Migration, y=ssf_coefsTCt)) + 
  geom_boxplot(outlier.shape = NA, color = "grey30", fill = "grey80") +
  geom_jitter(data = heat_plot_data_points, aes(), size=2, width=0.2, color = "grey8") +
  geom_point(data = heat_plot_highlight, aes(x=Migration, y=ssf_coefsTCt, fill = id), shape = 21, size = 3, color = "black")  + 
  geom_line(data = heat_plot_highlight, aes(group = id, color =id)) + 
  scale_fill_manual(values = id_colors) + 
  scale_color_manual(values = id_colors) + 
  labs(title = "Thermal uplift", x="Migration", y = "Coefficients" ) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position = "none")

ggarrange(tail_plot, cross_plot, heat_plot, inter_plot, ncol=2, nrow=2, legend = "none")# common.legend = TRUE, legend="right")

inter_plot_data <- plot_data3 %>% filter(var == "TailxLift" & stage == "juvenile")
inter_plot_highlight <- inter_plot_data %>% filter(id == "Jaana" | id == "Senta" | id == "Valentin" | id == "Mohammed" | id == "Lars")
inter_plot_data_points <- inter_plot_data %>% filter(id != "Jaana" & id != "Senta" & id != "Valentin" & id != "Mohammed" & id != "Lars")
inter_plot <- ggplot(inter_plot_data, aes(x=Migration, y=ssf_coefsTCt)) + 
  geom_boxplot(outlier.shape = NA, color = "grey30", fill = "grey80") +
  geom_jitter(data = inter_plot_data_points, aes(), size=2, width=0.2, color = "grey8") +
  geom_point(data = inter_plot_highlight, aes(x=Migration, y=ssf_coefsTCt, fill = id), shape = 21, size = 3, color = "black")  + 
  geom_line(data = inter_plot_highlight, aes(group = id, color =id)) + 
  scale_fill_manual(values = id_colors) + 
  scale_color_manual(values = id_colors) + 
  labs(title = "Wind support x uplift", x="Migration", y = "Coefficients" ) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position = "none")

ta_plot_data <- plot_data3[plot_data3$var == "Turn angle" & plot_data3$stage == "juvenile",]
ta_dots <- ggplot(ta_plot_data, aes(x=Migration, y=ssf_coefsTCt, color = id)) + 
  stat_summary(fun=median, geom="point", shape=18,size=3, color="black") +
  geom_jitter(aes(fill=id), size=2, width=0.2) +
  geom_line(aes(group = id, color =id)) +
  labs(title = "Turn angle", x="Migration", y = "Coefficients") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") #+facet_wrap(vars(var))
sl_plot_data <- plot_data3[plot_data3$var == "Step length" & plot_data3$stage == "juvenile",]
sl_dots <- ggplot(sl_plot_data, aes(x=Migration, y=ssf_coefsTCt, color = id)) + 
  stat_summary(fun=median, geom="point", shape=18,size=3, color="black") +
  geom_jitter(aes(fill=id), size=2, width=0.2) +
  geom_line(aes(group = id, color =id)) +
  labs(title = "Step length", x="Migration", y = "Coefficients") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") #+facet_wrap(vars(var))

## Density plots
ggplot(plot_data3[which(plot_data3$var == "Wind support"),], aes(x =  ssf_coefsTCt, group = Migration, fill = Migration)) + 
  geom_density(alpha = 0.4)+ 
  scale_fill_manual(values = id_colors) +
  labs(title = "Wind support", x="Coefficient", y = "Density")

supp <- as.data.frame(ssf_modeling_data$tail)
colnames(supp)[1] <- "value"
supp$variable <- "tail"

cross <- as.data.frame(ssf_modeling_data$cross)
colnames(cross)[1] <- "value"
cross$variable <- "cross"

heat <- as.data.frame(ssf_modeling_data$thermal)
colnames(heat)[1] <- "value"
heat$variable <- "thermal"

dens <- rbind(supp,cross)#,heat)

ggplot(dens, aes(x =  value, group = variable, fill = variable)) + 
  geom_density(alpha = 0.4)+ 
  scale_fill_manual(values = id_colors)

## Violin plots
violin_tail <- ggplot(tail_plot_data, aes(x=Migration, y=ssf_coefsTCt)) + 
  geom_violin(trim=FALSE, fill="forestgreen")+
  geom_boxplot(width=0.1, fill="white")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_line(aes(group = id, color =id)) +
  labs(title = "Wind support", x="", y = "Coefficients") +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette="YlGn")
violin_cross <- ggplot(cross_plot_data, aes(x=Migration, y=ssf_coefsTCt, fill = Migration)) + 
  geom_violin(trim=FALSE)+#, fill= "firebrick3")+
  geom_boxplot(width=0.1, fill="white")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_line(aes(group = id, color =id)) +
  labs(title = "Crosswind", x="Migration", y = "")+ 
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette="OrRd")
violin_thermal <- ggplot(thermal_plot_data, aes(x=Migration, y=ssf_coefsTCt)) + 
  geom_violin(trim=FALSE, fill="cornflowerblue")+
  geom_boxplot(width=0.1, fill="white")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_line(aes(group = id, color =id)) +
  labs(title = "Thermal uplift", x="", y = "")+
  geom_hline(yintercept = 0)+ 
  scale_fill_brewer(palette="PuBu")

ggarrange(violin_tail, violin_cross, violin_thermal,ncol=3, nrow=1, legend = "none")

# the same plot but faceted so that the axes are common
ggplot(plot_data3, aes(x=Migration, y=ssf_coefsTCt, color = id)) + 
  stat_summary(fun=median, geom="point", shape=18,size=3, color="black") +
  geom_jitter(aes(fill=id), size=2, width=0.2) +
  geom_line(aes(group = id, color =id)) +
  labs( x="Migration", y = "Coefficients") +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none") +
  facet_wrap(vars(var))

ggplot(heat_plot_data, aes(x = Migration, y = ssf_coefsTCt, colour = Variable, group = Variable )) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(id))

## create plots of the coefficients for all individuals
violins4 <- ggplot(A2_plot_data4, aes(x=ssf_coefs4, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)#+scale_fill_brewer(palette="Blues")

violins3 <- ggplot(plot_data3, aes(x=ssf_coefsTCt, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Coefficients", y = "Covariate")+ 
  #scale_y_discrete(limits=c("tailwind", "crosswind", "orographic"))+
  #theme(legend.position = "none")+
  geom_vline(xintercept = 0) +
  facet_wrap(vars(var), scales = "free")

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
fifth_autumn_violins <- ggplot(A3_plot_data4, aes(x=ssf_coefs4, y=var, fill=var)) + 
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

ggarrange(first_autumn_violins, second_autumn_violins,fifth_autumn_violins, fourth_autumn_violins, ncol=2, nrow=2,labels = c("1", "2","3","4"), legend = "none")# common.legend = TRUE, legend="right")

## create a single plot with all the individuals for each year
ggplot(plot_data4, aes(x=ssf_coefs4, y=Migration, fill = Migration)) +
  geom_violin(trim=FALSE)+
  #scale_fill_discrete(labels= c("First (n=28)","Second (n=5)","fifth (n=3)","Fourth (n=2)","Adult (n=5)"))+
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

