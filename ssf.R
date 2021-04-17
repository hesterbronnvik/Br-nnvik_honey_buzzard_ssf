
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

juveniles_experienced <- read.csv("5.juv.2A.29s.240.15-7017657342414562020.csv", stringsAsFactors = F)
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

adults <- read.csv("6.adu.A.120.15-7433769662941172309.csv", stringsAsFactors = F)
adults$timestamp <- as.POSIXct(strptime(adults$timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
adults$year <- format(as.POSIXct(adults$timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y")
colnames(adults)[8] <- "id"
adults$id_year <- paste(adults$id,adults$year,sep="_")
adults$Migration <- "Adult"


ssfdata <- rbind(juveniles,juveniles_experienced)#,adults)
ssfdata$case_ <- as.numeric(ssfdata$case_)


# remove Lisa and Sven, the outliers in coefs
ssfdata <- ssfdata %>% filter(id != "Sven" & id != "Lisa" & id != "Kirsi")

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
sum(thermalNA) # 3896 NA values
missing_thermals <- ssfdata[thermalNA!=F,] # select the observations missing values
coordinates(missing_thermals) <- ~ lon1 + lat1 # prepare for a map
proj4string(missing_thermals) <- wgs
mapView(missing_thermals, zcol = "id")


summary(ssfdata) # 118 NAs in heading

# drop rows with NA in thermal column to see how that affects the result
thermalNA <- is.na(ssfdata$thermal)
sum(thermalNA) # 3896 NA values
ssfdata <- ssfdata[thermalNA != T,]

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


# nest the data
SSF_results <- ssfdata %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_model4 = purrr::map(data, ssf_model4), # add each model as a column of lists 
         ssf_model3 = purrr::map(data, ssf_model3),
         ssf_model2 = purrr::map(data, ssf_model2),
         ssf_coefs4 = purrr::map(ssf_model4, coef), # extract the coefficients for each model
         ssf_coefs3 = purrr::map(ssf_model3, coef),
         ssf_coefs2 = purrr::map(ssf_model2, coef))



# extract AICs (extractAIC does not work on a list, so it cannot be applied above) for  multi-year data, this returns only the second year
##########################################################################################
ssf_models <- unnest(SSF_results, ssf_model4)
ssf_models <- do.call(rbind, lapply(unique(ssf_models$id), function(x){head(ssf_models[ssf_models$id == x,], 1)}))
ssf_models <- ssf_models %>% 
  mutate(ssf_AIC = purrr::map(ssf_model4, extractAIC))
ssf_AICs <- unnest(ssf_models, ssf_AIC)
ssf_AICs <- do.call(rbind, lapply(unique(ssf_AICs$id), function(x){tail(ssf_AICs[ssf_AICs$id == x,], 1)}))
ssf_AICs <- ssf_AICs %>% select(id,year,ssf_AIC)
head(ssf_AICs)

model4_AICs <- ssf_AICs
colnames(model4_AICs)[3] <- "AIC_4"

ssf_models <- unnest(SSF_results, ssf_model3)
ssf_models <- do.call(rbind, lapply(unique(ssf_models$id), function(x){head(ssf_models[ssf_models$id == x,], 1)}))
ssf_models <- ssf_models %>% 
  mutate(ssf_AIC = purrr::map(ssf_model3, extractAIC))
ssf_AICs <- unnest(ssf_models, ssf_AIC)
ssf_AICs <- do.call(rbind, lapply(unique(ssf_AICs$id), function(x){tail(ssf_AICs[ssf_AICs$id == x,], 1)}))
ssf_AICs <- ssf_AICs %>% select(id,year,ssf_AIC)

model3_AICs <- ssf_AICs


ssf_models <- unnest(SSF_results, ssf_model2)
ssf_models <- do.call(rbind, lapply(unique(ssf_models$id), function(x){head(ssf_models[ssf_models$id == x,], 1)}))
ssf_models <- ssf_models %>% 
  mutate(ssf_AIC = purrr::map(ssf_model2, extractAIC))
ssf_AICs <- unnest(ssf_models, ssf_AIC)
ssf_AICs <- do.call(rbind, lapply(unique(ssf_AICs$id), function(x){tail(ssf_AICs[ssf_AICs$id == x,], 1)}))
ssf_AICs <- ssf_AICs %>% select(id,year,ssf_AIC)

model2_AICs <- ssf_AICs

SSF_AICs <- cbind(model4_AICs,AIC_3 = (model3_AICs$ssf_AIC), AIC_2 = (model2_AICs$ssf_AIC))
head(SSF_AICs)

## extract AICs for the third year
ssf_models <- unnest(SSF_results, ssf_model4)
ssf_models <- ssf_models %>% filter(id != "Mohammed" & id != "Valentin")
ssf_models <- do.call(rbind, lapply(unique(ssf_models$id), function(x){tail(ssf_models[ssf_models$id == x,], 1)}))
ssf_models <- ssf_models %>% 
  mutate(ssf_AIC = purrr::map(ssf_model4, extractAIC))
ssf_AICs <- unnest(ssf_models, ssf_AIC)
ssf_AICs <- do.call(rbind, lapply(unique(ssf_AICs$id), function(x){tail(ssf_AICs[ssf_AICs$id == x,], 1)}))
ssf_AICs <- ssf_AICs %>% select(id,year,ssf_AIC)

model4_AICs <- ssf_AICs
model3_AICs <- ssf_AICs
model2_AICs <- ssf_AICs

colnames(model4_AICs)[3] <- "AIC_4"

SSF_AICs <- cbind(model4_AICs,AIC_3 = (model3_AICs$ssf_AIC), AIC_2 = (model2_AICs$ssf_AIC))
head(SSF_AICs)
##########################################################################################

## extract model coefficients for plotting (all years)
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
#plot_data4$Autumn <- "1"

A2_plot_data4 <- plot_data4[which(plot_data4$id_year == "Jaana_2013" | 
                                    plot_data4$id_year == "Lars_2013" | 
                                    plot_data4$id_year == "Mohammed_2018" | 
                                    plot_data4$id_year == "Senta_2015" | 
                                    plot_data4$id_year == "Valentin_2015"),]
A2_plot_data4$Autumn <- "2"
A3_plot_data4 <- plot_data4[which(plot_data4$id_year == "Jaana_2014" | 
                                    plot_data4$id_year == "Lars_2014" | 
                                    plot_data4$id_year == "Senta_2016"),]
A3_plot_data4$Autumn <- "3"
A4_plot_data4 <- plot_data4[which(plot_data4$id_year == "Lars_2015" |
                                    plot_data4$id_year == "Senta_2017"),]
A4_plot_data4$Autumn <- "4"

#plot_data_A1 <- plot_data4
# if plot_data_A1 stores the data from the first autumn and plot_data4 adults this gives all the coefs for model 4
#whole_plot_data <- rbind(plot_data_A1, A2_plot_data4,A3_plot_data4,A4_plot_data4, plot_data4)
#whole_plot_data <- whole_plot_data %>% select(id,year,ssf_coefs4,id_year,var,Autumn)
#whole_plot_data$Autumn <- factor(whole_plot_data$Autumn, levels = c("1","2","3","4"))

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
ggplot(whole_plot_data, aes(x=ssf_coefs4, y=Autumn, fill = Autumn)) +
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

