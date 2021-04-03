
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
stepLength <- 240 # minutes between bursts
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
                                name == "Jaana"| 
                                  name == "Julia" | name == "Kirsi"| name == "Lars"| name == "Lisa"|
                                name == "Matti"| name == "Miikka" | name == "Mohammed" | name == "Per"| 
                                name == "Piff"| name == "Puff" | name == "Roosa"| 
                                  name == "Rudolf" |
                                  name == "Senta"| 
                                name == "Sven" | name == "Tor"|
                                  name == "Ulla" | name == "Valentin"| name == "Venus" | name == "Viljo")
## all the adults
Autumnsl <- autumns %>% filter(name == "Aarne" | name == "Annika" |name == "Jouko" |name == "Kari" |name == "Mikko" |name == "Paivi")

## the multi-year juveniles
Autumnsl <- autumns %>% filter(name == "Lars" | name == "Mohammed" |name == "Valentin")# |name == "Senta" |name == "Jaana")


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
data_sp <- map_autumns
coordinates(data_sp) <- ~ long + lat
proj4string(data_sp) <- wgs
mapView(data_sp, zcol = "name")

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
write.csv(ssf.df, "28.juv.1A.240.15.csv", row.names = F)

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
ssfdata <- read.csv("28.juv.1A.240.15-1224970117268260264.csv", stringsAsFactors = F)
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
ssfdata$m <- ifelse(ssfdata$year == "2013" & ssfdata$id == "Heidi", "1", NA)

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

colnames(ssfdata)
#ssfdata[c(26:49)] <- NULL
all(complete.cases(ssfdata$order))

## derive the cross and tail wind components
colnames(ssfdata)[c(1,2,3,4,16,17,19,21)] <- c("lon1","lon2","lat1","lat2","u_wind","v_wind","orographic","thermal")


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

# nest the data
ind_by_year <- ssfdata %>% 
  group_by(id,year) %>% 
  nest()

# build the model function
ssf_model <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + scaled_orographic + strata(step_id_), data = df)
}

# model
ind_by_year <- ind_by_year %>% 
  mutate(ssf_model = purrr::map(data, ssf_model))

# apply the coefs
ind_by_year2 <- ind_by_year %>% 
  mutate(ssf_coefs = purrr::map(ssf_model, coef))


# extract AICs
ssf_models <- unnest(ind_by_year2, ssf_model)
ssf_models <- lapply(unique(ssf_models$id), function(x){head(ssf_models[ssf_models$id == x,], 1)}) %>% 
  reduce(rbind)
ssf_models <- ssf_models %>% 
  mutate(ssf_AIC = purrr::map(ssf_model, extractAIC))

## extract model coefficients for plotting
ssf_coefs <- unnest(ind_by_year2, ssf_coefs)


tail_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id), function(x){head(ssf_coefs[ssf_coefs$id == x,], 4)[1,]}))
cross_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id), function(x){head(ssf_coefs[ssf_coefs$id == x,], 4)[2,]}))
thermal_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id), function(x){head(ssf_coefs[ssf_coefs$id == x,], 4)[3,]}))
orographic_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id), function(x){tail(ssf_coefs[ssf_coefs$id == x,], 4)[4,]}))


tail_coefs$var <- "tailwind"
cross_coefs$var <- "crosswind"
thermal_coefs$var <- "thermal"
orographic_coefs$var <- "orographic"

plot_data <- rbind(tail_coefs, cross_coefs,thermal_coefs,orographic_coefs)

#############################################################################################


### Plot the wind coefficients
#############################################################################################

## create plots of the coefficients for all individuals
violins <- ggplot(plot_data, aes(x=ssf_coefs, y=var, fill=var)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title = "first autumn migration",x="Coefficients", y = "Covariate")+ 
  scale_y_discrete(limits=c("tailwind", "crosswind", "thermal", "orographic"))+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0)+
  scale_fill_brewer(palette="Blues")
 

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

