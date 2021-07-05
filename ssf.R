### Code for analyzing route selection in European Honey Buzzards
### Hester Brønnvik, 2021

# set working directory and load required packages
setwd("C:/Users/Tess Bronnvik/Desktop/Br-nnvik_honey_buzzard_ssf")
ssf_packs <- c("lubridate", "amt", "purrr","move", "mapview", "ggpubr", "survival", "sjPlot", "MASS", "tidyverse")
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
toleranceLength <- 15 # tolerance in minutes
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs") # map projection


### Prepare data for Movebank
#############################################################################################
# retrieve the full data set
full_data <- read.csv("./original_data/full_buzz_set.csv",stringsAsFactors = F, header = T)
# reformat the timestamps
full_data$timestamp <- as.POSIXct(strptime(full_data$dt, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")


# select the autumn migrations
all_autumns <- full_data[grep("autumn", full_data$phase),]

# find and remove duplicate observations
doubles <- all_autumns %>% dplyr::select(long, lat, name, timestamp) %>%
  duplicated
sum(doubles) # 59 duplicates
all_autumns <- all_autumns[doubles != TRUE,]

# get the 204 reads when the birds did not move and remove them from the data
resting_reads <- read.csv("resting.csv", stringsAsFactors = FALSE)
all_autumns <- all_autumns %>% mutate(id_ts = paste(name,timestamp, sep="_")) %>% 
  filter(!id_ts %in% resting_reads$id_ts)

# order all the data by timestamp
all_autumns <- all_autumns %>%
  arrange(timestamp)

# create a unique identifier for each migratory journey
all_autumns$id_year <- paste(all_autumns$name, all_autumns$yr, sep="_")

# make simple plots of lat/long to check for outliers
ggplot(all_autumns, aes(x=long, y=lat)) + geom_point()+ facet_wrap(~name, scales="free")
ggplot(all_autumns, aes(x=long, y=lat, color=as.factor(name))) + geom_point()

# make a map of individual tracks
data_sp <- all_autumns # store the data, then create a spatial object for plotting
coordinates(data_sp) <- ~ long + lat # set coordinates
proj4string(data_sp) <- wgs # set projection
mapView(data_sp, zcol = "id_year", burst = F, cex = 3, color = rainbow) # plot on a map


## The data need to be re-sampled at different rates. 
## Separate them into different objects so that each can be processed individually.

# the individuals with 1 hour sampling rates
autumns_one <- all_autumns %>% filter(name == "Anni" | name == "Per" |
                                        name == "Sven"| name == "Viljo" |
                                        name == "Jari" | 
                                        name == "Johannes")
# the individuals with 2 hour sampling rates
autumns_two <- all_autumns %>% filter(name == "Edit" | name == "Emma" | name == "Julia" | name == "Matti"| 
                                        name == "Miikka"| name == "Ulla" | name == "Aida" | 
                                        name == "Ella" | name == "Heidi" | name == "Kirsi" | name == "Gilda" |
                                        name == "Lisa"| name == "Aarne" | name == "Valentin")
# the individuals with 3 hour sampling rates
autumns_three <- all_autumns %>% filter(name == "Mohammed")
# the individuals with 4 hour sampling rates
autumns_four <- all_autumns %>% filter(name == "Annika" | name == "Jaana"| name == "Kari" | name == "Lars"| name == "Piff"| 
                                         name == "Puff" | name == "Roosa"| name == "Senta"| 
                                         name == "Tor" | name == "Tiina" |name == "Jouko" |
                                         name == "Mikko" |name == "Paivi" | name == "Hans" | name == "Venus" |
                                         name == "Rudolf")

# create one hourly tracks
autumn_tracks_1 <- lapply(split(autumns_one, autumns_one$id_year), function(x){
  # make a track object
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs)
  
  # resample to a consistent time between steps
  trk <- track_resample(trk, rate = minutes(60), tolerance = minutes(toleranceLength))
  
  # remove bursts with fewer than 3 steps to allow pythagorean headings
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
}) %>% reduce(rbind)

# create two hourly tracks
autumn_tracks_2 <- lapply(split(autumns_two, autumns_two$id_year), function(x){
  # make a track object
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs)
  
  # resample to a consistent time between steps
  trk <- track_resample(trk, rate = minutes(120), tolerance = minutes(toleranceLength))
  
  # remove bursts with fewer than 3 steps to allow pythagorean headings
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
}) %>% reduce(rbind)

# create three hourly tracks
autumn_tracks_3 <- lapply(split(autumns_three, autumns_three$id_year), function(x){
  # make a track object
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs)
  
  # resample to a consistent time between steps
  trk <- track_resample(trk, rate = minutes(180), tolerance = minutes(toleranceLength))
  
  # remove bursts with fewer than 3 steps to allow pythagorean headings
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
}) %>% reduce(rbind)

# create four hourly tracks
autumn_tracks_4 <- lapply(split(autumns_four, autumns_four$id_year), function(x){
  # make a track object
  trk <- mk_track(x,.x=long, .y=lat, .t=timestamp, id = name, crs = wgs)
  
  # resample to a consistent time between steps
  trk <- track_resample(trk, rate = minutes(240), tolerance = minutes(toleranceLength))
  
  # remove bursts with fewer than 3 steps to allow pythagorean headings
  trk <- filter_min_n_burst(trk, 3)
  
  # burst steps
  burst <- steps_by_burst(trk, keep_cols = "start")
  
  # create random steps using fitted gamma and von mises distributions and append
  rnd_stps <- burst %>% random_steps(n_control = stepNumber)
  
}) %>% reduce(rbind)

# combine those back into one data set
autumn_track <- rbind(autumn_tracks_1,autumn_tracks_2,autumn_tracks_3,autumn_tracks_4)

# plot the step lengths for observed and random steps
obs <- autumn_track[autumn_track$case_ == TRUE,] %>% dplyr::select(id, sl_) %>% 
  ggplot(aes(sl_, fill = factor(id))) + ggtitle("Observed") + geom_density(alpha = 0.4)
rand <- autumn_track[autumn_track$case_ == FALSE,] %>% dplyr::select(id, sl_) %>% 
  ggplot(aes(sl_, fill = factor(id))) + ggtitle("Alternative") + geom_density(alpha = 0.4)
ggarrange(obs, rand, ncol=2, nrow=1, legend = "none")

# calculate headings for each step
autumn_track <- autumn_track %>% rowwise() %>% 
  mutate(heading = NCEP.loxodrome.na(y1_, y2_, x1_, x2_)) %>% 
  ungroup() %>% add_row()

# make the columns meet Movebank's upload expectations
names(autumn_track)[c(2,4)] <-c("location-long", "location-lat") # select the end point for each step
autumn_track$timestamp<-autumn_track$t1_ # select the start time for each step
autumn_track$timestamp <- paste0(autumn_track$timestamp,".000" )

# split the data to conform to Movebank file size limits
half1 <- autumn_track[1:(0.5*nrow(autumn_track)),]
half2 <- autumn_track[(0.5*nrow(autumn_track)):nrow(autumn_track),]

# export each file for Movebank annotation
#write.csv(half1, "name_date.csv", row.names = FALSE)
#write.csv(half1, "name_date.csv", row.names = FALSE)

#############################################################################################

### Prepare data for models
#############################################################################################
# read in the annotated data
half1 <- read.csv("name_date.csv",stringsAsFactors = F, header = T)
half2 <- read.csv("name_date.csv",stringsAsFactors = F, header = T)
ssfdata <- rbind(half1,half2) %>% # reformat columns
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         year = format(as.POSIXct(timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y"),
         case_ = as.numeric(case_))

# check how many true steps survived filtering
print(ssfdata %>% group_by(id_year) %>% tally(case_ == 1), n= 75)

## derive the cross and tail wind components
colnames(ssfdata)[c(1,2,3,4,17,18,20,21)] <- c("lon1","lon2","lat1","lat2","u_wind","v_wind","thermal","sea_surface")

ssfdata$tail <- wind_support(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)
ssfdata$cross <- cross_wind(u = ssfdata$u_wind,v = ssfdata$v_wind, heading = ssfdata$heading)

# find and remove sea crossings
seacrossing_points <- !is.na(ssfdata$sea_surface) # get points over the sea (not NA sea surface)
sum(seacrossing_points)
ssfdata$id_year_burst_step <- paste(ssfdata$id,ssfdata$year,
                                    ssfdata$burst_,ssfdata$step_id_,sep = "_") # create a new column of each point for a unique vector
seacrossings <- ssfdata[seacrossing_points == T,] # make an object of the points over seas
ssfdata <- ssfdata %>% # remove all steps containing seacrossings
  filter(!id_year_burst_step %in% seacrossings$id_year_burst_step)

# map the data removing seacrossings to double-check accuracy
ssfdata_test <- ssfdata %>% distinct(lon1, lat1, .keep_all = TRUE)
coordinates(ssfdata_test) <- ~lon1 + lat1
proj4string(ssfdata_test) <- wgs
mapView(ssfdata_test)

## drop rows with NA in thermal column and then correct for new step numbers
# remove the missing values
thermalNA <- is.na(ssfdata$thermal)
ssfdata <- ssfdata[thermalNA != T,]

# select 101 steps from each step_id_. Choosing the first 101 allows protection of True steps.
ssfdata <- ssfdata %>% 
  group_by(id_year, burst_, step_id_) %>% 
  slice(1:101) %>% 
  ungroup()

# standardize the predictors
ssfdata <- ssfdata %>% 
  mutate(scaled_cross = abs(as.numeric(scale(cross, center = T, scale = T))), # take the absoulte value of crosswind so that it is strength and not direction of the effect
         scaled_tail = as.numeric(scale(tail, center = T, scale = T)),
         scaled_thermal = as.numeric(scale(thermal, center = T, scale = T)),)
#############################################################################################

### Build the model
#############################################################################################
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

# fit separate models for each journey, collect the coefficients, and calculate the AICs
SSF_results <- ssf_modeling_data %>% 
  group_by(id_year) %>% 
  nest() %>% 
  mutate(ssf_modelTCt = purrr::map(data, ssf_modelTCt),
         ssf_coefsTCt = purrr::map(ssf_modelTCt, coef),
         AIC_TCt = map_dbl(ssf_modelTCt, ~AIC(.)))

#############################################################################################

## Extract the model coefficients for plotting
##########################################################################################
# flatten the coefficient column
ssf_coefs <- unnest(SSF_results, ssf_coefsTCt)

# the coefficients for each variable are now in a single column, group them by their identities
tail_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id_year), function(x){head(ssf_coefs[ssf_coefs$id_year == x,], 3)[1,]}))
cross_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id_year), function(x){head(ssf_coefs[ssf_coefs$id_year == x,], 3)[2,]}))
thermal_coefs <- do.call(rbind, lapply(unique(ssf_coefs$id_year), function(x){tail(ssf_coefs[ssf_coefs$id_year == x,], 3)[3,]}))

# create a common label for the variables
tail_coefs$var <- "Wind support"
cross_coefs$var <- "Crosswind"
thermal_coefs$var <- "Thermal uplift"

# bind them into a single object again with their identites appended
plot_data <- rbind(tail_coefs, cross_coefs,thermal_coefs)

## add the life stage and migration number
# separate the id_year column so that stage is identified by individual
plot_data <- separate(plot_data, col = "id_year", c("id", "year"), remove = F)
plot_data$stage <- NA # create an empty column to store values
plot_data$stage[which(plot_data$id == "Aarne" | 
                        plot_data$id == "Annika" |
                        plot_data$id == "Jari" |
                        plot_data$id == "Johannes" |
                        plot_data$id == "Jouko" |
                        plot_data$id == "Mikko" |
                        plot_data$id == "Paivi" |
                        plot_data$id == "Tiina" |
                        plot_data$id == "Kari")] <- "adult" # label individuals of unknown age
plot_data$stage[which(is.na(plot_data$stage))] <- "juvenile" # label individuals tagged as fledglings

# add the number of each migration
plot_data$Migration <- NA 
plot_data$Migration[which(plot_data$id_year == "Jaana_2013" | # add year 2 info
                            plot_data$id_year == "Lars_2013" | 
                            plot_data$id_year == "Mohammed_2018" | 
                            plot_data$id_year == "Senta_2015" | 
                            plot_data$id_year == "Valentin_2015" | 
                            plot_data$id_year == "Aarne_2013" |
                            plot_data$id_year == "Annika_2012" |
                            plot_data$id_year == "Jari_2013" | 
                            plot_data$id_year == "Johannes_2014" | 
                            plot_data$id_year == "Jouko_2012" |
                            plot_data$id_year == "Kari_2012" |
                            plot_data$id_year == "Paivi_2014" |
                            plot_data$id_year == "Mikko_2012" |
                            plot_data$id_year == "Tiina_2012")] <- "2"

plot_data$Migration[which(plot_data$id_year == "Jaana_2014" | # add year 3 info
                            plot_data$id_year == "Lars_2014" | 
                            plot_data$id_year == "Senta_2016" |
                            plot_data$id_year == "Annika_2013" | 
                            plot_data$id_year == "Jari_2014" | 
                            plot_data$id_year == "Johannes_2015" |
                            plot_data$id_year == "Jouko_2013" |
                            plot_data$id_year == "Mikko_2013" | # missing
                            plot_data$id_year == "Paivi_2015" |
                            plot_data$id_year == "Mohammed_2019" |
                            plot_data$id_year == "Tiina_2013")] <- "3"

plot_data$Migration[which(plot_data$id_year == "Lars_2015" | # add year 4 info
                            plot_data$id_year == "Senta_2017" |
                            plot_data$id_year == "Annika_2014" |
                            plot_data$id_year == "Jouko_2014" |
                            plot_data$id_year == "Mikko_2014" |
                            plot_data$id_year == "Paivi_2016" |
                            plot_data$id_year == "Tiina_2014")] <- "4"

plot_data$Migration[which(plot_data$id_year == "Annika_2015" | # add year 5 info
                            plot_data$id_year == "Jouko_2015" | 
                            plot_data$id_year == "Mikko_2015" | 
                            plot_data$id_year == "Paivi_2017" | 
                            plot_data$id_year == "Tiina_2015")] <- "5" 
plot_data$Migration[which(is.na(plot_data$Migration))] <- "1" # everything else is year 1
#############################################################################################

### Plot the coefficients
#############################################################################################
## create a scatterplot showing each individual coef/year and connecting paths

# select colors that are distinct
id_colors <- c("#542344","#1478A3","#00A354","#F34213","#FFBF00")

# first, select the data from the stage and variable of interest, 
# then select the data from individuals that survived in multiple years, 
# remove the multiple year individuals from the set containing others, 
# plot a box plot for all the individuals, layer on points for birds surviving only one year, 
# layer on points for birds surviving in multiple years, connect those points with lines, 
# and finally plot all in panels with separate axes
tail_plot_data <- plot_data %>% filter(var == "Wind support" & stage == "juvenile")
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
cross_plot_data <- plot_data %>% filter(var == "Crosswind" & stage == "juvenile")
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
heat_plot_data <- plot_data %>% filter(var == "Thermal uplift" & stage == "juvenile")
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

ggarrange(tail_plot, cross_plot, heat_plot, ncol=3, nrow=1, legend = "none")
#############################################################################################

### Model selection
#############################################################################################
## For all-possible-regressions model selection, build 6 more model functions to cover all combinations of variables

# three with two variables
ssf_modelTC <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_tail + strata(step_id_), data = df)
}
ssf_modelTt <- function(df) {
  fit_clogit(case_ ~ scaled_tail + scaled_thermal + strata(step_id_), data = df)
}
ssf_modelCt <- function(df) {
  fit_clogit(case_ ~ scaled_cross + scaled_thermal + strata(step_id_), data = df)
}

# and three with single variables
ssf_modelC <- function(df) {
  fit_clogit(case_ ~ scaled_cross + strata(step_id_), data = df)
}
ssf_modelT <- function(df) {
  fit_clogit(case_ ~ scaled_tail + strata(step_id_), data = df)
}
ssf_modelt <- function(df) {
  fit_clogit(case_ ~ scaled_thermal + strata(step_id_), data = df)
}

# model the data again using these models
SSF_selections <- ssf_modeling_data %>% 
  group_by(id,year) %>% 
  nest() %>% 
  mutate(ssf_modelTCt = purrr::map(data, ssf_modelTCt), # add each model as a column of lists 
         ssf_modelTC = purrr::map(data, ssf_modelTC),
         ssf_modelCt = purrr::map(data, ssf_modelCt),
         ssf_modelTt = purrr::map(data, ssf_modelTt),
         ssf_modelC = purrr::map(data, ssf_modelC),
         ssf_modelT = purrr::map(data, ssf_modelT),
         ssf_modelt = purrr::map(data, ssf_modelt),
         tail_cross_thermal = map_dbl(ssf_modelTCt, ~AIC(.)), # extract AIC for each model
         tail_cross = map_dbl(ssf_modelTC, ~AIC(.)),
         cross_thermal = map_dbl(ssf_modelCt, ~AIC(.)),
         tail_thermal = map_dbl(ssf_modelTt, ~AIC(.)),
         cross = map_dbl(ssf_modelC, ~AIC(.)),
         tail = map_dbl(ssf_modelT, ~AIC(.)),
         thermal = map_dbl(ssf_modelt, ~AIC(.)))


## Extract deltaAICs

# assign migrations so that AICs can be taken for each
## add the life stage and migration number
SSF_selections$id_year <- paste(SSF_selections$id,SSF_selections$year,sep="_")
SSF_selections$stage <- NA # create an empty column to store values
SSF_selections$stage[which(SSF_selections$id == "Aarne" | 
                             SSF_selections$id == "Annika" |
                             SSF_selections$id == "Jari" |
                             SSF_selections$id == "Johannes" |
                             SSF_selections$id == "Jouko" |
                             SSF_selections$id == "Mikko" |
                             SSF_selections$id == "Paivi" |
                             SSF_selections$id == "Tiina" |
                             SSF_selections$id == "Kari")] <- "adult" # label individuals of unknown age
SSF_selections$stage[which(is.na(SSF_selections$stage))] <- "juvenile" # label individuals tagged as fledglings

# add the number of each migration
SSF_selections$Migration <- NA 
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
                                 SSF_selections$id_year == "Kari_2012" |
                                 SSF_selections$id_year == "Paivi_2014" |
                                 SSF_selections$id_year == "Mikko_2012" |
                                 SSF_selections$id_year == "Tiina_2012")] <- "2"

SSF_selections$Migration[which(SSF_selections$id_year == "Jaana_2014" | # add year 3 info
                                 SSF_selections$id_year == "Lars_2014" | 
                                 SSF_selections$id_year == "Senta_2016" |
                                 SSF_selections$id_year == "Annika_2013" | 
                                 SSF_selections$id_year == "Jari_2014" | 
                                 SSF_selections$id_year == "Johannes_2015" |
                                 SSF_selections$id_year == "Jouko_2013" |
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

SSF_selections$Migration[which(SSF_selections$id_year == "Annika_2015" | # add year 5 info
                                 SSF_selections$id_year == "Jouko_2015" | 
                                 SSF_selections$id_year == "Mikko_2015" | 
                                 SSF_selections$id_year == "Paivi_2017" | 
                                 SSF_selections$id_year == "Tiina_2015")] <- "5" 
SSF_selections$Migration[which(is.na(SSF_selections$Migration))] <- "1" # everything else is year 1

## calculate the average AIC for each migration
# begin by selecting which data to use (in this case excluding results for models of individuals of unknown age)
SSF_AICs <- SSF_selections %>% filter(stage == "juvenile" | stage == "adult" & Migration == "5")

# separate the migrations and calculate the average for each
SSF_AICs <- lapply(split(SSF_selections, SSF_selections$Migration), function(x){
  colMeans(x[11:16])
}) %>% 
  reduce(rbind) %>% 
  t() 
colnames(SSF_AICs) <- c("first_migration", "second_migration", "third_migration", "fourth_migration","adult_migration")

# set the order to display models in
model_order <- c("tail_cross_thermal", "tail_cross", "tail_thermal", "cross_thermal", "tail", "cross", "thermal")

# calculate the deltaAICs and order
SSF_AICs <- as.data.frame(SSF_AICs) %>%
  rownames_to_column() %>% 
  rename(variables = rowname) %>% 
  mutate(deltaAIC_1 = first_migration - min(first_migration),# calculate the change in AIC from the best model for each migration
         deltaAIC_2 = second_migration - min(second_migration),
         deltaAIC_3 = third_migration - min(third_migration),
         deltaAIC_4 = fourth_migration - min(fourth_migration),
         deltaAIC_5 = adult_migration - min(adult_migration))

# export the table as a .doc file
tab_df(SSF_AICs,
       alternate.rows = T, # this colors the rows
       title = "Model selection",
       file = "SSF_AICs.doc")

## For stepwise selection, separate the migrations and run each

## again, add the life stage and migration number
ssf_modeling_data$stage <- NA # create an empty column to store values
ssf_modeling_data$stage[which(ssf_modeling_data$id == "Aarne" | 
                                ssf_modeling_data$id == "Annika" |
                                ssf_modeling_data$id == "Jari" |
                                ssf_modeling_data$id == "Johannes" |
                                ssf_modeling_data$id == "Jouko" |
                                ssf_modeling_data$id == "Mikko" |
                                ssf_modeling_data$id == "Paivi" |
                                ssf_modeling_data$id == "Tiina" |
                                ssf_modeling_data$id == "Kari")] <- "adult" # label individuals of unknown age
ssf_modeling_data$stage[which(is.na(ssf_modeling_data$stage))] <- "juvenile" # label individuals tagged as fledglings

# add the number of each migration
ssf_modeling_data$Migration <- NA 
ssf_modeling_data$Migration[which(ssf_modeling_data$id_year == "Jaana_2013" | # add year 2 info
                                    ssf_modeling_data$id_year == "Lars_2013" | 
                                    ssf_modeling_data$id_year == "Mohammed_2018" | 
                                    ssf_modeling_data$id_year == "Senta_2015" | 
                                    ssf_modeling_data$id_year == "Valentin_2015" | 
                                    ssf_modeling_data$id_year == "Aarne_2013" |
                                    ssf_modeling_data$id_year == "Annika_2012" |
                                    ssf_modeling_data$id_year == "Jari_2013" | 
                                    ssf_modeling_data$id_year == "Johannes_2014" | 
                                    ssf_modeling_data$id_year == "Jouko_2012" |
                                    ssf_modeling_data$id_year == "Kari_2012" |
                                    ssf_modeling_data$id_year == "Paivi_2014" |
                                    ssf_modeling_data$id_year == "Mikko_2012" |
                                    ssf_modeling_data$id_year == "Tiina_2012")] <- "2"

ssf_modeling_data$Migration[which(ssf_modeling_data$id_year == "Jaana_2014" | # add year 3 info
                                    ssf_modeling_data$id_year == "Lars_2014" | 
                                    ssf_modeling_data$id_year == "Senta_2016" |
                                    ssf_modeling_data$id_year == "Annika_2013" | 
                                    ssf_modeling_data$id_year == "Jari_2014" | 
                                    ssf_modeling_data$id_year == "Johannes_2015" |
                                    ssf_modeling_data$id_year == "Jouko_2013" |
                                    ssf_modeling_data$id_year == "Mikko_2013" | # missing
                                    ssf_modeling_data$id_year == "Paivi_2015" |
                                    ssf_modeling_data$id_year == "Mohammed_2019" |
                                    ssf_modeling_data$id_year == "Tiina_2013")] <- "3"

ssf_modeling_data$Migration[which(ssf_modeling_data$id_year == "Lars_2015" | # add year 4 info
                                    ssf_modeling_data$id_year == "Senta_2017" |
                                    ssf_modeling_data$id_year == "Annika_2014" |
                                    ssf_modeling_data$id_year == "Jouko_2014" |
                                    ssf_modeling_data$id_year == "Mikko_2014" |
                                    ssf_modeling_data$id_year == "Paivi_2016" |
                                    ssf_modeling_data$id_year == "Tiina_2014")] <- "4"

ssf_modeling_data$Migration[which(ssf_modeling_data$id_year == "Annika_2015" | # add year 5 info
                                    ssf_modeling_data$id_year == "Jouko_2015" | 
                                    ssf_modeling_data$id_year == "Mikko_2015" | 
                                    ssf_modeling_data$id_year == "Paivi_2017" | 
                                    ssf_modeling_data$id_year == "Tiina_2015")] <- "5" 
ssf_modeling_data$Migration[which(is.na(ssf_modeling_data$Migration))] <- "1" # everything else is year 1

# select which data to use (in this case excluding individuals of unknown age)
ssf_modeling_data <- ssf_modeling_data %>% filter(stage == "adult" & Migration == "5" | stage == "juvenile") 

# using only data from each migration, model with all variables, then run stepwise selection
first_data <- ssf_modeling_data %>% filter(Migration == "1") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
first_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_), data = first_data)
stepAIC(first_model)

second_data <- ssf_modeling_data %>% filter(Migration == "2") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
second_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_), data = second_data)
stepAIC(second_model)

third_data <- ssf_modeling_data %>% filter(Migration == "3") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
third_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_), data = third_data)
stepAIC(third_model)

fourth_data <- ssf_modeling_data %>% filter(Migration == "4") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_)
fourth_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal  + strata(step_id_), data = fourth_data)
stepAIC(fourth_model)

fifth_data <- ssf_modeling_data %>% filter(Migration == "5") %>% dplyr::select(case_,scaled_tail,scaled_cross,scaled_thermal,step_id_) %>% filter(!is.na(scaled_tail))
fifth_model <- clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal  + strata(step_id_), data = fifth_data)
stepAIC(fifth_model)
##########################################################################################

### Cross-validation
#############################################################################################
# prune the data set to reduce its size
validation_data <- ssf_modeling_data %>% dplyr::select(id,case_,year,burst_,step_id_,id_year_burst_step,scaled_tail,scaled_cross,scaled_thermal)

# group the data by step so that partitioning happens on steps rather than slicing them
validation_data <- validation_data %>% group_by(id_year_burst_step) %>% 
  mutate(id_grouping = as.factor(id))

# partition data in an 80:20 split
data_partitioned <- partition(
  validation_data,
  p = 0.8,
  cat_col = "case_",
  list_out = T)


# fit the model to the training set and then calculate RMSE
valid_tibble <- validate(
  data_partitioned,
  formulas = "case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_)",
  family = "binomial", # use binary classification 
  REML = FALSE)

# save the result for later because it takes a long time to generate it
#save(valid_tibble, file ="validation_result.RData")
#############################################################################################
