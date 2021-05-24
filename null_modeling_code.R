### Null modeling code for Honey Buzzard GPS data
### Hester Brønnvik, Summer 2021


## set working directory and load required packages
setwd("C:/Users/Tess Bronnvik/Desktop/Br-nnvik_honey_buzzard_ssf")
HBNM_packs <- c("lubridate", "amt", "sf", "move", "ggpubr", "RNCEP", "mapview","tidyverse","viridis","hrbrthemes","patchwork", "groupdata2","cvms")
lapply(HBNM_packs, require, character.only = TRUE)

## acquire data
half1 <- read.csv("half_track_100s-8458281616483740304.csv",stringsAsFactors = F, header = T)
half2 <- read.csv("half_track2_100s-7832317087668080045.csv",stringsAsFactors = F, header = T)
annotated_data <- rbind(half1,half2) %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         year = format(as.POSIXct(timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y"),
         id_year = paste(id,year,sep="_"),
         case_ = as.numeric(case_))

colnames(annotated_data)[c(1,2,3,4,16,17,19,20,21)] <- c("lon1","lon2","lat1","lat2","u_wind","v_wind","orographic","thermal","sea_surface")

loc1 <- cbind(annotated_data$lon1,annotated_data$lat1)
loc2 <- cbind(annotated_data$lon2,annotated_data$lat2)
annotated_data$heading <- bearingRhumb(loc1,loc2)
annotated_data$tail <- wind_support(u = annotated_data$u_wind,v = annotated_data$v_wind, heading = annotated_data$heading)
annotated_data$cross <- cross_wind(u = annotated_data$u_wind,v = annotated_data$v_wind, heading = annotated_data$heading)

# find and remove sea crossings
seacrossing_points <- !is.na(annotated_data$sea_surface) # get points over the sea (not NA sea surface)
sum(seacrossing_points) # 186095 NA values
annotated_data_test$id_year_burst_step <- paste(annotated_datat$id,annotated_data$year,
                                         annotated_data$burst_,annotated_data$step_id_,sep = "_") # create a new column of each point for a unique vector
seacrossings <- annotated_data[seacrossing_points == T,] # make a df of the points over seas
annotated_data <- annotated_data %>% # remove all steps containing seacrossings
  filter(!id_year_burst_step %in% seacrossings$id_year_burst_step)


# add a label for each migration to the  data
annotated_data$id_year <- paste(annotated_data$name,annotated_data$yr,sep="_")
annotated_data$Migration <- NA # create an empty column to store values
annotated_data$Migration[which(annotated_data$id_year == "Jaana_2013" | # add year 2 info
                             annotated_data$id_year == "Lars_2013" | 
                             annotated_data$id_year == "Mohammed_2018" | 
                             annotated_data$id_year == "Senta_2015" | 
                             annotated_data$id_year == "Valentin_2015")] <- "2"
annotated_data$Migration[which(annotated_data$id_year == "Jaana_2014" | # add year 4 info
                             annotated_data$id_year == "Lars_2014" | 
                             annotated_data$id_year == "Senta_2016")] <- "3"
annotated_data$Migration[which(annotated_data$id_year == "Lars_2015" | # add year 4 info
                             annotated_data$id_year == "Senta_2017")] <- "4"
annotated_data$Migration[which(annotated_data$id_year == "Annika_2014" | # add adult info
                             annotated_data$id_year == "Annika_2015" | 
                             annotated_data$id_year == "Jouko_2014" | 
                              annotated_data$id_year == "Jouko_2015" | 
                             annotated_data$id_year == "Mikko_2014" | 
                             annotated_data$id_year == "Mikko_2015" | 
                             annotated_data$id_year == "Paivi_2016" | 
                             annotated_data$id_year == "Paivi_2017" |
                             annotated_data$id_year == "Tiina_2014" | 
                             annotated_data$id_year == "Tiina_2015")] <- "5 +" 
annotated_data$Migration[which(is.na(annotated_data$Migration))] <- "1" # everything else is year 1


