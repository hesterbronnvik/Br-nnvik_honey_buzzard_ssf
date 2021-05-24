### Null modeling code for Honey Buzzard GPS data
### Hester Brønnvik, Summer 2021


## set working directory and load required packages
setwd("C:/Users/Tess Bronnvik/Desktop/Br-nnvik_honey_buzzard_ssf")
HBNM_packs <- c("lubridate", "modelr", "sf", "move", "ggpubr", "RNCEP", "mapview","tidyverse","viridis","hrbrthemes","patchwork", "groupdata2","cvms")
lapply(HBNM_packs, require, character.only = TRUE)

# import required functions
source("wind_support_Kami.R")

# define values
no.perm <- 100


### Prepare data
############################################################################################
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
#annotated_data <- annotated_data %>% filter(!is.na(heading)) # remove locations where the bird did not move, 181 are true cases
annotated_data$tail <- wind_support(u = annotated_data$u_wind,v = annotated_data$v_wind, heading = annotated_data$heading)
annotated_data$cross <- cross_wind(u = annotated_data$u_wind,v = annotated_data$v_wind, heading = annotated_data$heading)

# find and remove sea crossings
seacrossing_points <- !is.na(annotated_data$sea_surface) # get points over the sea (not NA sea surface)
sum(seacrossing_points) # 191559 NA values
annotated_data$id_year_burst_step <- paste(annotated_data$id,annotated_data$year,
                                         annotated_data$burst_,annotated_data$step_id_,sep = "_") # create a new column of each point for a unique vector
seacrossings <- annotated_data[seacrossing_points == T,] # make a df of the points over seas
annotated_data <- annotated_data %>% # remove all steps containing seacrossings
  filter(!id_year_burst_step %in% seacrossings$id_year_burst_step)

# remove the rows missing thermal uplift values
thermalNA <- is.na(annotated_data$thermal)
sum(thermalNA) # 20,206
annotated_data <- annotated_data[thermalNA != T,]

# select data with missing values
Edit_annotated_data <- annotated_data[annotated_data$id == "Edit",] # burst 1 step 31 (True) 

#remove those data from the data set
annotated_data <- annotated_data %>% filter(id != "Edit")

# edit each ID with missing data to have the right number of steps remaining
Edit_annotated_data <- Edit_annotated_data[which(Edit_annotated_data$step_id_ != 31),] # remove step with 101 NAs

# add the IDs back to the data
annotated_data <- rbind(annotated_data,Edit_annotated_data)
rm(Edit_annotated_data)

# add a label for each migration to the  data
annotated_data$Migration <- NA # create an empty column to store values
annotated_data$Migration[which(annotated_data$id_year == "Jaana_2013" | # add year 2 info
                             annotated_data$id_year == "Lars_2013" | 
                             annotated_data$id_year == "Mohammed_2018" | 
                             annotated_data$id_year == "Senta_2015" | 
                             annotated_data$id_year == "Valentin_2015")] <- "2"
annotated_data$Migration[which(annotated_data$id_year == "Jaana_2014" | # add year 3 info
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


# select 50 steps from each step_id_. Choosing the first 50 allows protection of True steps.
annotated_data <- annotated_data %>% 
  mutate(id_year = paste(id,year,sep = "_")) %>% 
  group_by(id_year, burst_, step_id_) %>% 
  slice(1:51) %>% 
  ungroup()


### normalize the data
annotated_data <- annotated_data %>%  
  mutate(scaled_cross = abs(as.numeric(scale(cross, center = T, scale = T))),
         scaled_tail = as.numeric(scale(tail, center = T, scale = T)),
         scaled_thermal = as.numeric(scale(thermal, center = T, scale = T)),
         scaled_orographic = as.numeric(scale(orographic, center = T, scale = T))) %>%
  ungroup()
save(annotated_data, file = "prepped_data.RData")
############################################################################################

### Model
############################################################################################
load("prepped_data.RData")

## fit a model to each group
# build the model function
modelTCt <- function(df) {
  clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_), data = df)
}

# exclude certain data
modeling_data <- annotated_data %>% filter(#id != "Jouko" & # these never seem to converge
                                          id_year != "Jouko_2011" & id_year != "Jouko_2012" & id_year != "Jouko_2013" & 
                                          id_year != "Mikko_2011" & id_year != "Mikko_2012" & id_year != "Mikko_2013"&
                                          id_year != "Annika_2011" & id_year != "Annika_2012" & id_year != "Annika_2013" & 
                                          id_year != "Mikko_2011" & id_year != "Mikko_2012"  & id_year != "Mikko_2013" &
                                          id_year != "Paivi_2013" & id_year != "Paivi_2014"  & id_year != "Paivi_2015" & 
                                          id != "Jari" & id != "Aarne" & id != "Kari" & id != "Johannes" & id_year != "Tiina_2011" &
                                          id_year != "Tiina_2012" & id_year != "Tiina_2013"& # remove these for unknown age adults
                                          id_year != "Mohammed_2019") # remove these for thermal uplift issues

first_migration <- modeling_data[which(modeling_data$Migration == 1),]
adult_migration <- modeling_data[which(modeling_data$Migration == "5 +"),]

clogit_results <- modeling_data %>% 
  group_by(Migration) %>% 
  nest() %>% 
  mutate(modelTCt = purrr::map(data, modelTCt),
         coefsTCt = purrr::map(modelTCt, coef),
         AIC_TCt = map_dbl(modelTCt, ~AIC(.)))

# extract and record coefficients
clogit_results <- unnest(clogit_results, coefsTCt) %>% 
  mutate(variable = c("Wind support","Crosswind","Thermal uplift")) %>% #,"Wind support","Crosswind","Thermal uplift","Wind support","Crosswind","Thermal uplift","Wind support","Crosswind","Thermal uplift","Wind support","Crosswind","Thermal uplift")) %>% 
  arrange(variable)

## perform data stream permutation
randomization_dist <- permute(first_migration, no.perm, case_, .id = "permutation")

# extract coefficients
models <- purrr::map(randomization_dist$perm, modelTCt)
glanced <- purrr::map_df(models, broom::glance, .id = "permutation")
coefficients <- purrr::map(models, coef)
coefficients <- purrr::map_df(coefficients, broom::tidy, .id = "permutation")

# distribution of null permutation statistics
hist(glanced$statistic.log, breaks = 100)
# confirm these are roughly uniform p-values
hist(glanced$p.value.log, breaks = 100)

## compare
# test against the unpermuted model to get a permutation p-value
sig_tail <-  sum(abs(clogit_results$coefsTCt[which(clogit_results$variable == "Wind support")]) < abs(coefficients$x[which(coefficients$names == "scaled_tail")]))/no.perm
sig_cross <-  sum(abs(clogit_results$coefsTCt[which(clogit_results$variable == "Crosswind")]) < abs(coefficients$x[which(coefficients$names == "scaled_cross")]))/no.perm
sig_thermal <-  sum(abs(clogit_results$coefsTCt[which(clogit_results$variable == "Thermal uplift")]) < abs(coefficients$x[which(coefficients$names == "scaled_thermal")]))/no.perm


# plot
no.bin <- 25
tail_coefs <- coefficients %>% filter(names=="scaled_tail")
first_tail <- ggplot(data=tail_coefs, aes(x)) + 
  geom_vline(aes(xintercept=clogit_results$coefsTCt[3]),
             color="red", size=1) +
  geom_histogram(fill="black", bins = no.bin) +
  labs(title = "Wind support", x="Coefficient value", y="Frequency") +
  theme_classic()

cross_coefs <- coefficients %>% filter(names=="scaled_cross")
first_cross <- ggplot(data=cross_coefs, aes(x)) + 
  geom_vline(aes(xintercept=clogit_results$coefsTCt[1]),
             color="red", size=1) +
  geom_histogram(fill="black", bins = no.bin) +
  labs(title = "Crosswind", x="Coefficient value", y="Frequency") +
  theme_classic()

thermal_coefs <- coefficients %>% filter(names=="scaled_thermal")
first_thermal <- ggplot(data=thermal_coefs, aes(x)) + 
  geom_vline(aes(xintercept=clogit_results$coefsTCt[2]),
             color="red", size=1) +
  geom_histogram(fill="black", bins = no.bin) +
  labs(title = "Thermal uplift", x="Coefficient value", y="Frequency") +
  theme_classic()

ggarrange(first_tail, first_cross, first_thermal,ncol=3, nrow=1, legend = "none")
