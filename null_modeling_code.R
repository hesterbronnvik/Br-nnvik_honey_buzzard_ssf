### Null modeling code for Honey Buzzard GPS data
### Hester Brønnvik, Summer 2021


## set working directory and load required packages
setwd("C:/Users/Tess Bronnvik/Desktop/Br-nnvik_honey_buzzard_ssf")
HBNM_packs <- c("lubridate", "modelr", "survival", "move", "ggpubr", "RNCEP", "mapview","tidyverse","viridis","hrbrthemes","patchwork", "groupdata2","cvms")
lapply(HBNM_packs, require, character.only = TRUE)

# import required functions
source("wind_support_Kami.R")

# define values
no.perm <- 100


### Prepare data
############################################################################################
## acquire data
half1 <- read.csv("half_track1_again_150s-2811382202164615735.csv",stringsAsFactors = F, header = T)
half2 <- read.csv("half_track2_again_150s-2052759029422935039.csv",stringsAsFactors = F, header = T)
annotated_data <- rbind(half1,half2) %>% 
  rename("u_wind" = "ECMWF.ERA5.PL.U.Wind",
         "v_wind" = "ECMWF.ERA5.PL.V.wind") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         year = format(as.POSIXct(timestamp,format="%Y-%m-%d %H:%M:%S"),"%Y"),
         id_year = paste(id,year,sep="_"),
         case_ = as.numeric(case_),
         tail = wind_support(u_wind, v_wind, heading),
         cross = cross_wind(u_wind, v_wind, heading))

colnames(annotated_data)[c(1,2,3,4,21,22)] <- c("lon1","lon2","lat1","lat2","thermal","sea_surface")

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
sum(thermalNA) # 28,094
annotated_data <- annotated_data[thermalNA != T,]

# select data with missing values
Edit_data <- annotated_data[annotated_data$id == "Edit",] # burst 1 step 31 (True) 


#remove those data from the data set
annotated_data <- annotated_data %>% filter(id != "Edit")# & id != "Annika")

# edit each ID with missing data to have the right number of steps remaining
Edit_data <- Edit_data[which(Edit_data$step_id_ != 31),] # remove step with 101 NAs

# add the IDs back to the data
annotated_data <- rbind(annotated_data,Edit_data)


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


# select 101 steps from each step_id_. Choosing the first 101 allows protection of True steps.
annotated_data <- annotated_data %>% 
  mutate(id_year = paste(id,year,sep = "_")) %>% 
  group_by(id_year, burst_, step_id_) %>% 
  slice(1:101) %>% 
  ungroup()


### normalize the data
annotated_data <- annotated_data %>%  
  mutate(scaled_cross = abs(as.numeric(scale(cross, center = T, scale = T))),
         scaled_tail = as.numeric(scale(tail, center = T, scale = T)),
         scaled_thermal = as.numeric(scale(thermal, center = T, scale = T)),
         #scaled_orographic = as.numeric(scale(orographic, center = T, scale = T))
         ) %>%
  ungroup()
save(annotated_data, file = "prepped_data.RData")
############################################################################################

### Model
############################################################################################
load("prepped_data.RData")

annotated_data$age <- NA
annotated_data$age[which(annotated_data$id == "Aarne" | 
                               annotated_data$id == "Annika" |
                               annotated_data$id == "Jari" |
                               annotated_data$id == "Johannes" |
                               annotated_data$id == "Jouko" |
                               annotated_data$id == "Mikko" |
                               annotated_data$id == "Paivi" |
                               annotated_data$id == "Tiina" |
                               annotated_data$id == "Kari")] <- "adult"
annotated_data$age[which(is.na(annotated_data$age))] <- "juvenile"


annotated_data$Migration <- NA # create an empty column to store values
annotated_data$Migration[which(annotated_data$id_year == "Jaana_2013" | # add year 2 info
                                 annotated_data$id_year == "Lars_2013" | 
                                 annotated_data$id_year == "Mohammed_2018" | 
                                 annotated_data$id_year == "Senta_2015" | 
                                 annotated_data$id_year == "Valentin_2015" | 
                                 annotated_data$id_year == "Aarne_2013" |
                                 annotated_data$id_year == "Annika_2012" | 
                                 annotated_data$id_year == "Jari_2013" | 
                                 annotated_data$id_year == "Johannes_2014" | 
                                 annotated_data$id_year == "Jouko_2012" |
                                 annotated_data$id_year == "Paivi_2014" |
                                 annotated_data$id_year == "Mikko_2012" |
                                 annotated_data$id_year == "Tiina_2012")] <- "2"

annotated_data$Migration[which(annotated_data$id_year == "Jaana_2014" | # add year 3 info
                                 annotated_data$id_year == "Lars_2014" | 
                                 annotated_data$id_year == "Senta_2016" |
                                 annotated_data$id_year == "Annika_2013" | 
                                 annotated_data$id_year == "Jari_2014" | 
                                 annotated_data$id_year == "Johannes_2015" |
                                 annotated_data$id_year == "Mikko_2013" | # missing
                                 annotated_data$id_year == "Paivi_2015" |
                                 annotated_data$id_year == "Mohammed_2019" |
                                 annotated_data$id_year == "Tiina_2013")] <- "3"

annotated_data$Migration[which(annotated_data$id_year == "Lars_2015" | # add year 4 info
                                 annotated_data$id_year == "Senta_2017" |
                                 annotated_data$id_year == "Annika_2014" |
                                 annotated_data$id_year == "Jouko_2014" |
                                 annotated_data$id_year == "Mikko_2014" |
                                 annotated_data$id_year == "Paivi_2016" |
                                 annotated_data$id_year == "Tiina_2014")] <- "4"

annotated_data$Migration[which(annotated_data$id_year == "Annika_2014" | # add adult info
                                 annotated_data$id_year == "Annika_2015" |
                                 annotated_data$id_year == "Jouko_2015" | 
                                 annotated_data$id_year == "Mikko_2015" | 
                                 annotated_data$id_year == "Paivi_2017" | 
                                 annotated_data$id_year == "Tiina_2015")] <- "5 +" 
annotated_data$Migration[which(is.na(annotated_data$Migration))] <- "1" # everything else is year 1



## fit a model to each group
# build the model function
modelTCt <- function(df) {
  clogit(case_ ~ scaled_tail + scaled_cross + scaled_thermal + strata(step_id_), data = df)
}

# exclude certain data
modeling_data <- annotated_data %>% filter(
    id_year != "Mikko_2011" &# only one step seen, non-convergence
    id_year != "Mikko_2012" &# only two steps seen, claims NA/NaN/Inf but there are none
    id_year != "Jouko_2013" &# only one step seen, non-convergence
    id_year != "Jouko_2015" )# 4 steps seen, 2 of them have NA ta_ in the True step

first_migration <- modeling_data[which(modeling_data$Migration == "1" & modeling_data$age == "juvenile"),]
adult_migration <- modeling_data[which(modeling_data$age == "adult" & modeling_data$Migration == "5 +"),]

clogit_results <- first_migration %>% 
  group_by(Migration) %>% 
  nest() %>% 
  mutate(model_col = purrr::map(data, modelTCt),
         coefsTCt = purrr::map(model_col, coef),
         AIC_TCt = map_dbl(model_col, ~AIC(.)),
         summary_stats = purrr::map(model_col, ~summary(.x)))
clogit_results %>%  pull(summary_stats)
# extract and record coefficients
clogit_results <- unnest(clogit_results, coefsTCt) %>% 
  mutate(variable = c("Wind support","Crosswind","Thermal uplift")) %>% #,"Wind support","Crosswind","Thermal uplift","Wind support","Crosswind","Thermal uplift","Wind support","Crosswind","Thermal uplift","Wind support","Crosswind","Thermal uplift")) %>% 
  arrange(variable)

## perform data stream permutation
random_model <- data.frame()
for (i in 1:no.perm) {
  rando <- mosaic::resample(first_migration, replace = F, groups = id_year_burst_step, shuffled = "case_", fixed = c("id","year","burst_","step_id_","scaled_cross","scaled_tail","scaled_thermal","timestamp", "ta_", "sl_"))
  temp_model <- modelTCt(rando)
  temp_coefs <- summary(temp_model)$coefficients
  random_model <- rbind(temp_coefs, random_model)
}

# extract and record coefficients
random_tails <- random_model[grep("tail", rownames(random_model)),]
random_crosses <- random_model[grep("cross", rownames(random_model)),]
random_thermals <- random_model[grep("thermal", rownames(random_model)),]


## compare
# test against the unpermuted model to get a permutation p-value
sig_tail <-  sum(abs(clogit_results$coefsTCt[which(clogit_results$variable == "Wind support")]) < abs(random_tails$coef))/no.perm
sig_cross <-  sum(abs(clogit_results$coefsTCt[which(clogit_results$variable == "Crosswind")]) < abs(random_crosses$coef))/no.perm
sig_thermal <-  sum(abs(clogit_results$coefsTCt[which(clogit_results$variable == "Thermal uplift")]) < abs(random_thermals$coef))/no.perm


# plot
no.bin <- 25
first_tail <- ggplot(data=random_tails, aes(coef)) + 
  geom_vline(aes(xintercept=clogit_results$coefsTCt[3]),
             color="red", size=1) +
  geom_histogram(fill="black", bins = no.bin) +
  labs(title = "Wind support (p = 0.9)", x="Coefficient value", y="Frequency") +
  theme_classic()

first_cross <- ggplot(data=random_crosses, aes(coef)) + 
  geom_vline(aes(xintercept=clogit_results$coefsTCt[1]),
             color="red", size=1) +
  geom_histogram(fill="black", bins = no.bin) +
  labs(title = "Crosswind (p = 0.94)", x="Coefficient value", y="Frequency") +
  theme_classic()

first_thermal <- ggplot(data=random_thermals, aes(coef)) + 
  geom_vline(aes(xintercept=clogit_results$coefsTCt[2]),
             color="red", size=1) +
  geom_histogram(fill="black", bins = no.bin) +
  labs(title = "Thermal uplift (p = 0.34) ", x="Coefficient value", y="Frequency") +
  theme_classic()

ggarrange(first_tail, first_cross, first_thermal,ncol=3, nrow=1, legend = "none")
############################################################################################