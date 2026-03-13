
# full protocol for processing Sierra BioCube and Acoustics data
# Laura Berman
# created 30 April 2025
# last updated 30 Sep 2025

# load packages ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(lubridate)
library(sf)
library(sp)
library(mapview)
library(terra)
library(raster)
library(ggspatial)
library(corrplot)
library(caret)
library(factoextra)
library(FactoMineR)
library(ggrepel)
library(rnaturalearth)
library(ggdark)
library(units)
library(ggpmisc)
library(basemaps)
library(reshape2)
library(cowplot)
library(geosphere)
library(scales)
library(patchwork)

# working directory  -----------------------------------------------------------
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Draft 1")


################################################################################
# starting with ARU bioacoustic data 
################################################################################

# load ARU locations -----------------------------------------------------------

ARUlocations <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/focal_aru_metadata.csv") # 2100 locations
ARUlocations <- ARUlocations %>% distinct(survey_year, cell_unit, box, .keep_all = TRUE) # 2073 unique locations

# load detections --------------------------------------------------------------

det21 <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/2021_filtered_detections.csv") # 41,264,999 obs
det22 <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/2022_filtered_detections.csv") # 40,424,580 obs
det23 <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/2023_filtered_detections.csv") # 27,057,557 obs

det <- rbind(det21, det22) # 81,689,579 obs (correct)
det <- rbind(det, det23) # 108,747,136 obs (correct)

# names(det) = "filename" "cell_unit" "year" "survey_date" "month" "day" "hour" "min" "second" "relative_sec" "species_code" "logit_score"
# still missing common names and exact times
rm(det21, det22, det23)


# load common names ------------------------------------------------------------

commonNames <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/species_codes_thresholds.csv") # 241 species


# load effort ------------------------------------------------------------------

effort21 <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/2021_filtered_effort.csv") # 134,322 obs
effort22 <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/2022_filtered_effort.csv") # 127,002 obs
effort23 <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/2023_filtered_effort.csv") # 110,349 obs

effort <- rbind(effort21, effort22) # 261,324 obs (correct)
effort <- rbind(effort, effort23) # 371,673 obs (correct)
# effort includes rows for days with 0 effort
# nrow(subset(effort, effort$aru_hours==7)) = 66,349 obs = total full recording days

#effort$year <- substr(effort$date, 1, 4)
unit_effort <- effort %>%
  subset(aru_hours==7) %>%
  group_by(cell_unit) %>%
  dplyr::summarize(startDate = min(date),
            endDate = max(date),
            active_days = n())
# unit_effort = 2018 obs = unique cell units x year
# min(unit_effort$active_days) = 1
# max(unit_effort$active_days) = 37

# after removing year = 895 obs
# min(unit_effort$active_days) = 10
# max(unit_effort$active_days) = 108

rm(effort21, effort22, effort23)

write_rds(unit_effort, "unit_effort.rds")

# summarize detections by site -------------------------------------------------

siteDetections <- det %>% 
  group_by(cell_unit, species_code) %>%
  dplyr::summarise(n_detections = n())
# siteDetections = 113,003 obs = one row per site x year x species combination, so this count doesn't mean much
# after removing year = 59,626 = one row per site x species

# combine detections, effort, and common names ---------------------------------

siteDetections <- merge(siteDetections, commonNames) # still only 113,003 obs. correct # 59626
# merging with unit_effort loses some observations because we only include sites with full 7 hour days
siteDetections <- merge(siteDetections, unit_effort) # 112,805 obs # 59626


# check max dist between redeployments -----------------------------------------

# Function to compute max pairwise distance between lat/long coords
max_dist_within_cell <- function(df) {
  coords <- df[, c("long", "lat")]  # Note: geosphere expects (lon, lat)
  if (nrow(coords) < 2) return(0)   # No distance possible if < 2 points
  dists <- distm(coords)           # Calculate all pairwise distances (in meters)
  return(max(dists))
}

# Apply this function by cell_unit
max_dists_by_cell <- ARUlocations %>%
  group_by(cell_unit) %>%
  summarise(max_distance_m = max_dist_within_cell(cur_data()))

hist(max_dists_by_cell$max_distance_m)
mean(max_dists_by_cell$max_distance_m)

# combine years with the same cell_unit ----------------------------------------

ARUlocations <- ARUlocations %>% 
  group_by(cell_unit) %>%
  dplyr::summarise(long = mean(long),
                   lat = mean(lat))
#names(ARUlocations)[names(ARUlocations) == 'survey_year'] <- 'year'
siteDetections <- merge(siteDetections, ARUlocations)

# remove rows where active days < 30 -------------------------------------------

siteDetections <- subset(siteDetections, siteDetections$active_days > 30) # 58671 # 59342

# pivot wider ------------------------------------------------------------------
# [I will want active_days back later]
siteDetections$obs_per_day <- siteDetections$n_detections / siteDetections$active_days
siteDetections <- siteDetections[c("cell_unit", "common_name", "obs_per_day", "long", "lat")]
siteDetections <- siteDetections %>%
  group_by(cell_unit) %>%
  pivot_wider(names_from = common_name,
              values_from = obs_per_day) # 2016 obs. one per site.
siteDetections[is.na(siteDetections)] <- 0



################################################################################
# remove burned areas
################################################################################


# get burn years ---------------------------------------------------------------

# load fire data 2024
# fire data comes from https://gis.data.ca.gov/datasets/CALFIRE-Forestry::california-fire-perimeters-all/explore 
fire <- read_sf(dsn = "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Fire perimeters/California_Fire_Perimeters_2024") # 22,261 fires
# 1898-04-01 til 2023-12-22
mapview(fire)
names(fire)[names(fire) == 'YEAR_'] <- 'burnYear'
fire <- fire[c("burnYear", "geometry")]
fire$burnYear <- as.numeric(fire$burnYear)

# only keep fires 2018-2023
recentFires <- subset(fire, fire$burnYear >= 2018)
mapview(ARUlocations.sf) + mapview(recentFires)

# get all the years when ARU locations burned
sf_use_s2(FALSE)
ARUlocations.sf <- st_as_sf(ARUlocations, coords = c("long", "lat"), crs=4326, remove=FALSE)
fire <- st_transform(fire, crs = crs(ARUlocations.sf))
firePoints <- st_join(ARUlocations.sf, fire, multiple="all") # 2645 obs/burns (thats more than one burn per site)
mapview(firePoints, zcol="burnYear")
sf_use_s2(TRUE)

# keep only the most recent burn for each point
burnYear <- firePoints %>%
  group_by(cell_unit) %>% 
  slice_max(burnYear) %>% 
  ungroup() # 902
burnYear <- burnYear %>% distinct(cell_unit, .keep_all = TRUE) # 901 obs

rm(fire, firePoints)


# add burn year to siteDetections ----------------------------------------------

siteDetections <- merge(burnYear, siteDetections) # 889 obs, 106 vars

# remove points that burned ----------------------------------------------------

siteDetections <- subset(siteDetections, is.na(siteDetections$burnYear) | siteDetections$burnYear < 2018) # 649 obs

# get tahoe and yosemite boxes -------------------------------------------------

tahoe <- read_sf("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/tahoe box/tahoe_shp.shp")
yosemite <- read_sf("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/yosemite box/polygon_yosemite_box.shp")
mapview(tahoe) +
  mapview(yosemite)

# exclude ARUs not in boxes ----------------------------------------------------
mapview(siteDetections) + mapview(tahoe) + mapview(yosemite)
# not necessary

# save siteDetections ----------------------------------------------------------

write_csv(siteDetections, "siteDetections_20250501.csv")




################################################################################
# Get foliar traits
################################################################################


# buffer site locations
sites.buffer <- st_buffer(siteDetections, dist = 150) # 649
# a 150 m radius buffer = about 9 hectares, avg home range size of small passerines
sites.buffer <- sites.buffer[c("cell_unit", "long", "lat")]

# Nitrogen
NitrogenYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Nitrogen_30m_polished.tif")
NitrogenTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Nitrogen_30m_polished.tif")
x <- list(NitrogenYosemite, NitrogenTahoe)
Nitrogen <- do.call(merge, x)
mapview(Nitrogen) +
  mapview(sites.buffer)
rm(NitrogenTahoe, NitrogenYosemite)

x <- raster::extract(Nitrogen$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Nitrogen"
x$Nitrogen[x$Nitrogen == 0] <- NA
foliarTraits <- x # 649 x 5
mapview(Nitrogen) +
  mapview(foliarTraits, zcol="Nitrogen")


# LMA
LMAYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_LMA_30m_polished.tif")
LMATahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_LMA_30m_polished.tif")
x <- list(LMAYosemite, LMATahoe)
LMA <- do.call(merge, x)
mapview(LMA) +
  mapview(sites.buffer)
rm(LMATahoe, LMAYosemite)

x <- raster::extract(LMA$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "LMA"
x$LMA[x$LMA == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 x 6
mapview(LMA) +
  mapview(foliarTraits, zcol="LMA")


# Cellulose
CelluloseYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Cellulose_30m_polished.tif")
CelluloseTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Cellulose_30m_polished.tif")
x <- list(CelluloseYosemite, CelluloseTahoe)
Cellulose <- do.call(merge, x)
rm(CelluloseTahoe, CelluloseYosemite)

x <- raster::extract(Cellulose$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Cellulose"
x$Cellulose[x$Cellulose == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 


# Calcium
CalciumYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Calcium_30m_polished.tif")
CalciumTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Calcium_30m_polished.tif")
x <- list(CalciumYosemite, CalciumTahoe)
Calcium <- do.call(merge, x)
rm(CalciumTahoe, CalciumYosemite)

x <- raster::extract(Calcium$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Calcium"
x$Calcium[x$Calcium == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 

# Chlorophylls
ChlorophyllsYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Chlorophylls_30m_polished.tif")
ChlorophyllsTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Chlorophylls_30m_polished.tif")
x <- list(ChlorophyllsYosemite, ChlorophyllsTahoe)
Chlorophylls <- do.call(merge, x)
rm(ChlorophyllsTahoe, ChlorophyllsYosemite)

x <- raster::extract(Chlorophylls$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Chlorophylls"
x$Chlorophylls[x$Chlorophylls == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 


# Fiber
FiberYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Fiber_30m_polished.tif")
FiberTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Fiber_30m_polished.tif")
x <- list(FiberYosemite, FiberTahoe)
Fiber <- do.call(merge, x)
rm(FiberTahoe, FiberYosemite)

x <- raster::extract(Fiber$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Fiber"
x$Fiber[x$Fiber == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 


# Lignin
LigninYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Lignin_30m_polished.tif")
LigninTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Lignin_30m_polished.tif")
x <- list(LigninYosemite, LigninTahoe)
Lignin <- do.call(merge, x)
rm(LigninTahoe, LigninYosemite)

x <- raster::extract(Lignin$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Lignin"
x$Lignin[x$Lignin == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 


# NSC
NSCYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_NSC_30m_polished.tif")
NSCTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_NSC_30m_polished.tif")
x <- list(NSCYosemite, NSCTahoe)
NSC <- do.call(merge, x)
rm(NSCTahoe, NSCYosemite)

x <- raster::extract(NSC$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "NSC"
x$NSC[x$NSC == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649


# Phenolics
PhenolicsYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Phenolics_30m_polished.tif")
PhenolicsTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Phenolics_30m_polished.tif")
x <- list(PhenolicsYosemite, PhenolicsTahoe)
Phenolics <- do.call(merge, x)
rm(PhenolicsTahoe, PhenolicsYosemite)

x <- raster::extract(Phenolics$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Phenolics"
x$Phenolics[x$Phenolics == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649


# Phosphorus
PhosphorusYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Phosphorus_30m_polished.tif")
PhosphorusTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Phosphorus_30m_polished.tif")
x <- list(PhosphorusYosemite, PhosphorusTahoe)
Phosphorus <- do.call(merge, x)
rm(PhosphorusTahoe, PhosphorusYosemite)

x <- raster::extract(Phosphorus$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Phosphorus"
x$Phosphorus[x$Phosphorus == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649


# Potassium
PotassiumYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Potassium_30m_polished.tif")
PotassiumTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Potassium_30m_polished.tif")
x <- list(PotassiumYosemite, PotassiumTahoe)
Potassium <- do.call(merge, x)
rm(PotassiumTahoe, PotassiumYosemite)

x <- raster::extract(Potassium$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Potassium"
x$Potassium[x$Potassium == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649


# Starch
StarchYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Starch_30m_polished.tif")
StarchTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Starch_30m_polished.tif")
x <- list(StarchYosemite, StarchTahoe)
Starch <- do.call(merge, x)
rm(StarchTahoe, StarchYosemite)

x <- raster::extract(Starch$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Starch"
x$Starch[x$Starch == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649 x 5


# Sugar
SugarYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Sugar_30m_polished.tif")
SugarTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Sugar_30m_polished.tif")
x <- list(SugarYosemite, SugarTahoe)
Sugar <- do.call(merge, x)
rm(SugarTahoe, SugarYosemite)

x <- raster::extract(Sugar$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Sugar"
x$Sugar[x$Sugar == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649


# Sulfur
SulfurYosemite <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Sulfur_30m_polished.tif")
SulfurTahoe <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Sulfur_30m_polished.tif")
x <- list(SulfurYosemite, SulfurTahoe)
Sulfur <- do.call(merge, x)
rm(SulfurTahoe, SulfurYosemite)

x <- raster::extract(Sulfur$layer, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
x <- cbind(sites.buffer, x)
colnames(x)[colnames(x)=="x"] <- "Sulfur"
x$Sulfur[x$Sulfur == 0] <- NA
x <- as.data.frame(x)
foliarTraits <- merge(foliarTraits, x) # 649

# save foliar traits -----------------------------------------------------------

write_rds(foliarTraits, "foliarTraits_20250502.rds")
x <- as.data.frame(foliarTraits)
x$geometry <- NULL 
write_csv(x, "foliarTraits_20250502.csv")

# merge foliar traits with siteDetections --------------------------------------

foliarTraits <- as.data.frame(foliarTraits)
foliarTraits$geometry <- NULL 
siteDetections_foliarTraits <- merge(siteDetections, foliarTraits)

# subset to remove NA rows -----------------------------------------------------

siteDetections_foliarTraits <- siteDetections_foliarTraits[!is.na(siteDetections_foliarTraits$Nitrogen), ] #580

# save foliar traits x siteDetections ------------------------------------------

write_rds(siteDetections_foliarTraits, "siteDetections_foliarTraits_20250502.rds")
x <- as.data.frame(siteDetections_foliarTraits)
x$geometry <- NULL
write_csv(x, "siteDetections_foliarTraits_20250502.csv")


################################################################################
# Get biocube variables
################################################################################

# update buffer -----------------------------------------------------------------

# buffer site locations
sites.buffer <- st_buffer(siteDetections_foliarTraits, dist = 150) # 580
# a 150 m radius buffer = about 9 hectares, avg home range size of small passerines
sites.buffer <- sites.buffer[c("cell_unit", "long", "lat")]


# get biocube paths and files --------------------------------------------------

CAbiocube_files <- list.files("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/California biocube layers")
path_CA <- "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/California biocube layers/"


# get values for first var -----------------------------------------------------

x <- raster(paste(path_CA, CAbiocube_files[1], sep="")) # rasterize the file
x <- raster::extract(x, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean) # extract values in buffer around each ARU
x <- as.data.frame(cbind(sites.buffer, x)) # make it a dataframe
colnames(x)[colnames(x)=="x"] <- paste(CAbiocube_files[1]) # give it an informative name
x$geometry <- NULL 
BioCube_vars <- x # this is the one line that changes outside the loop


# loop the rest ----------------------------------------------------------(slow)

for(i in 2:length(CAbiocube_files)) {
  x <- raster(paste(path_CA, CAbiocube_files[i], sep=""))
  x <- raster::extract(x, sites.buffer, weights=TRUE, na.rm=TRUE, fun=mean)
  x <- as.data.frame(cbind(sites.buffer, x))
  colnames(x)[colnames(x)=="x"] <- paste(CAbiocube_files[i])
  x$geometry <- NULL 
  BioCube_vars <- merge(BioCube_vars, x)
}

# save it
write_rds(BioCube_vars, "BioCube_vars_20250520.rds")

# tidy variable names ----------------------------------------------------------

# load tidy names
tidy <- read_csv("tidyNames.csv")

# Create a named vector for renaming
rename_map <- setNames(tidy$VariableName, tidy$FileName)

# Rename only columns in BioCube that are in tidy$file
names(BioCube_vars)[names(BioCube_vars) %in% names(rename_map)] <- rename_map[names(BioCube_vars)[names(BioCube_vars) %in% names(rename_map)]]


# merge with detections and foliar traits --------------------------------------

siteDetections_foliarTraits_BioCube <- merge(siteDetections_foliarTraits, BioCube_vars)


# save -------------------------------------------------------------------------

write_rds(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.rds")
write_csv(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.csv")


################################################################################
# add USGS CONUS GAP 2011 landcover type
################################################################################

# load data
YosemiteGAP <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Yosemite_GAP_2011.tif")
TahoeGAP <- raster("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Tahoe_GAP_2011.tif")
YosemiteGAP$landcover[YosemiteGAP$landcover==0] <- NA
TahoeGAP$landcover[TahoeGAP$landcover==0] <- NA
x <- list(YosemiteGAP, TahoeGAP)
GAP <- do.call(merge, x)
mapview(GAP)

# site locations
sites <- siteDetections_foliarTraits_BioCube[c("cell_unit", "long", "lat")]

# extract
landcover_ID <- raster::extract(GAP$layer, sites, na.rm=TRUE)
siteDetections_foliarTraits_BioCube <- as.data.frame(cbind(siteDetections_foliarTraits_BioCube, landcover_ID))

# add key names
coverNames <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/USGS_CONUS_GAP_2011_key.csv")
coverNames <- coverNames[c("landcover_ID", "landcover_type")]
siteDetections_foliarTraits_BioCube <- merge(siteDetections_foliarTraits_BioCube, coverNames)
siteDetections_foliarTraits_BioCube$landcover_ID <- NULL


################################################################################
#  remove species with insufficient data
################################################################################

# at how many sites is each species detected? ----------------------------------
x <- siteDetections_foliarTraits_BioCube[,c(5:96)]
x <- sort(colSums(x != 0))
speciesSites <- data.frame(species = names(x), totalSites = as.numeric(x))

# remove columns for species <= 30
siteDetections_foliarTraits_BioCube$Tree.Swallow <- NULL
siteDetections_foliarTraits_BioCube$California.Thrasher <- NULL
siteDetections_foliarTraits_BioCube$Song.Sparrow <- NULL
siteDetections_foliarTraits_BioCube$Wild.Turkey <- NULL
siteDetections_foliarTraits_BioCube$Brown.headed.Cowbird <- NULL
siteDetections_foliarTraits_BioCube$Vaux.s.Swift <- NULL
siteDetections_foliarTraits_BioCube$Yellow.breasted.Chat <- NULL


# how many times in total was a species detected? ------------------------------

# add back active days
x <- merge(siteDetections_foliarTraits_BioCube, unit_effort[c("cell_unit", "active_days")])
x <- sort(colSums(x[ , 5:96] * x$active_days, na.rm = TRUE))
speciesDetections <- data.frame(species = names(x), totalDetections = as.numeric(x))

# remove species with < 300 detections
siteDetections_foliarTraits_BioCube$Brewer.s.Sparrow <- NULL
siteDetections_foliarTraits_BioCube$Sooty.Grouse <- NULL

# 92 species

# save -------------------------------------------------------------------------

write_rds(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.rds")
write_csv(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.csv")

# 580 obs, 207 var

################################################################################
#  which points are too close to each other?
################################################################################

# Compute full distance matrix
dists <- st_distance(sites)
# Replace self-distances with NA (or Inf)
diag(dists) <- NA
# Find the minimum distance to another point for each point
min_dists <- apply(dists, 1, min, na.rm = TRUE)
sites$min_dist <- min_dists

# add back effort
x <- merge(sites, unit_effort)

# Too close = less than 500 m apart

# pick which too-close sites are better to exclude
# C5196_U1 & C5196_U2. remove U1.
# C3824_U1 & C3824_U2. remove U2.
# C6940_U1 & C6940_U2. remove U2.

# remove them
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCube[siteDetections_foliarTraits_BioCube$cell_unit != "C5196_U1", ]
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCube[siteDetections_foliarTraits_BioCube$cell_unit != "C3824_U2", ]
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCube[siteDetections_foliarTraits_BioCube$cell_unit != "C6940_U2", ]

# 577 obs, 207 var

# save updated version

write_rds(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.rds")
write_csv(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.csv")

# how many unique cell units? --------------------------------------------------
unique(str_extract(siteDetections_foliarTraits_BioCube$cell_unit, ".*(?=_)"))
length(unique(str_extract(siteDetections_foliarTraits_BioCube$cell_unit, ".*(?=_)"))) # 285

# average distance between ARUs in the same cell? ------------------------------
sites <- siteDetections_foliarTraits_BioCube[c("cell_unit", "long", "lat")]
sites$cell <- str_extract(sites$cell_unit, ".*(?=_)")
sites_sf <- st_as_sf(sites, coords = c("long", "lat"), crs = 4326) # Use WGS84
sites_sf <- st_transform(sites_sf, 3857)
avg_distances <- sites_sf %>%
  group_by(cell) %>%
  summarise(count= n(),
            avg_dist = {
    geom <- geometry
    if (n() < 2) {
      NA_real_  # Not enough points to calculate distance
    } else {
      dists <- st_distance(geom)
      dists <- as.numeric(dists)
      # Remove self-distances and duplicates
      dists <- dists[lower.tri(matrix(NA, n(), n()))]
      mean(dists, na.rm = TRUE)
    }
  })
mean(avg_distances$avg_dist, na.rm=TRUE) #1206 meters

################################################################################
#  make a nice map
################################################################################

# update sites
sites <- siteDetections_foliarTraits_BioCube[c("cell_unit", "long", "lat")]
sites <- st_as_sf(sites, coords = c("long", "lat"), crs=4326, remove=FALSE)

# Create bounding box
ext <- st_as_sfc(st_bbox(c(
  xmin = -121, xmax = -118.5,
  ymin = 36.5, ymax = 39.5
), crs = 4326))

# Transform to EPSG:3857 (Web Mercator)
ext_3857 <- st_transform(ext, 3857)
sites_3857 <- st_transform(sites, 3857)

# set default basemap settings
basemaps::set_defaults(ext = ext_3857)

#plot it
ggplot() +
  basemap_gglayer(ext_3857, map_service = "esri", map_type = "world_dark_gray_base") +
  scale_fill_identity() + 
  geom_sf(data = sites_3857, color = "navy", size = 1, show.legend = FALSE, alpha=0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = st_bbox(ext_3857)[c("xmin", "xmax")],
           ylim = st_bbox(ext_3857)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")


ggplot() +
  basemap_gglayer(ext_3857, map_service = "esri", map_type = "world_physical_map") +
  scale_fill_identity() + 
  geom_sf(data = sites_3857, color = "navy", size = 1, show.legend = FALSE, alpha=0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = st_bbox(ext_3857)[c("xmin", "xmax")],
           ylim = st_bbox(ext_3857)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")

ggsave("point_map.pdf", height=7, width=5)

################################################################################
# check for colinearity
################################################################################

# 90% cutoff

cor_matrix <- cor(siteDetections_foliarTraits_BioCube[,97:182], use = "complete.obs")
corrplot(cor_matrix,
         method = "color",
         order = "hclust",         
         hclust.method = "complete", 
         addrect = 3,             
         tl.cex = 0.8)

high_corr_matrix <- cor_matrix
high_corr_matrix[abs(high_corr_matrix) < 0.95 | diag(ncol(high_corr_matrix)) == 1] <- 0  
corrplot(high_corr_matrix,
         method = "color",
         order = "hclust",
         hclust.method = "complete", 
         tl.cex = 0.8)

high_corr_pairs <- as.data.frame(as.table(cor_matrix)) %>%
  filter(Var1 != Var2, abs(Freq) > 0.9)

# in Cali, there are dry summers and wet winters, so some bioclim vars are duplicate.
# bio16 and bio19 are 0.999 identical. precip of wettest/coldest quarter. remove bio16(wettest)?
# bio09 and bio10 are 0.999 identical. temp of dryest/warmest quarter. remove bio09(dryest)?
# ggd5 and gdd10 are 0.999 identical. remove gdd5?
# bio08 and bio11 are 0.998 identical. temp of wettest/coldest quarter. remove bio11(coldest)?
# swe and ngd0 are 0.98 identical. remove ngd0
# swe and scd 0.98 identical

# All these are > 95% correlated with DTM
# bio01, bio05, bio06, bio08, bio10, gdd10, ngd10, ngd5, annualPET, PETDryQuart, PETseas, PETWarmQuart, thermInd

# remove
siteDetections_foliarTraits_BioCube$bio16 <- NULL
siteDetections_foliarTraits_BioCube$bio09 <- NULL
siteDetections_foliarTraits_BioCube$gdd5 <- NULL
siteDetections_foliarTraits_BioCube$bio11 <- NULL
siteDetections_foliarTraits_BioCube$bio01 <- NULL
siteDetections_foliarTraits_BioCube$bio05 <- NULL
siteDetections_foliarTraits_BioCube$bio06 <- NULL
siteDetections_foliarTraits_BioCube$bio08 <- NULL
siteDetections_foliarTraits_BioCube$bio10 <- NULL
siteDetections_foliarTraits_BioCube$gdd0 <- NULL
siteDetections_foliarTraits_BioCube$gdd10 <- NULL
siteDetections_foliarTraits_BioCube$ngd10 <- NULL
siteDetections_foliarTraits_BioCube$ngd5 <- NULL
siteDetections_foliarTraits_BioCube$annualPET <- NULL
siteDetections_foliarTraits_BioCube$PETDryQuart <- NULL
siteDetections_foliarTraits_BioCube$PETseas <- NULL
siteDetections_foliarTraits_BioCube$PETWarmQuart <- NULL
siteDetections_foliarTraits_BioCube$ngd0 <- NULL
siteDetections_foliarTraits_BioCube$thermInd <- NULL
siteDetections_foliarTraits_BioCube$TRI <- NULL
siteDetections_foliarTraits_BioCube$scd <- NULL
siteDetections_foliarTraits_BioCube$bio19 <- NULL
siteDetections_foliarTraits_BioCube$Starch <- NULL

# save updated version
write_rds(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.rds")
write_csv(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.csv")


################################################################################
# Niche Space PCA
################################################################################

# convert to count data --------------------------------------------------------
names(siteDetections_foliarTraits_BioCube)
roundedCounts <- siteDetections_foliarTraits_BioCube %>% mutate(across(5:96, ~ round(.)))  # double check these are species columns. maybe multiplying by 7 (1 detection per week) to keep less vocal species # round(. * 7)
roundedCounts$geometry <- NULL

# name vars --------------------------------------------------------------------
species <- colnames(roundedCounts)[5:96]
spatVars <- colnames(roundedCounts)[97:182]

# calculate PCA ----------------------------------------------------------------
traitsPCA <-  prcomp(roundedCounts[,spatVars], center=TRUE, scale = TRUE)
# test <- PCA(roundedCounts[,spatVars])
# add PCs to rounded dataframe
roundedPCA <- cbind(roundedCounts, traitsPCA$x)


################################################################################
# PCA loadings (supp figures)
################################################################################


# what are the PCA loadings? ---------------------------------------------------
loadings <- traitsPCA$rotation[, 1:4]  # PC1 loadings
loadings <- as.data.frame(loadings)
loadings$variable <- rownames(loadings)
rownames(loadings) <- NULL

# get biocube category
tidy <- read_csv("tidyNames.csv")
tidy$fileName <- NULL
loadings <- merge(loadings, tidy, by.x= "variable", by.y= "VariableName")
loadings$FileName <- NULL

# plot loadings - supplemental figure 1
ggplot(loadings, aes(x=reorder(variable, abs(PC1)), y=abs(PC1), fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + coord_flip() +
  xlab("Habitat Variable") +
  ylab("PC1 loading (absolute value)") +
  theme_minimal()
ggsave("Loading_PC1_20250624.pdf", height=10, width=7)

# supplemental figure 2
ggplot(loadings, aes(x=reorder(variable, abs(PC2)), y=abs(PC2), fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + coord_flip() +
  xlab("Habitat Variable") +
  ylab("PC2 loading (absolute value)") +
  theme_minimal()
ggsave("Loading_PC2_20250624.pdf", height=10, width=7)

# supplemental figure 3
ggplot(loadings, aes(x=reorder(variable, abs(PC3)), y=abs(PC3), fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + coord_flip() +
  xlab("Habitat Variable") +
  ylab("PC3 loading (absolute value)") +
  theme_minimal()
ggsave("Loading_PC3_20260313.pdf", height=10, width=7)

# supplemental figure 4
ggplot(loadings, aes(x=reorder(variable, abs(PC4)), y=abs(PC4), fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + coord_flip() +
  xlab("Habitat Variable") +
  ylab("PC4 loading (absolute value)") +
  theme_minimal()
ggsave("Loading_PC4_20260313.pdf", height=10, width=7)

# visualize eigenvectors - supplemental figure 5
fviz_pca_var(traitsPCA, axes = c(1, 2))
ggsave("Eigenvectors_PC12.pdf", width=10, height=10)

# supplemental figure 6
fviz_pca_var(traitsPCA, axes = c(3, 4))
ggsave("Eigenvectors_PC34.pdf", width=10, height=10)

fviz_contrib(traitsPCA, choice="var", axes = 1, sort.val ="asc") + 
  coord_flip() + labs(title = "Contribution of variables to PC1 (%)") +
  ylab("Habitat")

fviz_contrib(traitsPCA, choice="var", axes = 2, sort.val ="asc") + 
  coord_flip() + labs(title = "Contribution of variables to PC2 (%)") +
  ylab("Habitat")

################################################################################
# render landcover labels
################################################################################

# the plot
landcover_pca_plot <- fviz_pca_ind(traitsPCA,
                                   label = "none",
                                   habillage = siteDetections_foliarTraits_BioCube$landcover_type) +
  guides(color = guide_legend(ncol = 1),
         fill = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1)) +
  geom_text(aes(x = -7, y = -8, label = "Oak"), inherit.aes = FALSE) +
  geom_text(aes(x =  6, y = -2, label = "Pine"), inherit.aes = FALSE) +
  geom_text(aes(x =  0, y = 1, label = "Mixed Conifer"), inherit.aes = FALSE) +
  ggtitle(NULL)
landcover_pca_plot

# keep it for later as a supp figure
ggsave("landcover_pca_20250624.PDF", plot=landcover_pca_plot, height=7, width=12)

# plot pc3 and pc4 - supplemental figure 7
fviz_pca_ind(traitsPCA,
             axes = c(2, 3),
             label = "none",
             habillage = siteDetections_foliarTraits_BioCube$landcover_type) +
  guides(color = guide_legend(ncol = 1),
         fill = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1)) +
  ggtitle(NULL)

make_pca_plot <- function(ax) {
  fviz_pca_ind(traitsPCA,
               axes = ax,
               label = "none",
               habillage = siteDetections_foliarTraits_BioCube$landcover_type) +
    ggtitle(paste0("PC", ax[1], " vs PC", ax[2])) +
    theme(legend.position = "none")
}

p1 <- make_pca_plot(c(4,1))
p2 <- make_pca_plot(c(3,1))
p3 <- make_pca_plot(c(4,2))
p4 <- make_pca_plot(c(3,2))
p5 <- make_pca_plot(c(4,3))

(p1 | p2) /
  (p3 | p4) /
  (p5 | patchwork::plot_spacer())
ggsave("landcover_pca_PC34panels_20260313.pdf", height=10, width=10)

# extract data
landcover_pca_scores <- as.data.frame(landcover_pca_plot$data)

centroids <- aggregate(cbind(x, y) ~ Groups, data = landcover_pca_scores, FUN = mean)

group_sizes <- landcover_pca_scores %>%
  count(Groups) %>%
  rename(size = n)

landcover_centroids <- left_join(centroids, group_sizes, by = "Groups")

# pick which landcover labels to keep ------------------------------------------
# these are all types which occur 10+ times in the data
landcover_centroids <- landcover_centroids %>%
  filter(Groups %in% c("Mediterranean California Mesic Mixed Conifer Forest and Woodland",
                       "Mediterranean California Dry-Mesic Mixed Conifer Forest and Woodland",
                       #"California Lower Montane Blue Oak-Foothill Pine Woodland and Savanna",
                       "California Central Valley Mixed Oak Savanna",
                       "Mediterranean California Red Fir Forest",
                       "California Montane Jeffrey Pine-(Ponderosa Pine) Woodland",
                       "Harvested Forest - Northwestern Conifer Regeneration",
                       "Sierra Nevada Subalpine Lodgepole Pine Forest and Woodland"))

# rename with shorter names
landcover_centroids$Groups <- as.character(landcover_centroids$Groups)
landcover_centroids$Groups[landcover_centroids$Groups == "Mediterranean California Mesic Mixed Conifer Forest and Woodland"] <- "Mesic Mixed Conifer Forest and Woodland"
landcover_centroids$Groups[landcover_centroids$Groups == "Mediterranean California Dry-Mesic Mixed Conifer Forest and Woodland"] <- "Dry-Mesic Mixed Conifer Forest and Woodland"
landcover_centroids$Groups[landcover_centroids$Groups == "California Lower Montane Blue Oak-Foothill Pine Woodland and Savanna"] <- "Lower Montane Blue Oak-Foothill Pine Woodland and Savanna"
landcover_centroids$Groups[landcover_centroids$Groups == "California Central Valley Mixed Oak Savanna"] <- "Mixed Oak Savanna"
landcover_centroids$Groups[landcover_centroids$Groups == "Mediterranean California Red Fir Forest"] <- "Red Fir Forest"
landcover_centroids$Groups[landcover_centroids$Groups == "California Montane Jeffrey Pine-(Ponderosa Pine) Woodland"] <- "Montane Jeffrey Pine-(Ponderosa Pine) Woodland"
landcover_centroids$Groups[landcover_centroids$Groups == "Sierra Nevada Subalpine Lodgepole Pine Forest and Woodland"] <- "Subalpine Lodgepole Pine Forest and Woodland"

################################################################################
# test plot for errored species polygons
################################################################################

# some species weren't plotting correctly because not enough points
# Solution #1) 1 = 1 detection per 10 days (x10 increase) (NOPE)
# Solution #2) add jitter to all points (non-identical coordinates)

species <- names(roundedPCA)[5:96]
species
species1 <- "Swainson.s.Thrush"
species1 <- "American.Robin"
species1 <- "Downy.Woodpecker" 
species1 <- "Lawrence.s.Goldfinch" 
species1 <- "Eurasian.Collared.Dove" 


get_density_threshold <- function(kde, prob = 0.90) {
  z_vals <- sort(as.vector(kde$z), decreasing = TRUE)
  cum_probs <- cumsum(z_vals) / sum(z_vals)
  threshold_index <- which(cum_probs >= prob)[1]
  z_thresh <- z_vals[threshold_index]
  return(z_thresh)
}

uncounted_data1 <- tidyr::uncount(roundedPCA, get(species1))
uncounted_data1$PC1 <- uncounted_data1$PC1 + rnorm(nrow(uncounted_data1), mean = 0, sd = 0.5)
uncounted_data1$PC2 <- uncounted_data1$PC2 + rnorm(nrow(uncounted_data1), mean = 0, sd = 0.5)
ggplot(uncounted_data1, aes(x=PC1, y=PC2)) +  geom_density2d()
d <- MASS::kde2d(x = uncounted_data1$PC1, y = uncounted_data1$PC2, lims = c(15, -20, -15, 15), n = 200)
dens <- data.frame(expand.grid(x = d$x, y = d$y), z = as.vector(d$z))
z_thresh <- get_density_threshold(d, prob = 0.90)
ggplot() +
  geom_contour_filled(data = dens, aes(x = x, y = y, z = z),
                      breaks = c(z_thresh, max(dens$z)))

ggplot(roundedPCA, aes(x=PC1, y=PC2, color=get(species1))) + geom_point() + scale_color_viridis_c()
ggplot(uncounted_data1, aes(x=PC1, y=PC2)) + geom_point(alpha=0.05)

################################################################################
# Species Niche Polygons
################################################################################


# make a function to set a specific % of points inside each contour ------------
get_density_threshold <- function(kde, prob = 0.90) {
  z_vals <- sort(as.vector(kde$z), decreasing = TRUE)
  cum_probs <- cumsum(z_vals) / sum(z_vals)
  threshold_index <- which(cum_probs >= prob)[1]
  z_thresh <- z_vals[threshold_index]
  return(z_thresh)
}


# Make function: calculate niche polygons --------------------------------------
calculate_polygons <- function(species1) {
  # Species count data
  uncounted_data1 <- tidyr::uncount(roundedPCA, get(species1))
  
  # add jitter
  uncounted_data1$PC1 <- uncounted_data1$PC1 + rnorm(nrow(uncounted_data1), mean = 0, sd = 0.5)
  uncounted_data1$PC2 <- uncounted_data1$PC2 + rnorm(nrow(uncounted_data1), mean = 0, sd = 0.5)

  # Try KDE and plot extraction
  p1dat <- tryCatch({
    # 2d kernel density estimation
    d <- MASS::kde2d(uncounted_data1$PC1, uncounted_data1$PC2, lims=c(15,-20,-15,15), n=200)    # lims() here prevents polygon edges from getting cut off. n sets resolution. 
    dens <- data.frame(expand.grid(x = d$x, y = d$y), z = as.vector(d$z))
    #plot it
    z_thresh <- get_density_threshold(d, prob = 0.90) # set 90% threshold
    p1 <- ggplot() +
      geom_contour_filled(data = dens, aes(x = x, y = y, z = z),
                          breaks = c(z_thresh, max(dens$z)))
    #get build
    ggplot_build(p1)$data[[1]]
  }, error = function(e) {
    warning(paste("Skipping species", species1, "- KDE failed:", e$message))
    return(NULL)
  })
  
  # If KDE failed
  if (is.null(p1dat)) return(NULL)
  
  # Close each polygon by repeating the first coordinate
  p1_closed <- p1dat %>%
    group_by(subgroup, level) %>%
    reframe(
      x = c(x, x[1]),
      y = c(y, y[1])
    )
  
  # Split into list by subgroup and level
  p1_split <- p1_closed %>%
    group_by(subgroup, level) %>%
    group_split()
  
  # Create sf POLYGON objects
  polygons <- lapply(p1_split, function(df) {
    coords <- cbind(df$x, df$y)
    st_polygon(list(coords))
  })
  
  # Extract levels for each polygon
  levels_vec <- sapply(p1_split, function(df) unique(df$level))
  
  # Combine into sf object
  p1_sf <- st_sf(
    species = rep(species1, length(levels_vec)),
    level = levels_vec,
    geometry = st_sfc(polygons)
  )
  
  # Union polygons by level and cast as MULTIPOLYGON
  p1_multi <- p1_sf %>%
    group_by(species, level) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    st_cast("MULTIPOLYGON")
  
  # Return the sf object
  return(p1_multi)
}


# Make a dataframe with all species polygons -----------------------------------

speciesPolygons <- map_dfr(species, ~ calculate_polygons(.x), .id = "species1")
speciesPolygons$species1 <- NULL
speciesPolygons$level <-NULL

# save it ----------------------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Draft 1")
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Shiny/Avispace4.0")

write_rds(speciesPolygons, "speciesPolygons_jittered_detPerDay_20250624.rds")


################################################################################
# plot 2 species niche overlap
################################################################################

# pick my favorite eigenvectors ------------------------------------------------
selected_vars <- c("Nitrogen", "LMA", "Phenolics", "Cellulose", "Chlorophylls", "Lignin")

# Get the eigenvalues (variances)
eig_values <- traitsPCA$sdev^2

# Scale var_data by sqrt of eigenvalues
var_data <- as.data.frame(traitsPCA$rotation)
var_data$PC1 <- var_data$PC1 * sqrt(eig_values[1])
var_data$PC2 <- var_data$PC2 * sqrt(eig_values[2])
var_data$label <- rownames(var_data)

# pick your species ------------------------------------------------------------
species <- colnames(roundedCounts)[5:96]
sort(species)
Species1 <- "Acorn.Woodpecker"
Species2 <- "California.Scrub.Jay"

# get polygons -----------------------------------------------------------------
species1_poly <- speciesPolygons %>% filter(species == Species1)
species2_poly <- speciesPolygons %>% filter(species == Species2)

# Ensure valid geometries
species1_poly <- st_make_valid(species1_poly)
species2_poly <- st_make_valid(species2_poly)

# Calculate intersection
overlap_poly <- st_intersection(species1_poly, species2_poly)

# Plot with overlap ------------------------------------------------------------
fviz_pca_ind(traitsPCA,
             label = "none",
             alpha.ind = 0.1) +
  geom_sf(data = species1_poly,  inherit.aes = FALSE, alpha = 0.6, color = NA, aes(fill = Species1)) +
  geom_sf(data = species2_poly,  inherit.aes = FALSE, alpha = 0.6, color = NA, aes(fill = Species2)) +
  geom_sf(data = overlap_poly, inherit.aes = FALSE, alpha = 0.6, color = NA, fill = "#8849a5") +
  scale_fill_manual(values = c("#25b9be", "#e41a1c")) +
  #geom_segment(data = var_data[var_data$label %in% selected_vars,],
   #            aes(x = 0, y = 0, xend = PC1 * 14, yend = PC2 * 14),  
    #           arrow = arrow(length = unit(0.2, "cm")),
     #          color = "#4682b4") +
  #geom_text(data = var_data[var_data$label %in% selected_vars,],
   #         aes(x = PC1 * 14, y = PC2 * 14, label = label),
    #        color = "#4682b4", vjust = -0.8) +
  #geom_text_repel(data = landcover_centroids, aes(x = x, y = y, label = Groups), 
   #               color = "black", point.padding = 0) +
  theme_minimal() +
  coord_sf() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  labs(title = "")

ggsave("test.pdf", height=7, width=7)


# dark version plot overlap
fviz_pca_ind(traitsPCA,
             label = "none",
             alpha.ind = 0.1,
             col.ind="white") +
  geom_sf(data = species1_poly,  inherit.aes = FALSE, alpha = 0.6, color = NA, aes(fill = Species1)) +
  geom_sf(data = species2_poly,  inherit.aes = FALSE, alpha = 0.6, color = NA, aes(fill = Species2)) +
  geom_sf(data = overlap_poly, inherit.aes = FALSE, alpha = 0.6, color = NA, fill = "#8849a5") +
  scale_fill_manual(values = c("#25b9be", "#e41a1c")) +
  geom_segment(data = var_data[var_data$label %in% selected_vars,],
               aes(x = 0, y = 0, xend = PC1 * 14, yend = PC2 * 14),  
               arrow = arrow(length = unit(0.2, "cm")),
               color = "#73b1e6") +
  geom_text(data = var_data[var_data$label %in% selected_vars,],
            aes(x = PC1 * 14, y = PC2 * 14, label = label),
            color = "#73b1e6", vjust = -0.8) +
  geom_text_repel(data = landcover_centroids, aes(x = x, y = y, label = Groups), 
                  color = "white", point.padding = 0) +
  dark_theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  labs(title = "")


# make a not warped version ----------------------------------------------------

species <- colnames(roundedCounts)[5:96]
sort(species)
Species1 <- "Yellow.rumped.Warbler" 
Species2 <- "Red.breasted.Nuthatch" 

# get polygons -----------------------------------------------------------------
species1_poly <- speciesPolygons %>% filter(species == Species1)
species2_poly <- speciesPolygons %>% filter(species == Species2)

# Ensure valid geometries
species1_poly <- st_make_valid(species1_poly)
species2_poly <- st_make_valid(species2_poly)

# Calculate intersection
overlap_poly <- st_intersection(species1_poly, species2_poly)

# Function to convert sf POLYGON to data frame for ggplot
sf_to_df <- function(sf_obj) {
  sf::st_geometry(sf_obj) %>%
    sf::st_cast("POLYGON") %>%
    lapply(sf::st_coordinates) %>%
    lapply(as.data.frame) %>%
    bind_rows(.id = "group")
}

species1_df <- sf_to_df(species1_poly)
species2_df <- sf_to_df(species2_poly)
overlap_df  <- sf_to_df(overlap_poly)

fviz_pca_ind(traitsPCA,
             label = "none",
             alpha.ind = 0.1) +
  geom_polygon(data = species1_df, aes(x = X, y = Y, group = group), fill = "#25b9be", alpha = 0.6) +
  geom_polygon(data = species2_df, aes(x = X, y = Y, group = group), fill = "#e41a1c", alpha = 0.6) +
  geom_polygon(data = overlap_df, aes(x = X, y = Y, group = group), fill = "#8849a5", alpha = 0.6) +
  annotate("text", x = -7, y = -8, label = "Oak") +
  annotate("text", x =  6, y = -2, label = "Pine") +
  annotate("text", x =  0, y = 1, label = "Mixed Conifer") +
  coord_fixed(ratio = 17.2/27.4) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  labs(title = "")

ggsave("OverlapPoly_YRWarb_RBNut.pdf", height=5, width=7)

################################################################################
# Calculate all pairwise overlaps
################################################################################

# define species combinations --------------------------------------------------
pairwise_combinations <- expand_grid(species1 = speciesPolygons$species,
                                     species2 = speciesPolygons$species) %>% filter(species1 < species2)

# make fxn to calculate overlaps
calculate_overlaps <- function(sp1, sp2) {
  # get polygons -----------------------------------------------------------------
  species1_poly <- speciesPolygons %>% filter(species == sp1)
  species2_poly <- speciesPolygons %>% filter(species == sp2)
  
  # Ensure valid geometries
  species1_poly <- st_make_valid(species1_poly)
  species2_poly <- st_make_valid(species2_poly)
  
  # Calculate intersection
  overlap_poly <- st_intersection(species1_poly, species2_poly)
  
  # calculate overlap area
  overlap_area <- st_area(overlap_poly)
  # Return the sf object
  if (nrow(overlap_poly) == 0) {return(0)}
  return(overlap_area)
}


# Calculate intersection areas -------------------------------------------------
pairwise_combinations$overlap_area <- mapply(calculate_overlaps,
                                             pairwise_combinations$species1,
                                             pairwise_combinations$species2)


# save the dataframe -----------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Draft 1")
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Shiny/Avispace4.0")

write.csv(pairwise_combinations, "speciesNicheOverlapAreas_jittered_detPerWeek_20250624.csv")



################################################################################
# Species niche Areas & Percent overlap
################################################################################

# species areas
speciesAreas <- data.frame(species = speciesPolygons$species,
                           area = st_area(speciesPolygons))

# save it 
write_rds(speciesAreas, "speciesNicheAreas_20250624.rds")

# add areas to pairwise overlap
percentOverlap <- merge(pairwise_combinations, speciesAreas, by.x = "species1", by.y = "species")
names(percentOverlap)[names(percentOverlap) == 'area'] <- 'area1'
percentOverlap <- merge(percentOverlap, speciesAreas, by.x = "species2", by.y = "species")
names(percentOverlap)[names(percentOverlap) == 'area'] <- 'area2'

# calculate % overlap
percentOverlap$percentOverlap <- percentOverlap$overlap_area / (percentOverlap$area1 + percentOverlap$area2 - percentOverlap$overlap_area)

# save it 
write_rds(percentOverlap, "percentOverlap_20250624.rds")



################################################################################
# calculate Range overlaps
################################################################################


# toy - single species pair ----------------------------------------------------

# read in data
x <- st_read("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps/Yellow-rumped Warbler/yerwar_range_2022.gpkg")
mapview(x)
y <- st_read("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps/Red-breasted Nuthatch/rebnut_range_2022.gpkg")
mapview(y)

# keep only range in breeding season
x <- subset(x, x$season=="resident" | x$season=="breeding")
y <- subset(y, y$season=="resident" | y$season=="breeding")
mapview(x)
mapview(y)

# what is the CRS?
st_crs(x) #WGS84, EPSG4326
st_crs(y)

# lat/lon doesn't give accurate area measurements
# switch to Albers Equal Area 5070
x <- st_transform(x, crs = 5070)
y <- st_transform(y, crs = 5070)

# measure area
st_area(x)
st_area(y)

# get overlap area
z <- st_intersection(x,y)
mapview(z)
st_area(z)

# get percent overlap
w <- st_area(z) / (st_area(x) + st_area(y) - st_area(z))
w

# visualize
mapview(x) +
  mapview(y)


# LOOP - all species pairs -----------------------------------------------------

# temporarily chage wd
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps")

# Get all folder names (aka species names) in the current working directory
species_paths <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# Create an empty list to store results
species_data <- list()  

# loop through to load all .gpkg files
for (species in species_paths) {
  files <- list.files(path = species, pattern = ".*_range_(2022|2023)\\.gpkg$", full.names = TRUE)
  
  if (length(files) > 0) {
    species_data[[species]] <- st_read(files[1])  # Store in the list with species name as key
    print(paste("Loaded:", files[1]))
  } else {
    print(paste("No matching file found for species:", species))
  }
}

# Initialize an empty data frame to store results
rangeOverlaps <- data.frame(
  species1 = character(),
  species2 = character(),
  area1 = numeric(),
  area2 = numeric(),
  overlap_area = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all pairwise combinations of species
for (i in 1:(length(species_data) - 1)) {
  for (j in (i + 1):length(species_data)) {
    
    species1 <- str_extract(names(species_data)[i], "[^/]+$")
    species2 <- str_extract(names(species_data)[j], "[^/]+$")
    
    poly1 <- species_data[[names(species_data)[i]]]
    poly2 <- species_data[[names(species_data)[j]]]
    
    # Keep only breeding + resident
    poly1 <- subset(poly1, season %in% c("resident", "breeding"))
    poly2 <- subset(poly2, season %in% c("resident", "breeding"))
    
    # Make valid and reproject
    poly1 <- st_transform(st_make_valid(poly1), 5070)
    poly2 <- st_transform(st_make_valid(poly2), 5070)
    
    # Try intersection and area calculation
    tryCatch({
      area1 <- sum(st_area(poly1))
      area2 <- sum(st_area(poly2))
      overlap <- st_intersection(poly1, poly2)
      overlap_area <- sum(st_area(overlap))
      
      rangeOverlaps <- rbind(rangeOverlaps, data.frame(
        species1 = species1,
        species2 = species2,
        area1 = as.numeric(area1),
        area2 = as.numeric(area2),
        overlap_area = as.numeric(overlap_area)
      ))
    }, error = function(e) {
      message(paste("Failed:", species1, "-", species2, "|", e$message)) # E.C. Dove gives some errors...
    })
  }
}


# calculate range overlap percent
rangeOverlaps$overlap_percent <- rangeOverlaps$overlap_area / (rangeOverlaps$area1 + rangeOverlaps$area2 - rangeOverlaps$overlap_area)

# go back to main wd
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Draft 1")

# write csv

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Shiny/Avispace4.0")

write.csv(rangeOverlaps, "rangeOverlaps_20250624.csv")
write_rds(rangeOverlaps, "rangeOverlaps_20250624.rds")

# make dot version -------------------------------------------------------------
dotVclean_nameMatches <-  read.csv("dotVclean_nameMatches.csv")

rangeOverlapsDot <- merge(rangeOverlaps, dotVclean_nameMatches, by.x = "species1", by.y = "cleanName")
rangeOverlapsDot$species1 <- NULL
names(rangeOverlapsDot)[names(rangeOverlapsDot) == "dotName"] <- "species1"

rangeOverlapsDot <- merge(rangeOverlapsDot, dotVclean_nameMatches, by.x = "species2", by.y = "cleanName")
rangeOverlapsDot$species2 <- NULL
names(rangeOverlapsDot)[names(rangeOverlapsDot) == "dotName"] <- "species2"

# save dot version -------------------------------------------------------------

write.csv(rangeOverlapsDot, "rangeOverlaps_dotNames_20250624.csv")
write_rds(rangeOverlapsDot, "rangeOverlaps_dotNames_20250624.rds")


################################################################################
# Range Overlap Maps
################################################################################


# read in data
x <- st_read("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps/Yellow-rumped Warbler/yerwar_range_2022.gpkg")
mapview(x)
y <- st_read("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps/Red-breasted Nuthatch/rebnut_range_2022.gpkg")
mapview(y)
z <- st_read("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps/Purple Finch/purfin_range_2022.gpkg")
mapview(z)


# keep only range in breeding season
x <- subset(x, x$season=="resident" | x$season=="breeding")
y <- subset(y, y$season=="resident" | y$season=="breeding")
z <- subset(z, z$season=="resident" | z$season=="breeding")

# get intersection
w <- st_intersection(x, y)

# Get North America map
north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")

# Define a bounding box (xmin, ymin, xmax, ymax)
bbox <- st_bbox(c(xmin = -140, xmax = -50, ymin = 5, ymax = 70), crs = st_crs(north_america))

# Clip to the bounding box
north_america_clipped <- st_crop(north_america, bbox)
w_clipped <- st_crop(w, bbox)
x_clipped <- st_crop(x, bbox)
y_clipped <- st_crop(y, bbox)
z_clipped <- st_crop(z, bbox)

# plot it
ggplot() +
  geom_sf(data = north_america_clipped, fill = "grey40", color = "grey40") +
  geom_sf(data = y_clipped, fill = "#25b9be", color = "#25b9be", size = 1, alpha=0.6) +
  geom_sf(data = x_clipped, fill = "#e41a1c", color = "#e41a1c", size = 1, alpha=0.6) +
  geom_sf(data = w_clipped, fill = "#8849a5", color = "#8849a5", size = 1, alpha=0.6) +
  #geom_sf(data = z_clipped, fill = "#00ba39", color = "#00ba39", size = 1, alpha=0.6) +
  theme_void() +  
  theme(legend.position = "none")

ggsave("test_20250522.pdf", height = 5, width = 5)

# plot dark version
ggplot() +
  geom_sf(data = north_america_clipped, fill = "grey40", color = "grey40") +
  geom_sf(data = x_clipped, fill = "#25b9be", color = "#25b9be", size = 1, alpha=0.6) +
  geom_sf(data = y_clipped, fill = "#e41a1c", color = "#e41a1c", size = 1, alpha=0.6) +
  geom_sf(data = w_clipped, fill = "#8849a5", color = "#8849a5", size = 1, alpha=0.6) +
  #geom_sf(data = z_clipped, fill = "#00ba39", color = "#00ba39", size = 1, alpha=0.6) +
  dark_theme_void() +  
  theme(legend.position = "none")


# make a shiny-ready flexible version of the range map -------------------------

# pick your players
species <- colnames(roundedCounts)[5:96]
sort(species)
Species1 <- "California.Scrub.Jay"
Species2 <- "Woodhouse.s.Scrub.Jay"

# find the paths
map_paths <- list.files("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/RangeMaps/", 
                        pattern = ".*_range_(2022|2023)\\.gpkg$", full.names = TRUE, recursive = TRUE)

# get the maps
Species1_map <- st_read(grep(Species1, map_paths, value=TRUE))
Species2_map <- st_read(grep(Species2, map_paths, value=TRUE))

# only breeding season please
Species1_map_breeding <- subset(Species1_map, Species1_map$season=="resident" | Species1_map$season=="breeding")
Species2_map_breeding <- subset(Species2_map, Species2_map$season=="resident" | Species2_map$season=="breeding")

# get intersection
intersection_map <- st_intersection(Species1_map_breeding, Species2_map_breeding)

# Get North America map
north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")

# Define a bounding box (xmin, ymin, xmax, ymax)
bbox <- st_bbox(c(xmin = -140, xmax = -50, ymin = 5, ymax = 70), crs = st_crs(north_america))

# Clip to the bounding box
north_america_clipped <- st_crop(north_america, bbox)
Species1_map_clipped <- st_crop(Species1_map_breeding, bbox)
Species2_map_clipped <- st_crop(Species2_map_breeding, bbox)
intersection_map_clipped <- st_crop(intersection_map, bbox)

# plot it
ggplot() +
  geom_sf(data = north_america_clipped, fill = "grey40", color = "grey40") +
  geom_sf(data = Species1_map_clipped, fill = "#25b9be", color = "#25b9be", size = 1, alpha=0.6) +
  geom_sf(data = Species2_map_clipped, fill = "#e41a1c", color = "#e41a1c", size = 1, alpha=0.6) +
  geom_sf(data = intersection_map_clipped, fill = "#8849a5", color = "#8849a5", size = 1, alpha=0.6) +
  theme_void() +  
  theme(legend.position = "none")

# save it
ggsave("RangeOverlap_OakTit_RNSap.pdf", height=7, width=7)

################################################################################
# Combine range overlaps and niche overlaps
################################################################################

# ensure column names are informative ------------------------------------------
NicheOverlap <- percentOverlap[, c("species1", "species2", "percentOverlap")]
names(NicheOverlap)[names(NicheOverlap) == 'percentOverlap'] <- 'niche_overlap_percent'

rangeOverlaps <- rangeOverlaps[, c("species1", "species2", "overlap_percent")]
names(rangeOverlaps)[names(rangeOverlaps) == "overlap_percent"] <- "range_overlap_percent" 


# match species name formats ---------------------------------------------------
dotVclean_nameMatches <-  read.csv("dotVclean_nameMatches.csv")

NicheOverlap <- merge(NicheOverlap, dotVclean_nameMatches, by.x = "species1", by.y = "dotName")
NicheOverlap$species1 <- NULL
names(NicheOverlap)[names(NicheOverlap) == "cleanName"] <- "species1"

NicheOverlap <- merge(NicheOverlap, dotVclean_nameMatches, by.x = "species2", by.y = "dotName")
NicheOverlap$species2 <- NULL
names(NicheOverlap)[names(NicheOverlap) == "cleanName"] <- "species2"

# merge ------------------------------------------------------------------------

# Merge when species1 = species1 and species2 = species2
merge1 <- merge(NicheOverlap, rangeOverlaps, by = c("species1", "species2"))

# Merge when species1 = species2 and species2 = species1
merge2 <- merge(NicheOverlap, rangeOverlaps, by.x = c("species1", "species2"), by.y = c("species2", "species1"))

# Combine both merged dataframes
NicheRangeOverlap <- rbind(merge1, merge2)

################################################################################
# plot range vs niche overlap
################################################################################

# mark congeneric pairs --------------------------------------------------------

# load scientific names
SciNames <- read.csv("~/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/species_scientificNames.csv")
SciNames <- SciNames[, c("Species", "genus")]

# Join to get genus for species1 and species2
NicheRangeOverlap <- NicheRangeOverlap %>% left_join(SciNames, by = c("species1" = "Species"))
NicheRangeOverlap <- NicheRangeOverlap %>% rename("genus_species1" = "genus")

NicheRangeOverlap <- NicheRangeOverlap %>% left_join(SciNames, by = c("species2" = "Species"))
NicheRangeOverlap <- NicheRangeOverlap %>% rename("genus_species2" = "genus")

# Create the 'congeneric' column to check if species1 and species2 belong to the same genus
NicheRangeOverlap <- NicheRangeOverlap %>%
  mutate(congeneric = if_else(genus_species1 == genus_species2, "Same Genus", "Different Genus"))

# drop units from range
NicheRangeOverlap$range_overlap_percent <- drop_units(NicheRangeOverlap$range_overlap_percent)

# add sister species -----------------------------------------------------------

NicheRangeOverlap$sisterSpecies <- ifelse(NicheRangeOverlap$species1 == "California Scrub-Jay" & NicheRangeOverlap$species2 == "Woodhouse's Scrub-Jay" |
                                            NicheRangeOverlap$species1 == "Red-breasted Sapsucker" & NicheRangeOverlap$species2 == "Red-naped Sapsucker" |
                                            NicheRangeOverlap$species2 == "Plumbeous Vireo" & NicheRangeOverlap$species1 == "Cassin's Vireo" |
                                            NicheRangeOverlap$species2 == "Purple Finch" & NicheRangeOverlap$species1 == "Cassin's Finch",
                                          TRUE, FALSE)

# save it ----------------------------------------------------------------------
write.csv(NicheRangeOverlap, "NicheRangeOverlap_cleanNames_20250624.csv")
write_rds(NicheRangeOverlap, "NicheRangeOverlap_cleanNames_20250624.rds")

# make dot names version -------------------------------------------------------

NicheRangeOverlapDot <- merge(NicheRangeOverlap, dotVclean_nameMatches, by.x = "species1", by.y = "cleanName")
NicheRangeOverlapDot$species1 <- NULL
names(NicheRangeOverlapDot)[names(NicheRangeOverlapDot) == "dotName"] <- "species1"

NicheRangeOverlapDot <- merge(NicheRangeOverlapDot, dotVclean_nameMatches, by.x = "species2", by.y = "cleanName")
NicheRangeOverlapDot$species2 <- NULL
names(NicheRangeOverlapDot)[names(NicheRangeOverlapDot) == "dotName"] <- "species2"

# save dot version -------------------------------------------------------------

write.csv(NicheRangeOverlapDot, "NicheRangeOverlap_dotNames_20250624.csv")
write_rds(NicheRangeOverlapDot, "NicheRangeOverlap_dotNames_20250624.rds")

# plot it ----------------------------------------------------------------------



# light
ggplot(NicheRangeOverlap, aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_smooth(aes(fill = congeneric), method = "lm", alpha = 0.2) +  
  theme_minimal() +
  scale_color_manual(values = c("purple", "darkgreen")) +
  scale_fill_manual(values = c("purple", "darkgreen")) +  
  geom_point(data=subset(NicheRangeOverlap, NicheRangeOverlap$sisterSpecies==TRUE), color="black") +
  theme(legend.title = element_blank()) +
  stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x, parse = TRUE) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)")
ggsave("NicheRangeOverlap_light_20250624.pdf", height=5, width=7)

# dark
ggplot(NicheRangeOverlap, aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_smooth(aes(fill = congeneric), method = "lm", alpha = 0.3) +  
  dark_theme_minimal() +
  scale_color_manual(values = c("purple", "#A0C878")) +
  scale_fill_manual(values = c("purple", "#A0C878")) +  
  geom_point(data=subset(NicheRangeOverlap, NicheRangeOverlap$sisterSpecies==TRUE), color="white") +
  theme(legend.title = element_blank()) +
  stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x, parse = TRUE) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)")
ggsave("NicheRangeOverlap_dark_20250515.pdf", height=5, width=7)



ggplot(NicheRangeOverlap, aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_point() +
  facet_grid(~congeneric)
ggsave("NicheRangeOverlap_points_20250624.pdf", height=5, width=7)


################################################################################
# ANCOVA
################################################################################

m1 <- lm(range_overlap_percent ~ niche_overlap_percent * congeneric, data = NicheRangeOverlap)
summary(m1)
anova(m1)

################################################################################
# Figure 2 niche areas 
################################################################################

ggplot(NicheAreas, aes(x = reorder(species, area), y = area)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip the coordinates to make it a horizontal bar chart
  labs(
    x = "Species",
    y = "Niche Area",
    title = "Species Niche Areas Sorted by Size"
  ) +
  theme_minimal()

ggsave("NicheAreas_20250709.pdf", height=12, width=7)

################################################################################
# Figure 3
################################################################################

# Heirarchical Matrix

# fill out the gaps in overlaps ------------------------------------------------
selfOverlap <- data_frame(species1 = speciesAreas$species,
                          species2 = speciesAreas$species,
                          niche_overlap_percent = 1)

backwards <- data_frame(species1 = NicheRangeOverlapDot$species2,
                      species2 = NicheRangeOverlapDot$species1,
                      niche_overlap_percent = NicheRangeOverlapDot$niche_overlap_percent)

forwards <- data_frame(species1 = NicheRangeOverlapDot$species1,
                      species2 = NicheRangeOverlapDot$species2,
                      niche_overlap_percent = NicheRangeOverlapDot$niche_overlap_percent)

alloverlaps <- rbind(selfOverlap, backwards, forwards)

# fill gaps
alloverlaps <- alloverlaps %>%
  complete(species1, species2, fill = list(niche_overlap_percent = 0))

# make an overlap matrix
overlap_matrix <- acast(alloverlaps, species1 ~ species2, value.var = "niche_overlap_percent", fill = 0)

# heirarchical clustering
hc <- hclust(as.dist(1 - overlap_matrix))
species_order <- hc$labels[hc$order]

# Reorder species1 and species2 based on desired order
alloverlaps$species1 <- factor(alloverlaps$species1, levels = species_order)
alloverlaps$species2 <- factor(alloverlaps$species2, levels = species_order)

# plot percent overlap heatmap
ggplot(alloverlaps, aes(x=species1, y=species2, fill=niche_overlap_percent)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  scale_fill_viridis_c() +
  coord_trans(y = "reverse") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))

# bar plot ---------------------------------------------------------------------

# calculate species centroids
speciesCentroids <- st_centroid(speciesPolygons)
spPC <- as.data.frame(st_coordinates(speciesCentroids))
speciesCentroids <- cbind(speciesCentroids, spPC)

# merge with areas
speciesAreas_centroids <- merge(speciesAreas, speciesCentroids)
  
# reorder
speciesAreas_centroids$species <- factor(speciesAreas_centroids$species, levels = species_order)

# set color scheme
speciesAreas_centroids$color <- rgb(scales::rescale(speciesAreas_centroids$X), scales::rescale(speciesAreas_centroids$Y), 0.5)

# plot bar plot
ggplot(speciesAreas_centroids, aes(y = species, x = area, fill=X)) +
  geom_bar(stat="summary") +
  scale_y_discrete(limits = rev(levels(speciesAreas_centroids$species))) +
  xlab("Niche Area") +
  scale_fill_viridis_c(option="A") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

# save Heatmap plot
heatmap_plot <- ggplot(alloverlaps, aes(x=species1, y=species2, fill=niche_overlap_percent)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  scale_fill_viridis_c() +
  coord_trans(y = "reverse") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0), legend.position = "none") +
  xlab("") +
  ylab("")


# save Bar plot for NicheArea2
bar_plot <- ggplot(speciesAreas_centroids, aes(y = species, x = area, fill=X)) +
  geom_bar(stat="summary") +
  scale_y_discrete(limits = rev(levels(speciesAreas_centroids$species))) +
  xlab("Niche Area") +
  scale_fill_viridis_c(option="A") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")

# Align plots
combined_plot <- plot_grid(
  heatmap_plot,
  bar_plot,
  align = "h",
  axis = "tblr",
  rel_widths = c(8, 1)
)

combined_plot



################################################################################
# Species details supp table
################################################################################

# niche area
# total # detections
# total # sites

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Draft 2")


speciesDetails <- merge(speciesAreas, speciesDetections)
colnames(speciesDetails)[colnames(speciesDetails) == 'area'] <- 'nicheArea'
speciesDetails <- merge(speciesDetails, speciesSites)
write_rds(speciesDetails, "speciesDetails_20250708.rds")
write_csv(speciesDetails, "speciesDetails_20250708.csv")





