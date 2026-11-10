
# full protocol for processing Sierra BioCube and Acoustics data
# Laura Berman
# created 30 April 2025
# last updated 04 November 2026

# load packages ----------------------------------------------------------------

library(tidyverse)
library(vroom)
library(sf)
library(mapview)
library(raster)
library(ggspatial)
library(corrplot)
library(factoextra)
library(ggrepel)
library(rnaturalearth)
#library(ggdark)
library(units)
library(ggpmisc)
library(basemaps)
library(reshape2)
library(cowplot)
library(geosphere)
library(scales)
library(patchwork)
library(MCMCglmm)
library(coda)
library(mgcv)

# Table of contents ------------------------------------------------------------
# 1) Processing ARU bioacoustic data - 73
# 2) Species complexes and cross-identification - 146
# 3) Adjust detection rates for day-of-year and year - 243
# 4) Site buffers - 473
# 5) Remove burned areas - 500
# 6) Get 1 km habitat variables - 567
# 7) add USGS CONUS GAP 2011 landcover type - 676
# 8) which points are too close to each other? - 708
# 9) remove species with insufficient data - 778
# 10) make a nice map - 830
# 11) check for colinearity - 868
# 12) Niche Space PCA - 978
# 13) PCA loadings (supp figures) - 998
# 14) render landcover labels - 1067
# 15) test plot for errored species polygons - 1150
# 16) Species Niche Polygons - 1189
# 17) plot 2 species niche overlap - 1288
# 18) Calculate all pairwise overlaps - 1416
# 19) Species niche Areas & Percent overlap - 1460
# 20) calculate Range overlaps - 1486
# 21) Range Overlap Maps - 1635
# 22) Combine range overlaps and niche overlaps - 1740
# 23) plot range vs niche overlap - 1774
# 24) Mixed effects model - 1866
# 25) permutation test - 2020
# 26) Figure 4 niche areas - 2063
# 27) Figure 5 - matrix - 2079
# 28) Species details supp table - 2218


# working directory  -----------------------------------------------------------
# base_dir holds the large input data (detections, mosaics, range maps, etc.)
# that are too big for the GitHub repository. Change this one line to run
# the protocol on another machine.
base_dir <- "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics"

setwd(file.path(base_dir, "Draft 1"))


################################################################################
# 1) Processing ARU bioacoustic data 
################################################################################

# load ARU locations -----------------------------------------------------------

ARUlocations <- read.csv(file.path(base_dir, "laura_berman_diurnal_birds/data-clean/focal_aru_metadata.csv")) # 2100 locations
ARUlocations <- ARUlocations %>% distinct(survey_year, cell_unit, box, .keep_all = TRUE) # 2073 unique locations

# load detections --------------------------------------------------------------

# only cell_unit, species_code, survey_date, filename and relative_sec are
# used downstream, so read just those columns (filename and relative_sec
# identify the BirdNET segment, for the cross-identification check in 1.1). vroom reads and row-binds all three files in one call.
# species_code is read as integer to match the species_codes_thresholds key.
# full columns: "filename" "cell_unit" "year" "survey_date" "month" "day" "hour" "min" "second" "relative_sec" "species_code" "logit_score"
det <- vroom(file.path(base_dir, "laura_berman_diurnal_birds/data-clean",
                       c("2021_filtered_detections.csv",    # 41,264,999 obs
                         "2022_filtered_detections.csv",    # 40,424,580 obs
                         "2023_filtered_detections.csv")),  # 27,057,557 obs
             delim = ",",
             col_types = cols_only(cell_unit = col_character(),
                                   survey_date = col_character(),
                                   filename = col_character(),
                                   relative_sec = col_double(),
                                   species_code = col_integer()))

# nrow(det) should be 108,747,136 # correct


# load common names ------------------------------------------------------------
# this file also has BirdNET confidence thresholds
# and includes species not actually in the dataset
# a threshold of 100 is a placeholder for species who never occur

commonNames <- read.csv(file.path(base_dir, "laura_berman_diurnal_birds/data-clean/species_codes_thresholds.csv")) # 241 species


# load effort ------------------------------------------------------------------

# only cell_unit, date and aru_hours are used downstream. date is kept as
# character (as read.csv returned it) so startDate/endDate in unit_effort
# keep the same format.
effort <- vroom(file.path(base_dir, "laura_berman_diurnal_birds/data-clean",
                          c("2021_filtered_effort.csv",    # 134,322 obs
                            "2022_filtered_effort.csv",    # 127,002 obs
                            "2023_filtered_effort.csv")),  # 110,349 obs
                delim = ",",
                col_types = cols_only(cell_unit = col_character(),
                                      date = col_character(),
                                      aru_hours = col_double()))
# nrow(effort) should be 371,673 # correct
# effort includes rows for days with 0 effort
# nrow(subset(effort, effort$aru_hours==7)) = 66,349 obs = total full recording days (subsampling records 7 hrs per day)
# nrow(subset(effort, aru_hours != 0 & aru_hours != 7)) = 3,975 partial days. We'll exclude these

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

write_rds(unit_effort, "unit_effort.rds")

################################################################################
# 2) Species complexes and cross-identification
################################################################################

# cross-identification check ---------------------------------------------------
# 05 Nov 2026. Candidate groups of species whose vocalizations BirdNET may
# confuse. If two members of a group are detected in the same BirdNET segment
# (same site, same file, same relative_sec), it is likely one vocalization
# labelled as both. These counts are reported in the supplement and justify
# which species are combined below.

crossIDgroups <- data.frame(member = c("Vireo cassinii", "Vireo plumbeus",
                                       "Sphyrapicus ruber", "Sphyrapicus nuchalis",
                                       "Sphyrapicus thyroideus"),
                            group = c("Vireo", "Vireo",
                                      "Sphyrapicus", "Sphyrapicus", "Sphyrapicus"))
crossIDgroups$species_code <- commonNames$species_code[match(crossIDgroups$member, commonNames$scientific_name)]
stopifnot(!anyNA(crossIDgroups$species_code))

memberDet <- det %>%
  filter(species_code %in% crossIDgroups$species_code) %>%
  dplyr::select(cell_unit, filename, relative_sec, survey_date, species_code) %>%
  mutate(group = crossIDgroups$group[match(species_code, crossIDgroups$species_code)]) %>%
  group_by(group, cell_unit, filename, relative_sec) %>%
  mutate(n_members = n_distinct(species_code)) %>%
  ungroup()

# per species: % of its detections that share a segment with another group member
crossID <- memberDet %>%
  group_by(group, species_code) %>%
  dplyr::summarise(total_detections = n(),
                   shared_detections = sum(n_members > 1),
                   .groups = "drop") %>%
  mutate(pct_shared = 100 * shared_detections / total_detections) %>%
  left_join(commonNames[c("species_code", "scientific_name", "common_name")], by = "species_code")
crossID

# per species pair: number of segments where both were detected
crossIDpairs <- memberDet %>%
  dplyr::select(group, cell_unit, filename, relative_sec, species_code) %>%
  inner_join(., ., by = c("group", "cell_unit", "filename", "relative_sec"),
             suffix = c("_1", "_2"), relationship = "many-to-many") %>%
  filter(species_code_1 < species_code_2) %>%
  dplyr::count(group, species_code_1, species_code_2, name = "shared_segments") %>%
  mutate(species_1 = commonNames$common_name[match(species_code_1, commonNames$species_code)],
         species_2 = commonNames$common_name[match(species_code_2, commonNames$species_code)])
crossIDpairs

write_csv(crossID, "crossIdentification_species_20261105.csv")
write_csv(crossIDpairs, "crossIdentification_pairs_20261105.csv")
rm(memberDet)

# species complexes ------------------------------------------------------------
# Based on the check above, only Red-breasted and Red-naped Sapsucker share a
# substantial share of segments (~25% and ~35% of their detections).
# Williamson's Sapsucker (0.5%) and Cassin's/Plumbeous Vireo (0.1% / 2.6%)
# are kept separate.
# The complex is added as an extra taxon. The member species are also kept
# with their own (unmerged) detections, so later comparisons can use either.

complexNames <- data.frame(common_name = "Red-breasted/Red-naped Sapsucker",
                           scientific_name = "Sphyrapicus ruber/nuchalis")
complexNames$species_code <- max(commonNames$species_code) + seq_len(nrow(complexNames))

complexMembers <- data.frame(member = c("Sphyrapicus ruber", "Sphyrapicus nuchalis"),
                             common_name = "Red-breasted/Red-naped Sapsucker")
complexMembers$member_code <- commonNames$species_code[match(complexMembers$member, commonNames$scientific_name)]
complexMembers$complex_code <- complexNames$species_code[match(complexMembers$common_name, complexNames$common_name)]
stopifnot(!anyNA(complexMembers$member_code))

complexDet <- det %>%
  filter(species_code %in% complexMembers$member_code) %>%
  dplyr::select(cell_unit, filename, relative_sec, survey_date, species_code) %>%
  mutate(complex_code = complexMembers$complex_code[match(species_code, complexMembers$member_code)])

# per complex: detections, segments, and duplicates removed
crossIDcomplex <- complexDet %>%
  group_by(complex_code, cell_unit, filename, relative_sec) %>%
  dplyr::summarise(n_members = n_distinct(species_code), n = n(), .groups = "drop") %>%
  group_by(complex_code) %>%
  dplyr::summarise(member_detections = sum(n),
                   segments = n(),
                   shared_segments = sum(n_members > 1)) %>%
  mutate(duplicates_removed = member_detections - segments,
         pct_removed = 100 * duplicates_removed / member_detections) %>%
  left_join(complexNames[c("species_code", "common_name")], by = c("complex_code" = "species_code"))
crossIDcomplex
write_csv(crossIDcomplex, "crossIdentification_complexes_20261105.csv")

# keep one detection per segment, labelled as the complex
complexDet <- complexDet %>%
  distinct(complex_code, cell_unit, filename, relative_sec, .keep_all = TRUE) %>%
  mutate(species_code = complex_code)

# add the complex to the species name table
commonNames <- bind_rows(commonNames, complexNames)

################################################################################
# 3) Adjust detection rates for day-of-year and year
################################################################################
# 04 Nov 2026 - reviewer request. Recorders were deployed on different dates
# and in different years, so a raw site average partly reflects when a site
# was sampled. For each species we fit the seasonal (doy) and annual curve of
# vocal activity on full 7 h recording days, then rescale each site's
# observed detections by the activity expected on the days it was sampled.
#
# model (per species, one row per site-day):
#   n_det ~ s(doy) + year + s(cell_unit, bs = "re"),  negative binomial
# site is a random effect so the curve is not confounded with where sites
# were sampled early vs late.
#
# adjusted rate for site s:
#   obs_per_day = sum(n_det_s) / sum(g(day_s)) * mean(g(reference days))
# where g() is the population-level curve (site effect excluded). This is
# the observed / expected ratio scaled to a common calendar. Sites with zero
# detections stay at zero. A site sampled on every reference day in every
# year keeps its raw mean daily count.

# full 7 h recording days ------------------------------------------------------

fullDays <- effort %>%
  filter(aru_hours == 7) %>%
  mutate(date = as.Date(date),
         year = factor(format(date, "%Y")),
         doy  = as.numeric(format(date, "%j"))) %>%
  dplyr::select(cell_unit, date, year, doy)
# nrow(fullDays) should be 66,349

# daily counts per site x species ----------------------------------------------
# count first (fast on character dates), then parse the much smaller table.
# all species as detected, plus the de-duplicated complex counts.

dailyCounts <- bind_rows(
  det %>%
    dplyr::count(cell_unit, species_code, survey_date, name = "n_det"),
  complexDet %>%
    dplyr::count(cell_unit, species_code, survey_date, name = "n_det")) %>%
  mutate(date = as.Date(lubridate::parse_date_time(survey_date,
                                                   orders = c("Ymd", "mdY", "dmY")))) %>%
  dplyr::select(-survey_date)
rm(det, complexDet)

# check: dates parsed and matched to effort days
sum(is.na(dailyCounts$date)) # should be 0
# share of detections that fall on full 7 h days (only these are used now)
inFull <- paste(dailyCounts$cell_unit, dailyCounts$date) %in%
  paste(fullDays$cell_unit, fullDays$date)
sum(dailyCounts$n_det[inFull]) / sum(dailyCounts$n_det) #98% full days

dailyCounts <- dailyCounts[inFull, ]
rm(inFull)


# which species have enough data to fit a curve? -------------------------------
# fit cut-offs: > 30 sites and >= 300 detections. Section 9 uses the same
# cut-offs on the final site set of each buffer version.

speciesTotals <- dailyCounts %>%
  group_by(species_code) %>%
  dplyr::summarise(n_sites = n_distinct(cell_unit),
                   n_detections = sum(n_det))
fitSpecies <- speciesTotals$species_code[speciesTotals$n_sites > 30 &
                                           speciesTotals$n_detections >= 300]
length(fitSpecies)

# reference calendar -----------------------------------------------------------
# days covered by most of the sampling (5th - 95th percentile of doy), all years

refDoy <- seq(floor(quantile(fullDays$doy, 0.05)), ceiling(quantile(fullDays$doy, 0.95)))
refGrid <- expand.grid(doy = refDoy, year = levels(fullDays$year))
refGrid$year <- factor(refGrid$year, levels = levels(fullDays$year))
refGrid$cell_unit <- factor(fullDays$cell_unit[1], levels = sort(unique(fullDays$cell_unit)))

# fit one species --------------------------------------------------------------

fit_phenology <- function(sp) {
  d <- fullDays %>%
    left_join(dailyCounts[dailyCounts$species_code == sp, c("cell_unit", "date", "n_det")],
              by = c("cell_unit", "date")) %>%
    mutate(n_det = replace_na(n_det, 0L),
           cell_unit = factor(cell_unit, levels = levels(refGrid$cell_unit)))
  m <- bam(n_det ~ s(doy, k = 10) + year + s(cell_unit, bs = "re"),
           family = nb(), data = d, discrete = TRUE)

  # population-level curve (site effect excluded) on each sampled day
  d$g <- predict(m, newdata = d, type = "response", exclude = "s(cell_unit)")
  gRef <- mean(predict(m, newdata = refGrid, type = "response", exclude = "s(cell_unit)"))

  rates <- d %>%
    group_by(cell_unit) %>%
    dplyr::summarise(n_detections = sum(n_det),
                     active_days = n(),
                     expected = sum(g)) %>%
    mutate(cell_unit = as.character(cell_unit),
           species_code = sp,
           obs_per_day_raw = n_detections / active_days,
           obs_per_day = n_detections / expected * gRef) %>%
    filter(n_detections > 0)

  # curve for the supplementary figure (link scale CI, back-transformed)
  pr <- predict(m, newdata = refGrid, type = "link", se.fit = TRUE, exclude = "s(cell_unit)")
  curve <- refGrid %>%
    dplyr::select(doy, year) %>%
    mutate(species_code = sp,
           fit = exp(pr$fit),
           lwr = exp(pr$fit - 1.96 * pr$se.fit),
           upr = exp(pr$fit + 1.96 * pr$se.fit))

  sm <- summary(m)
  info <- data.frame(species_code = sp,
                     doy_edf = sm$s.table["s(doy)", "edf"],
                     doy_p = sm$s.table["s(doy)", "p-value"],
                     dev_expl = sm$dev.expl,
                     t(setNames(coef(m)[grep("^year", names(coef(m)))],
                                grep("^year", names(coef(m)), value = TRUE))))
  list(rates = rates, curve = curve, info = info)
}

# fit all species --------------------------------------------------------------

phenoFits <- map(setNames(fitSpecies, fitSpecies),
                 ~ tryCatch(fit_phenology(.x),
                            error = function(e) { message("species ", .x, ": ", conditionMessage(e)); NULL }),
                 .progress = TRUE)
failed <- names(phenoFits)[map_lgl(phenoFits, is.null)]
failed # species that did not fit fall back to the raw rate below

phenoRates <- map_dfr(phenoFits, "rates")
phenoCurves <- map_dfr(phenoFits, "curve")
phenoInfo <- map_dfr(phenoFits, "info")

# raw rates for species not fitted (too few data; removed later in section 9)

rawRates <- dailyCounts %>%
  filter(!species_code %in% phenoRates$species_code) %>%
  group_by(cell_unit, species_code) %>%
  dplyr::summarise(n_detections = sum(n_det), .groups = "drop") %>%
  merge(unit_effort[c("cell_unit", "active_days")]) %>%
  mutate(obs_per_day_raw = n_detections / active_days,
         obs_per_day = obs_per_day_raw)

siteDetections <- bind_rows(phenoRates[c("cell_unit", "species_code", "n_detections", "obs_per_day_raw", "obs_per_day")],
                            rawRates[c("cell_unit", "species_code", "n_detections", "obs_per_day_raw", "obs_per_day")])

# how much did the adjustment change site values?
cor(phenoRates$obs_per_day_raw, phenoRates$obs_per_day, method = "spearman") #0.99 spearman's correlation
summary(phenoRates$obs_per_day / phenoRates$obs_per_day_raw)

phenoInfo <- merge(phenoInfo, commonNames[c("species_code", "common_name")])
phenoInfo <- merge(phenoInfo, speciesTotals)
write_rds(phenoInfo, "phenologyModels_20261104.rds")
write_csv(phenoInfo, "phenologyModels_20261104.csv")
write_rds(phenoCurves, "phenologyCurves_20261104.rds")

# supplementary figure: species phenology curves -------------------------------
# predicted daily detections at an average site, by year

phenoCurves <- merge(phenoCurves, commonNames[c("species_code", "common_name")])

ggplot(phenoCurves, aes(x = as.Date(doy - 1, origin = "2021-01-01"), y = fit,
                        color = year, fill = year)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ common_name, scales = "free_y", ncol = 8) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(x = "Date", y = "Predicted detections per day (average site)",
       color = "Year", fill = "Year") +
  theme_minimal(base_size = 7) +
  theme(strip.text = element_text(size = 6),
        legend.position = "bottom")
ggsave("FigureS_phenologyCurves_20261104.pdf", height = 14, width = 11)

# combine detections, effort, and common names ---------------------------------

siteDetections <- merge(siteDetections, commonNames)
siteDetections <- merge(siteDetections, unit_effort)

# raw detection counts on 7 h days per site x species (used in section 9)
siteCounts <- siteDetections[c("cell_unit", "common_name", "n_detections")]
siteCounts$species <- make.names(siteCounts$common_name) # matches column names later on
write_rds(siteCounts, "siteCounts_20261105.rds")


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

hist(max_dists_by_cell$max_distance_m) # mostly very small
mean(max_dists_by_cell$max_distance_m) # avg 27m between redeployments. Close enough.

# combine years with the same cell_unit ----------------------------------------

ARUlocations <- ARUlocations %>% 
  group_by(cell_unit) %>%
  dplyr::summarise(long = mean(long),
                   lat = mean(lat))
#names(ARUlocations)[names(ARUlocations) == 'survey_year'] <- 'year'
siteDetections <- merge(siteDetections, ARUlocations)

# remove rows where active days < 30 -------------------------------------------

siteDetections <- subset(siteDetections, siteDetections$active_days > 30) # 59944 obs

# pivot wider ------------------------------------------------------------------
# [I will want active_days back later]
# obs_per_day is now the day-of-year and year adjusted rate from section 3
# (unadjusted rate kept in obs_per_day_raw)
siteDetections <- siteDetections[c("cell_unit", "common_name", "obs_per_day", "long", "lat")]
siteDetections <- siteDetections %>%
  group_by(cell_unit) %>%
  pivot_wider(names_from = common_name,
              values_from = obs_per_day) # 889 obs. one per site.
siteDetections[is.na(siteDetections)] <- 0



################################################################################
# 4) Site buffers
################################################################################
# 06 Nov 2026. Habitat variables and fire are summarized within a buffer around
# each ARU. 150 m is the main analysis (about 9 ha, the average home range of a
# small passerine). 50 m and 300 m versions are made for the buffer
# sensitivity test (reviewer request). From here on, each step is run on all
# three versions.
# Each site also carries the date of its last full recording day, used in
# section 5 so that a fire only counts if it started before that day.

ARUlocations.sf <- st_as_sf(ARUlocations, coords = c("long", "lat"), crs = 4326, remove = FALSE)
ARUlocations.sf <- merge(ARUlocations.sf,
                         data.frame(cell_unit = unit_effort$cell_unit,
                                    endDate = as.Date(unit_effort$endDate)),
                         by = "cell_unit")

bufferSizes <- c(50, 150, 300) # m
siteBuffers <- lapply(setNames(bufferSizes, bufferSizes),
                      function(b) st_buffer(st_transform(ARUlocations.sf, 3310), dist = b))
# California Albers (EPSG 3310), so buffer distances are in meters

mapview(siteBuffers[["300"]], col.regions = "red") +
  mapview(siteBuffers[["150"]]) +
  mapview(siteBuffers[["50"]], col.regions = "yellow")


################################################################################
# 5) Remove burned areas
################################################################################


# get fire perimeters ----------------------------------------------------------
# CAL FIRE, California Fire Perimeters (all), downloaded 2024
# https://gis.data.ca.gov/datasets/CALFIRE-Forestry::california-fire-perimeters-all/explore
fire <- read_sf(dsn = file.path(base_dir, "Fire perimeters/California_Fire_Perimeters_2024")) # 22,261 fires
# 1898-04-01 til 2023-12-22

# keep fires that started after the 2018 hyperspectral imagery (21-22 June 2018)
# 6 recent fires have no ALARM_DATE; for these use 21 June of YEAR_ (conservative)
fire <- fire %>%
  mutate(fireDate = as.Date(ALARM_DATE),
         fireDate = if_else(is.na(fireDate) & !is.na(YEAR_),
                            as.Date(paste0(YEAR_, "-06-21")), fireDate)) %>%
  filter(!is.na(fireDate), fireDate >= as.Date("2018-06-21")) %>%
  dplyr::select(fireDate)
fire <- st_make_valid(st_transform(fire, 3310)) # California Albers, meters
fire <- st_cast(st_zm(fire), "MULTIPOLYGON") # one geometry type, so mapview can draw it
# 2,091 fires

# burned sites for each buffer size --------------------------------------------
# a site burned if any part of its buffer burned between the imagery and the
# site's last full recording day. burnYear = year of the most recent such fire.

burnedSites <- lapply(siteBuffers, function(buf) {
  st_join(buf[c("cell_unit", "endDate")], fire, left = FALSE) %>%
    st_drop_geometry() %>%
    filter(fireDate <= endDate) %>%
    group_by(cell_unit) %>%
    dplyr::summarise(burnYear = max(as.numeric(format(fireDate, "%Y"))))
})
sapply(burnedSites, nrow) # number of burned sites for 50, 150 and 300 m buffers

# map only the fires that touch a 300 m buffer (faster to draw)
fireNear <- fire[lengths(st_intersects(fire, siteBuffers[["300"]])) > 0, ]
mapview(siteBuffers[["150"]]) + mapview(fireNear, col.regions = "red")

# three versions of siteDetections, one per buffer size ------------------------
# each version drops the sites that burned within its own buffer

siteDetectionsBuffers <- lapply(names(siteBuffers), function(b) {
  sitePoints <- merge(ARUlocations.sf["cell_unit"], burnedSites[[b]],
                      by = "cell_unit", all.x = TRUE)
  x <- merge(sitePoints, siteDetections, by = "cell_unit")
  # keep the original column order (cell_unit, long, lat, burnYear, species...)
  x <- x[c("cell_unit", "long", "lat", "burnYear",
           setdiff(names(x), c("cell_unit", "long", "lat", "burnYear", "geometry")))]
  subset(x, is.na(x$burnYear))
})
names(siteDetectionsBuffers) <- names(siteBuffers)
sapply(siteDetectionsBuffers, nrow) # sites remaining for 50, 150 and 300 m buffers
# 50m 659 sites
# 150m 644 sites
# 300m 627 sites

# save -------------------------------------------------------------------------

for (b in names(siteDetectionsBuffers)) {
  write_csv(siteDetectionsBuffers[[b]], paste0("siteDetections_buffer", b, "m_20261106.csv"))
}




################################################################################
# 6) Get 1 km habitat variables
################################################################################
# 10 Nov 2026. Equal-grain version: every habitat variable is at the BioCube
# resolution (30 arc-seconds, ~1 km). Foliar traits (mean and SD of the 30 m
# pixels in each 1 km cell) and canopy height were resampled to the BioCube grid
# in resampling.R and saved with the other layers in California_biocube_1km_layers.
# This replaces the earlier 30 m trait extraction (see git history).

# RUNNING LEFT OFF HERE ######

# get tahoe and yosemite boxes -------------------------------------------------

tahoe <- read_sf(file.path(base_dir, "tahoe box/tahoe_shp.shp"))
yosemite <- read_sf(file.path(base_dir, "yosemite box/polygon_yosemite_box.shp"))
mapview(tahoe) +
  mapview(yosemite)

# exclude ARUs not in boxes ----------------------------------------------------
mapview(siteDetectionsBuffers[["150"]]) + mapview(tahoe) + mapview(yosemite)
# not necessary


# buffers for the sites kept after fire removal --------------------------------
# one set per buffer size (sections 4-5), in WGS84 like the 1 km layers

bufferPolys <- lapply(names(siteBuffers), function(b)
  siteBuffers[[b]] %>%
    filter(cell_unit %in% siteDetectionsBuffers[[b]]$cell_unit) %>%
    st_transform(4326) %>%
    dplyr::select(cell_unit, long, lat))
names(bufferPolys) <- names(siteBuffers)


# 1 km layers ------------------------------------------------------------------

path_CA <- file.path(base_dir, "California_biocube_1km_layers")
CAbiocube_files <- list.files(path_CA, pattern = "\\.tif$")
length(CAbiocube_files) # 124 = 94 BioCube layers + 14 trait means + 14 trait SDs + 2 canopy height

# check one trait against the buffers
mapview(raster(file.path(path_CA, "CA_Physiology_Nitrogen_mean_1km.tif"))) +
  mapview(bufferPolys[["150"]])


# area-weighted mean of each layer within each buffer --------------------------
# Each layer is loaded once and extracted for all three buffer sizes.
# terra::extract with exact = TRUE weights each 1 km cell by the share of the
# buffer it covers (same values as raster::extract with weights = TRUE, faster).
# Results tables have no geometry; the extraction uses the polygons in bufferPolys.

vectPolys <- lapply(bufferPolys, terra::vect)
BioCube_varsBuffers <- lapply(bufferPolys, st_drop_geometry)

for (f in CAbiocube_files) {
  r <- terra::rast(file.path(path_CA, f))[[1]]
  for (b in names(vectPolys)) {
    BioCube_varsBuffers[[b]][[f]] <- terra::extract(r, vectPolys[[b]], fun = mean,
                                                    exact = TRUE, na.rm = TRUE)[, 2]
  }
  message(f, " done")
}
rm(r)

# save the raw extraction right away (original file names as columns).
# To skip the extraction on a rerun, load this instead of running the loop above:
# BioCube_varsBuffers <- readRDS("BioCube_varsBuffers_raw_20261110.rds")
write_rds(BioCube_varsBuffers, "BioCube_varsBuffers_raw_20261110.rds")


# tidy variable names ----------------------------------------------------------
# tidyNames.csv includes the resampled trait and canopy height layers (added by resampling.R)

tidy <- read_csv("tidyNames.csv")
rename_map <- setNames(tidy$VariableName, tidy$FileName)

BioCube_varsBuffers <- lapply(BioCube_varsBuffers, function(x) {
  hit <- names(x) %in% names(rename_map)
  names(x)[hit] <- rename_map[names(x)[hit]]
  x
})
sapply(BioCube_varsBuffers, function(x) sum(!names(x) %in% c("cell_unit", "long", "lat", rename_map))) # should be 0 (all layers renamed)


# merge with detections and remove sites with no foliar trait data -------------
# Nitrogen is NA where the 1 km cell(s) under a buffer have no valid 30 m trait pixels

siteDetections_foliarTraits_BioCubeBuffers <- lapply(names(bufferPolys), function(b) {
  x <- merge(siteDetectionsBuffers[[b]],
             BioCube_varsBuffers[[b]][setdiff(names(BioCube_varsBuffers[[b]]), c("long", "lat"))],
             by = "cell_unit")
  x[!is.na(x$Nitrogen), ]
})
names(siteDetections_foliarTraits_BioCubeBuffers) <- names(bufferPolys)
sapply(siteDetections_foliarTraits_BioCubeBuffers, nrow) # sites for 50, 150, 300 m (30 m traits: 584, 578, 570)
sapply(siteDetections_foliarTraits_BioCubeBuffers, function(x) sum(is.na(st_drop_geometry(x)))) # NA values


# save -------------------------------------------------------------------------

for (b in names(bufferPolys)) {
  write_rds(BioCube_varsBuffers[[b]], paste0("BioCube_vars_buffer", b, "m_20261110.rds"))
  write_rds(siteDetections_foliarTraits_BioCubeBuffers[[b]],
            paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261110.rds"))
  write_csv(st_drop_geometry(siteDetections_foliarTraits_BioCubeBuffers[[b]]),
            paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261110.csv"))
}


################################################################################
# 7) add USGS CONUS GAP 2011 landcover type
################################################################################
# landcover type at each ARU location (point value), added to each buffer version

# load data
YosemiteGAP <- raster(file.path(base_dir, "Yosemite_GAP_2011.tif"))
TahoeGAP <- raster(file.path(base_dir, "Tahoe_GAP_2011.tif"))
YosemiteGAP$landcover[YosemiteGAP$landcover==0] <- NA
TahoeGAP$landcover[TahoeGAP$landcover==0] <- NA
x <- list(YosemiteGAP, TahoeGAP)
GAP <- do.call(merge, x)
mapview(GAP)

# key names
coverNames <- read.csv(file.path(base_dir, "USGS_CONUS_GAP_2011_key.csv"))
coverNames <- coverNames[c("landcover_ID", "landcover_type")]

# extract and add key names
siteDetections_foliarTraits_BioCubeBuffers <- lapply(siteDetections_foliarTraits_BioCubeBuffers, function(x) {
  sites <- x[c("cell_unit", "long", "lat")]
  x <- as.data.frame(x)
  x$landcover_ID <- raster::extract(GAP$layer, sites, na.rm = TRUE)
  x <- merge(x, coverNames, by = "landcover_ID")
  x$landcover_ID <- NULL
  x
})
sapply(siteDetections_foliarTraits_BioCubeBuffers, nrow) # sites with a landcover type, 50, 150, 300 m
# 584 obs 50m
# 578 obs 150m
# 570 obs 300m

################################################################################
# 8) which points are too close to each other?
################################################################################
# 06 Nov 2026. Automated, run separately for each buffer version.
# When two sites are < 500 m apart, one is dropped: the site with more missing
# values across variables, or if tied, the site with fewer full recording days.
# Pairs are handled from closest to farthest, so a site is not dropped for a
# neighbour that has already been dropped.
# (Manual choice before, 150 m: dropped C5196_U1, C3824_U2, C6940_U2.)

drop_close_sites <- function(x, min_dist = 500) {
  pts <- st_as_sf(x[c("cell_unit", "long", "lat")], coords = c("long", "lat"), crs = 4326)
  d <- st_distance(pts)
  units(d) <- NULL
  idx <- which(d < min_dist & upper.tri(d), arr.ind = TRUE)
  if (nrow(idx) == 0) return(list(x = x, pairs = NULL))
  pairs <- data.frame(site1 = x$cell_unit[idx[, 1]], site2 = x$cell_unit[idx[, 2]],
                      dist_m = round(d[idx]))
  pairs <- pairs[order(pairs$dist_m), ]
  info <- data.frame(cell_unit = x$cell_unit,
                     nMissing = rowSums(is.na(x[setdiff(names(x), "geometry")])),
                     active_days = unit_effort$active_days[match(x$cell_unit, unit_effort$cell_unit)])
  pairs$dropped <- NA_character_
  dropped <- character(0)
  for (i in seq_len(nrow(pairs))) {
    a <- info[info$cell_unit == pairs$site1[i], ]
    b <- info[info$cell_unit == pairs$site2[i], ]
    if (a$cell_unit %in% dropped || b$cell_unit %in% dropped) next
    dropA <- if (a$nMissing != b$nMissing) a$nMissing > b$nMissing else a$active_days < b$active_days
    pairs$dropped[i] <- if (dropA) a$cell_unit else b$cell_unit
    dropped <- c(dropped, pairs$dropped[i])
  }
  list(x = x[!x$cell_unit %in% dropped, ], pairs = pairs)
}

closeSites <- lapply(siteDetections_foliarTraits_BioCubeBuffers, drop_close_sites)
lapply(closeSites, `[[`, "pairs") # pairs < 500 m apart and the site dropped from each
siteDetections_foliarTraits_BioCubeBuffers <- lapply(closeSites, `[[`, "x")
sapply(siteDetections_foliarTraits_BioCubeBuffers, nrow) # sites kept, 50, 150, 300 m
# 581 obs 50m
# 575 obs 150m
# 567 obs 300m

# how many unique cells? --------------------------------------------------------
sapply(siteDetections_foliarTraits_BioCubeBuffers,
       function(x) length(unique(str_extract(x$cell_unit, ".*(?=_)")))) # 285 before (150 m)
# 289 obs 50m
# 285 obs 150m
# 279 obs 300m

# average distance between ARUs in the same cell? ------------------------------
# California Albers (EPSG 3310) so distances are true meters
# (the earlier version used Web Mercator, EPSG 3857, which overstates
# distances at this latitude by about 27%)
avg_cell_distance <- function(x) {
  s <- st_as_sf(x[c("cell_unit", "long", "lat")], coords = c("long", "lat"), crs = 4326) %>%
    st_transform(3310) %>%
    mutate(cell = str_extract(cell_unit, ".*(?=_)"))
  d <- sapply(split(s, s$cell), function(g) {
    if (nrow(g) < 2) return(NA_real_)
    m <- as.numeric(st_distance(g))
    mean(m[lower.tri(matrix(NA, nrow(g), nrow(g)))])
  })
  mean(d, na.rm = TRUE)
}
sapply(siteDetections_foliarTraits_BioCubeBuffers, avg_cell_distance) # meters (1206 before, Web Mercator)
# 951m 50m
# 950m 150m
# 948m 300m

################################################################################
# 9) remove species with insufficient data
################################################################################
# 06 Nov 2026. Automated, run separately for each buffer version.
# Remove species detected at 30 or fewer sites or with fewer than 300 total detections
# (same cut-offs as the curve fitting in section 3)
# (raw counts on 7 h days, siteCounts from section 3, at the sites in each version).

count_species <- function(x) {
  spCols <- intersect(names(x), unique(siteCounts$common_name))
  dets <- siteCounts %>%
    filter(cell_unit %in% x$cell_unit, common_name %in% spCols) %>%
    group_by(common_name) %>%
    dplyr::summarise(totalDetections = sum(n_detections))
  out <- data.frame(common_name = spCols, totalSites = colSums(x[spCols] != 0))
  out <- merge(out, dets, by = "common_name", all.x = TRUE)
  out$totalDetections[is.na(out$totalDetections)] <- 0
  out$keep <- out$totalSites > 30 & out$totalDetections >= 300
  out[order(out$totalSites), ]
}

speciesCounts <- lapply(siteDetections_foliarTraits_BioCubeBuffers, count_species)
lapply(speciesCounts, function(s) s[!s$keep, ]) # species removed from each version
sapply(speciesCounts, function(s) sum(s$keep)) # species kept, 50, 150, 300 m
# all buffer sizes removed: California Thrasher, Tree Swallow, Song Sparrow, 
# Wild Turkey, Brown-headed Cowbird, Vaux's Swift, Yellow-breasted Chat, 
# Brewer's Sparrow, Lawrence's Goldfinch, Sooty Grouse
# all buffer sizes kept 92 species

siteDetections_foliarTraits_BioCubeBuffers <- lapply(names(speciesCounts), function(b) {
  drop <- speciesCounts[[b]]$common_name[!speciesCounts[[b]]$keep]
  siteDetections_foliarTraits_BioCubeBuffers[[b]][setdiff(names(siteDetections_foliarTraits_BioCubeBuffers[[b]]), drop)]
})
names(siteDetections_foliarTraits_BioCubeBuffers) <- names(speciesCounts)

# save -------------------------------------------------------------------------

for (b in names(siteDetections_foliarTraits_BioCubeBuffers)) {
  write_csv(speciesCounts[[b]], paste0("speciesCounts_buffer", b, "m_20261106.csv"))
  write_rds(siteDetections_foliarTraits_BioCubeBuffers[[b]],
            paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261106.rds"))
  x <- siteDetections_foliarTraits_BioCubeBuffers[[b]]
  x$geometry <- NULL
  write_csv(x, paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261106.csv"))
}

# 150 m version, used for the site map in section 10
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCubeBuffers[["150"]]
speciesSites <- speciesCounts[["150"]][c("common_name", "totalSites")]
speciesDetections <- speciesCounts[["150"]][c("common_name", "totalDetections")]


################################################################################
# 10) make a nice map
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

# plot it (Figure 2A)
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
# 11) check for colinearity
################################################################################
# 06 Nov 2026. Variables are selected by name (all foliar traits and all BioCube
# layers). Correlations are checked on the 150 m version, where the choices of
# which variable to keep were made; the same removal list is then applied to
# all three buffer versions so the PCA uses identical variables.

# A 0.95 threshold only loses us two foliar traits. A 0.9 threshold would lose more.

corThreshold <- 0.95

# all environmental variables
spatVarsAll <- setdiff(names(BioCube_varsBuffers[["150"]]), c("cell_unit", "long", "lat"))
length(spatVarsAll) # 124 = 94 BioCube layers + 14 trait means + 14 trait SDs + 2 canopy height
# (30 m version: 14 traits + 96 BioCube layers = 110 vars)

# pairs of variables with |r| above the threshold (each pair listed once)
find_high_cor <- function(x, vars, thr = corThreshold) {
  m <- cor(x[vars], use = "complete.obs")
  m[lower.tri(m, diag = TRUE)] <- NA
  as.data.frame(as.table(m)) %>%
    filter(!is.na(Freq), abs(Freq) > thr) %>%
    dplyr::rename(var1 = Var1, var2 = Var2, r = Freq) %>%
    arrange(desc(abs(r)))
}

# check correlations for the current set of variables -------------------------
# Variables are removed one at a time: add a name to removeVars (below), then rerun
# check_cor(removeVars). The data are not changed until the end of this section,
# so the plots and the pair list can be regenerated as often as needed.
# Variables already dropped from the table are skipped, so this also works after
# the removal step has been run.

check_cor <- function(removeVars = c(), thr = corThreshold, b = "150") {
  x <- siteDetections_foliarTraits_BioCubeBuffers[[b]]
  vars <- intersect(setdiff(spatVarsAll, removeVars), names(x))
  m <- cor(x[vars], use = "complete.obs")
  corrplot(m,
           method = "color",
           order = "hclust",
           hclust.method = "complete",
           addrect = 3,
           tl.cex = 0.8)
  h <- m
  h[abs(h) < thr | diag(ncol(h)) == 1] <- 0
  corrplot(h,
           method = "color",
           order = "hclust",
           hclust.method = "complete",
           tl.cex = 0.8)
  find_high_cor(x, vars, thr) # pairs above the threshold, highest first
}

check_cor() # all variables (150 m)

# in Cali, there are dry summers and wet winters, so some bioclim vars are duplicate.
# bio16 and bio19 are 0.999 identical. precip of wettest/coldest quarter. remove bio16(wettest)?
# bio09 and bio10 are 0.999 identical. temp of dryest/warmest quarter. remove bio09(dryest)?
# gdd5 and gdd10 are 0.999 identical. remove gdd5?
# bio08 and bio11 are 0.998 identical. temp of wettest/coldest quarter. remove bio11(coldest)?
# All these are > 95% correlated with DTM: "bio01", "bio05", "bio10", "bio08", "PETDryQuart", "bio06", "annualPET", "gdd10", "ngd5", "PETseas", "PETWarmQuart", "thermInd"
# gdd10 and ngd10, and gg0 removing ngd10 and gdd0
# Slope and TRI, remove TRI
# swe and scd 0.98 identical, remove swe
# bio12 and bio19 0.98, keep annualPrecip(12), remove 19
# LMA and starch 0.98, remove starch (goodbye one foliar trait)
# ngdo, scd 0.97, remove ngd0

# remove (same list for all three buffer versions)
removeVars <- c("bio16", "bio09", "gdd5", "bio11", "bio01", "bio05", "ngd10", 
                "gdd0", "bio10", "bio08", "PETDryQuart", "bio06", "annualPET", 
                "TRI", "swe", "gdd10", "ngd5", "PETseas", "PETWarmQuart", 
                "thermInd", "bio19", "Starch", "ngd0", "Sulfur", "DTM_std", "bio18")
check_cor(removeVars) # rerun after each addition to removeVars
removeVars[!removeVars %in% spatVarsAll] # should be empty (catches typos in names)

siteDetections_foliarTraits_BioCubeBuffers <- lapply(siteDetections_foliarTraits_BioCubeBuffers,
                                                     function(x) x[setdiff(names(x), removeVars)])
spatVars <- setdiff(spatVarsAll, removeVars)
length(spatVars) # variables used in the PCA
# length 87 # now 84

# check: no remaining pairs above the threshold in any buffer version
lapply(siteDetections_foliarTraits_BioCubeBuffers, find_high_cor, vars = spatVars) # should all be empty
# 300m has 3 more vars above the 0.95 threshold
# LMA and Sulfur -0.952, remove sulfur
# DTM_std and Slope, 0.952, remove DTM_std
# bio17 and bio18 0.951, remove bio18

# save updated version
for (b in names(siteDetections_foliarTraits_BioCubeBuffers)) {
  write_rds(siteDetections_foliarTraits_BioCubeBuffers[[b]],
            paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261106.rds"))
  x <- siteDetections_foliarTraits_BioCubeBuffers[[b]]
  x$geometry <- NULL
  write_csv(x, paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261106.csv"))
}

# TODO 06 Nov 2026: sections 12 onward still use a single table. Until they are
# converted, they use the 150 m version.
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCubeBuffers[["150"]]


# check for NA values ----------------------------------------------------------
lapply(siteDetections_foliarTraits_BioCubeBuffers, function(x) {
  n <- colSums(is.na(x[spatVars])); n[n > 0]
})
# no NAs. good to go. 

################################################################################
# 12) Niche Space PCA
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
# 13) PCA loadings (supp figures)
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
# 14) render landcover labels
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
# 15) test plot for errored species polygons
################################################################################

# some species weren't plotting correctly because not enough points
# Solution #1) 1 = 1 detection per 10 days (x10 increase) (NOPE)
# Solution #2) add jitter to all points (non-identical coordinates) (this worked)

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
# 16) Species Niche Polygons
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
# (KDE fails for Downy.Woodpecker and Lawrence.s.Goldfinch)

speciesPolygons <- map_dfr(species, ~ calculate_polygons(.x), .id = "species1")
speciesPolygons$species1 <- NULL
speciesPolygons$level <-NULL

# save it ----------------------------------------------------------------------

setwd(file.path(base_dir, "Draft 1"))
setwd(file.path(base_dir, "AviSpace_GitHub/AviSpace4.0"))

write_rds(speciesPolygons, "speciesPolygons_jittered_detPerDay_20250624.rds")


################################################################################
# 17) plot 2 species niche overlap
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


# make a not warped version (aspect ratio should match PC%)---------------------

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
# 18) Calculate all pairwise overlaps
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

setwd(file.path(base_dir, "Draft 1"))
setwd(file.path(base_dir, "AviSpace_GitHub/AviSpace4.0"))

write.csv(pairwise_combinations, "speciesNicheOverlapAreas_jittered_detPerWeek_20250624.csv")



################################################################################
# 19) Species niche Areas & Percent overlap
################################################################################

# species areas
speciesAreas <- data.frame(species = speciesPolygons$species,
                           area = st_area(speciesPolygons))

# save it 
write_rds(speciesAreas, "speciesNicheAreas_20250624.rds")
#speciesAreas <- readRDS("speciesNicheAreas_20250624.rds")

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
# 20) calculate Range overlaps
################################################################################


# toy - single species pair ----------------------------------------------------

# read in data
x <- st_read(file.path(base_dir, "RangeMaps/Yellow-rumped Warbler/yerwar_range_2022.gpkg"))
mapview(x)
y <- st_read(file.path(base_dir, "RangeMaps/Red-breasted Nuthatch/rebnut_range_2022.gpkg"))
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
setwd(file.path(base_dir, "RangeMaps"))

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
setwd(file.path(base_dir, "Draft 1"))

# write csv

setwd(file.path(base_dir, "AviSpace_GitHub/AviSpace4.0"))

write.csv(rangeOverlaps, "rangeOverlaps_20250624.csv")
write_rds(rangeOverlaps, "rangeOverlaps_20250624.rds")


# make dot version -------------------------------------------------------------
dotVclean_nameMatches <-  read.csv(file.path(base_dir, "dotVclean_nameMatches.csv"))

rangeOverlapsDot <- merge(rangeOverlaps, dotVclean_nameMatches, by.x = "species1", by.y = "cleanName")
rangeOverlapsDot$species1 <- NULL
names(rangeOverlapsDot)[names(rangeOverlapsDot) == "dotName"] <- "species1"

rangeOverlapsDot <- merge(rangeOverlapsDot, dotVclean_nameMatches, by.x = "species2", by.y = "cleanName")
rangeOverlapsDot$species2 <- NULL
names(rangeOverlapsDot)[names(rangeOverlapsDot) == "dotName"] <- "species2"

# save dot version -------------------------------------------------------------

write.csv(rangeOverlapsDot, "rangeOverlaps_dotNames_20250624.csv")
write_rds(rangeOverlapsDot, "rangeOverlaps_dotNames_20250624.rds")
rangeOverlaps <- readRDS("rangeOverlaps_dotNames_20250623.rds")

################################################################################
# 21) Range Overlap Maps
################################################################################


# read in data
x <- st_read(file.path(base_dir, "RangeMaps/Yellow-rumped Warbler/yerwar_range_2022.gpkg"))
mapview(x)
y <- st_read(file.path(base_dir, "RangeMaps/Red-breasted Nuthatch/rebnut_range_2022.gpkg"))
mapview(y)
z <- st_read(file.path(base_dir, "RangeMaps/Purple Finch/purfin_range_2022.gpkg"))
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
map_paths <- list.files(file.path(base_dir, "RangeMaps/"), 
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
# 22) Combine range overlaps and niche overlaps
################################################################################

# ensure column names are informative ------------------------------------------
NicheOverlap <- percentOverlap[, c("species1", "species2", "percentOverlap")]
names(NicheOverlap)[names(NicheOverlap) == 'percentOverlap'] <- 'niche_overlap_percent'

rangeOverlaps <- rangeOverlaps[, c("species1", "species2", "overlap_percent")]
names(rangeOverlaps)[names(rangeOverlaps) == "overlap_percent"] <- "range_overlap_percent" 


# match species name formats ---------------------------------------------------
dotVclean_nameMatches <-  read.csv(file.path(base_dir, "dotVclean_nameMatches.csv"))

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
# 23) plot range vs niche overlap
################################################################################

# mark congeneric pairs --------------------------------------------------------

# load scientific names
SciNames <- read.csv(file.path(base_dir, "species_scientificNames.csv"))
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
#write.csv(NicheRangeOverlap, "NicheRangeOverlap_cleanNames_20250624.csv")
#write_rds(NicheRangeOverlap, "NicheRangeOverlap_cleanNames_20250624.rds")
NicheRangeOverlap <- readRDS(file.path(base_dir, "NicheRangeOverlap_cleanNames_20250624.rds"))

# make dot names version -------------------------------------------------------

NicheRangeOverlapDot <- merge(NicheRangeOverlap, dotVclean_nameMatches, by.x = "species1", by.y = "cleanName")
NicheRangeOverlapDot$species1 <- NULL
names(NicheRangeOverlapDot)[names(NicheRangeOverlapDot) == "dotName"] <- "species1"

NicheRangeOverlapDot <- merge(NicheRangeOverlapDot, dotVclean_nameMatches, by.x = "species2", by.y = "cleanName")
NicheRangeOverlapDot$species2 <- NULL
names(NicheRangeOverlapDot)[names(NicheRangeOverlapDot) == "dotName"] <- "species2"

# save dot version -------------------------------------------------------------

#write.csv(NicheRangeOverlapDot, "NicheRangeOverlap_dotNames_20250624.csv")
#write_rds(NicheRangeOverlapDot, "NicheRangeOverlap_dotNames_20250624.rds")
NicheRangeOverlapDot <- readRDS("NicheRangeOverlap_dotNames_20250624.rds")

# remove potentially cross-identified sister species ---------------------------

# remove sapsucker sister species
NicheRangeOverlapDot <- NicheRangeOverlapDot[!((NicheRangeOverlapDot$species1 == "Red.breasted.Sapsucker" & NicheRangeOverlapDot$species2 == "Red.naped.Sapsucker") | 
                                                 (NicheRangeOverlapDot$species1 == "Red.naped.Sapsucker" & NicheRangeOverlapDot$species2 == "Red.breasted.Sapsucker")), ]
# remove vireo sister species
NicheRangeOverlapDot <- NicheRangeOverlapDot[!((NicheRangeOverlapDot$species1 == "Cassin.s.Vireo" & NicheRangeOverlapDot$species2 == "Plumbeous.Vireo") | 
                                                 (NicheRangeOverlapDot$species1 == "Plumbeous.Vireo" & NicheRangeOverlapDot$species2 == "Cassin.s.Vireo")), ]
# plot the updated data (linear model)
ggplot(NicheRangeOverlapDot, aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_smooth(aes(fill = congeneric), method = "lm", alpha = 0.2) +  
  theme_minimal() +
  scale_color_manual(values = c("purple", "darkgreen")) +
  scale_fill_manual(values = c("purple", "darkgreen")) +  
  geom_point(data=subset(NicheRangeOverlapDot, NicheRangeOverlapDot$sisterSpecies==TRUE), color="black") +
  theme(legend.title = element_blank()) +
  stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x, parse = TRUE) + # ggpmisc
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)")
ggsave("NicheRangeOverlap_light_20260706.pdf", height=5, width=7)

# make supp figure with points (linear model)
ggplot(NicheRangeOverlapDot, aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_point(alpha=0.2) +
  scale_color_manual(values = c("purple", "darkgreen")) +
  scale_fill_manual(values = c("purple", "darkgreen")) +
  facet_grid(~congeneric) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x, parse = TRUE) +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)") +
  geom_smooth(aes(fill = congeneric), method = "lm", alpha = 0.2) +
  theme_minimal()
ggsave("NicheRangeOverlap_points_20260706.pdf", height=5, width=7)


################################################################################
# 24) Mixed effects model
################################################################################

# Mixed effects model - accounts for species non-independence (the same species appear in multiple pairs)

# --- Prep ---
# make sure mixed membership includes the same full set of species 
all_species <- sort(unique(c(as.character(NicheRangeOverlapDot$species1),
                             as.character(NicheRangeOverlapDot$species2))))
NicheRangeOverlapDot$species1   <- factor(NicheRangeOverlapDot$species1, levels = all_species)
NicheRangeOverlapDot$species2   <- factor(NicheRangeOverlapDot$species2, levels = all_species)
NicheRangeOverlapDot$congeneric <- factor(NicheRangeOverlapDot$congeneric)

# --- Prior: one multiple-membership random effect (G) + residual (R) ---
# Parameter-expanded prior on the random effect helps mixing.
prior <- list(
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))
)

# --- Model: lm fixed structure + species identity as an MM random effect ---
chains <- lapply(1:4, function(i) {
  set.seed(i)
  MCMCglmm(range_overlap_percent ~ niche_overlap_percent * congeneric,
           random = ~ mm(species1 + species2),
           family = "gaussian", data = NicheRangeOverlapDot, prior = prior,
           nitt = 260000, burnin = 60000, thin = 200, verbose = FALSE)
})
gelman.diag(mcmc.list(lapply(chains, function(m) m$Sol)))  

write_rds(chains, "MixedEffectsModel_4chains_20260713.RDS")

# --- Bundle the 4 chains ---
sol_list <- mcmc.list(lapply(chains, function(m) m$Sol))   # fixed effects
vcv_list <- mcmc.list(lapply(chains, function(m) m$VCV))   # variance components

# --- Convergence checks across all 4 chains ---
gelman.diag(sol_list)      # R-hat: point est + upper CI (want ~1.00)
gelman.plot(sol_list)      # shows R-hat settling to 1 as iterations accrue
plot(sol_list)             # trace + density; want all 4 chains overlapping in one band
effectiveSize(sol_list)    # pooled ESS across chains (~4,000 if each was ~1,000)
autocorr.diag(sol_list)    # autocorrelation should drop off fast

# --- Pool chains for the reported posterior (inference) ---
sol_all <- do.call(rbind, lapply(chains, function(m) m$Sol))
vcv_all <- do.call(rbind, lapply(chains, function(m) m$VCV))
apply(sol_all, 2, mean)                 # posterior means
t(apply(sol_all, 2, quantile, c(.025, .975)))   # 95% credible intervals


# ---- Build fitted lines from the POOLED posterior (all 4 chains) ----
sol  <- sol_all                                          # pooled draws, ~4,000 samples
lev  <- levels(NicheRangeOverlapDot$congeneric)          # lev[1] = reference (heterogeneric)
xseq <- seq(min(NicheRangeOverlapDot$niche_overlap_percent),
            max(NicheRangeOverlapDot$niche_overlap_percent), length.out = 100)

# pMCMC computed directly from the pooled draws (same definition MCMCglmm uses)
p_from_draws <- function(x) 2 * max(0.5 / length(x), min(mean(x > 0), mean(x < 0)))

dslope <- sol[, "niche_overlap_percent:congenericSame Genus"]
pmc    <- p_from_draws(dslope)
eff_label <- sprintf("Difference in slopes = %.2f  [%.2f, %.2f],  pMCMC = %.3f",
                     mean(dslope), quantile(dslope, .025), quantile(dslope, .975), pmc)

# per-group posterior draws of intercept (b0) and slope (b1)
coefs <- function(g) {
  if (g == lev[1]) {
    list(b0 = sol[, "(Intercept)"],
         b1 = sol[, "niche_overlap_percent"])
  } else {
    list(b0 = sol[, "(Intercept)"] + sol[, paste0("congeneric", g)],
         b1 = sol[, "niche_overlap_percent"] +
           sol[, paste0("niche_overlap_percent:congeneric", g)])
  }
}

pred_df <- do.call(rbind, lapply(lev, function(g) {
  cf    <- coefs(g)
  preds <- sapply(xseq, function(xi) cf$b0 + cf$b1 * xi)   # draws x xgrid
  fit   <- colMeans(preds)
  se    <- apply(preds, 2, sd)
  data.frame(niche_overlap_percent = xseq, congeneric = g,
             fit = fit, lower = fit - 1.96 * se, upper = fit + 1.96 * se)
}))
pred_df$congeneric <- factor(pred_df$congeneric, levels = lev)

# slope labels (posterior mean + 95% credible interval)
lab_df <- do.call(rbind, lapply(lev, function(g) {
  b1 <- coefs(g)$b1
  data.frame(congeneric = factor(g, levels = lev),
             label = sprintf("slope = %.2f [%.2f, %.2f]",
                             mean(b1), quantile(b1, 0.025), quantile(b1, 0.975)))
}))
lab_df$x <- 0  # left side
lab_df$y <- c(0.43, 0.40)[seq_along(lev)]   # stacked
eff_df <- data.frame(x = min(xseq), y = 0.37, label = eff_label)

# ---- Plot it ----
# Figure 8
ggplot(NicheRangeOverlapDot,
       aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_ribbon(data = pred_df,
              aes(y = fit, ymin = lower, ymax = upper, fill = congeneric),
              alpha = 0.2, color = NA) +
  geom_line(data = pred_df, aes(y = fit), linewidth = 1) +
  #geom_point(data = subset(NicheRangeOverlapDot, sisterSpecies == TRUE), color = "black") +
  geom_text(data = lab_df, aes(x = x, y = y, label = label),
            hjust = 0, show.legend = FALSE) +
  geom_text(data = eff_df, aes(x = x, y = y, label = label),
            hjust = 0, color = "black", inherit.aes = FALSE) +
  scale_color_manual(values = c("purple", "darkgreen")) +
  scale_fill_manual(values = c("purple", "darkgreen")) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)")
ggsave("NicheRangeOverlap_points_20260713.pdf", height=5, width=7)

# supp figure 
ggplot(NicheRangeOverlapDot,
       aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_ribbon(data = pred_df,
              aes(y = fit, ymin = lower, ymax = upper, fill = congeneric),
              alpha = 0.2, color = NA) +
  geom_point(alpha = 0.2) +
  geom_line(data = pred_df, aes(y = fit), linewidth = 1) +
  geom_text(data = lab_df, aes(x = x, y = y+0.3, label = label),
            hjust = 0, show.legend = FALSE) +
  scale_color_manual(values = c("purple", "darkgreen")) +
  scale_fill_manual(values = c("purple", "darkgreen")) +
  facet_grid(~ congeneric) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)")
ggsave("NicheRangeOverlap_points_20260706.pdf", height=5, width=7)

ggplot(subset(NicheRangeOverlapDot, NicheRangeOverlapDot$congeneric=="Same Genus"),
       aes(x = niche_overlap_percent, y = range_overlap_percent, color = sisterSpecies)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("darkgreen", "black")) +
  facet_grid(~ congeneric) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Niche Overlap (%)") +
  ylab("Range Overlap (%)")

################################################################################
# 25) permutation test
################################################################################

# ---- 1. Species -> genus map, straight from the genus columns ----
sp_gen <- unique(rbind(
  data.frame(species = as.character(NicheRangeOverlapDot$species1),
             genus   = as.character(NicheRangeOverlapDot$genus_species1)),
  data.frame(species = as.character(NicheRangeOverlapDot$species2),
             genus   = as.character(NicheRangeOverlapDot$genus_species2))
))
genus_of <- setNames(sp_gen$genus, sp_gen$species)

# ---- 2. Statistic: interaction (slope difference) from a given species->genus map ----
interaction_coef <- function(gmap) {
  cong_status <- gmap[as.character(NicheRangeOverlapDot$species1)] ==
    gmap[as.character(NicheRangeOverlapDot$species2)]
  fit <- lm(range_overlap_percent ~ niche_overlap_percent * cong_status,
            data = NicheRangeOverlapDot)
  unname(coef(fit)["niche_overlap_percent:cong_statusTRUE"])
}

obs <- interaction_coef(genus_of)                     # observed slope difference

# ---- 3. Null: shuffle genus labels across species (preserves genus-size distribution) ----
set.seed(1)
nperm   <- 5000
species <- names(genus_of)
null <- replicate(nperm, {
  permuted <- setNames(sample(unname(genus_of)), species)
  interaction_coef(permuted)
})

# ---- 4. Two-sided permutation p-value ----
p_perm <- (sum(abs(null) >= abs(obs)) + 1) / (nperm + 1)
c(observed_slope_diff = obs, p_perm = p_perm)

# optional: see where the observed value falls
hist(null, breaks = 40, main = "Permutation null: slope difference",
     xlab = "interaction coefficient"); abline(v = obs, col = "red", lwd = 2)



################################################################################
# 26) Figure 4 niche areas 
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
# 27) Figure 5 - matrix
################################################################################

# Heirarchical Matrix

# bring back sister species ----------------------------------------------------
NicheRangeOverlapDot <- readRDS("NicheRangeOverlap_dotNames_20250624.rds")
speciesAreas <- readRDS("speciesNicheAreas_20250624.rds")

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
hc <- hclust(as.dist(1 - overlap_matrix), method = "complete")
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
ggsave("Figure5_matrix.PDF", height=15, width=15)


# a reviewer asked for a simplified version with fewer species -----------------

# subset to fewer species ------------------------------------------------------
spSubset <- c("Pacific.Wren", "Brown.Creeper", "Red.breasted.Nuthatch", "Red.breasted.Sapsucker", "Red.naped.Sapsucker", "White.headed.Woodpecker",
              "Townsed.s.Solitaire", "Yellow.rumped.Warbler", "Hammond.s.Flycatcher", "Hermit.Warbler", "Hermit.Thrush", "Golden.crowned.Kinglet",
              "Mountain.Chickadee", "Olive.sided.Flycatcher", "Cassin.s.Finch", "Red.Crossbill", "Green.Tailed.Towhee", "Dusky.Flycatcher",
              "California.Quail", "California.Towhee", "Oak.Titmouse", "Blue.gray.Gnatcatcher", "Bewick.s.Wren", "Wrentit", "Nuttall.s.Woodpecker",
              "Mourning.Dove", "California.Scrub.Jay", "Woodhouse.s.Scrub.Jay", "Black.Phoebe", "Acorn.Woodpecker", "Bushtit",
              "Black.headed.Grosbeak", "Plumbeous.Vireo", "Purple.Finch", "Cassin.s.Vireo", "Nashville.Warbler", "Common.Raven", "American.Robin",
              "Hairy.Woodpecker", "Northern.Flicker", "Steller.s.Jay", "Dark.eyed.Junco", "Western.Tanager", "Western.Wood.Pewee")

NicheRangeOverlapDot <- NicheRangeOverlapDot[
  NicheRangeOverlapDot$species1 %in% spSubset &
    NicheRangeOverlapDot$species2 %in% spSubset,
]

speciesAreas <- speciesAreas[
  speciesAreas$species %in% spSubset,
]
# then rerun earlier script 

combined_plot <- plot_grid(
  heatmap_plot,
  bar_plot,
  align = "h",
  axis = "tblr",
  rel_widths = c(9, 2)
)

combined_plot

ggsave("Figure5_matrix_sup.PDF", height=8, width=9)

################################################################################
# 28) Species details supp table
################################################################################

# niche area
# total # detections
# total # sites
# 06 July 2026 - Connor has added birdnet validation data to this table. do not overwrite.

setwd(file.path(base_dir, "Draft 2_PNAS"))


speciesDetails <- merge(speciesAreas, speciesDetections)
colnames(speciesDetails)[colnames(speciesDetails) == 'area'] <- 'nicheArea'
speciesDetails <- merge(speciesDetails, speciesSites)
write_rds(speciesDetails, "speciesDetails_20250708.rds")
write_csv(speciesDetails, "speciesDetails_20250708.csv")





