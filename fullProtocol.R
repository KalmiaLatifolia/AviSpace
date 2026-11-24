
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
library(ks)      # weighted kernel density (section 15)

# Table of contents ------------------------------------------------------------
# 1) Processing ARU bioacoustic data - 72
# 2) Species complexes and cross-identification - 146
# 3) Adjust detection rates for day-of-year and year - 259
# 4) Site buffers - 497
# 5) Remove burned areas - 524
# 6) Get 1 km habitat variables - 591
# 7) add USGS CONUS GAP 2011 landcover type - 741
# 8) which points are too close to each other? - 773
# 9) remove species with insufficient data - 843
# 10) make a nice map - 891
# 11) check for colinearity - 934
# 12) Niche Space PCA - 1045
# 13) PCA loadings (supp figures) - 1086
# 14) render landcover labels - 1165
# 15) Species acoustic-use space: weighted 4-D kernel density - 1254
# 16) plot 2 species acoustic-use space overlap (test plot) - 1407
# 17) Species niche areas and pairwise overlaps: PC1 x PC2, PC3 x PC4 and 4-D - 1502
# 18) Range areas and pairwise range overlaps - 1571
# 19) Range overlap maps (two species) - 1769
# 20) Combine range overlaps and niche overlaps - 1871
# 21) Mixed effects models - 1967
# 22) Leave-one-genus-out models (main dataset) - 2182
# 23) permutation test - 2261
# 24) Figure 4 acoustic-use space volumes - 2320
# 25) Figure 5 - matrix - 2345
# 26) Species details supp table - 2494


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

# 23 Nov 2026: scrub-jays added. Woodhouse's Scrub-Jay is detected at 321 sites,
# none inside its eBird breeding range, at about 3-6% of the California
# Scrub-Jay rate away from that range. Note this check only finds segments
# labelled as both species; a call labelled only as the wrong species is not
# counted, so it can understate one-way confusion.
# 24 Nov 2026: all congeneric species are now checked, not just the vireos,
# sapsuckers and scrub-jays: every detected species is grouped by genus (first
# word of the scientific name), and genera with 2+ detected species are kept.
# BirdNET's non-bird labels (e.g. "Human vocal") have the same scientific and
# common name and are left out.
detectedCodes <- unique(det$species_code)
crossIDgroups <- commonNames %>%
  filter(species_code %in% detectedCodes, scientific_name != common_name) %>%
  transmute(member = scientific_name, group = word(scientific_name, 1), species_code) %>%
  group_by(group) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  as.data.frame()
crossIDgroups %>% dplyr::count(group) # genera checked and number of detected species in each
stopifnot(!anyNA(crossIDgroups$species_code))

memberDet <- det %>%             # slow
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
  filter(n_members > 1) %>% # only shared segments can form pairs (keeps the join small)
  dplyr::select(group, cell_unit, filename, relative_sec, species_code) %>%
  inner_join(., ., by = c("group", "cell_unit", "filename", "relative_sec"),
             suffix = c("_1", "_2"), relationship = "many-to-many") %>%
  filter(species_code_1 < species_code_2) %>%
  dplyr::count(group, species_code_1, species_code_2, name = "shared_segments") %>%
  mutate(species_1 = commonNames$common_name[match(species_code_1, commonNames$species_code)],
         species_2 = commonNames$common_name[match(species_code_2, commonNames$species_code)])
crossIDpairs

write_csv(crossID, "crossIdentification_species_20261124.csv")
write_csv(crossIDpairs, "crossIdentification_pairs_20261124.csv")
rm(memberDet)

# species complexes ------------------------------------------------------------
# Based on the check above, Red-breasted and Red-naped Sapsucker share a
# substantial share of segments (~25% and ~35% of their detections).
# Williamson's Sapsucker (0.5%) and Cassin's/Plumbeous Vireo (0.1% / 2.6%)
# are kept separate.
# California Scrub-Jay is fine, but Woodhouse's is probbaly mostly mis-IDs 
# and should be dropped.
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

# two pages with wider panels so labels are legible (20 Nov 2026)
phenoSpecies <- sort(unique(phenoCurves$common_name))
phenoPages <- split(phenoSpecies, ceiling(seq_along(phenoSpecies) / ceiling(length(phenoSpecies) / 2)))
pdf("FigureS_phenologyCurves_20261120.pdf", height = 8.5, width = 7.5)
for (sp in phenoPages) {
  print(ggplot(filter(phenoCurves, common_name %in% sp),
               aes(x = as.Date(doy - 1, origin = "2021-01-01"), y = fit,
                   color = year, fill = year)) +
          geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
          geom_line(linewidth = 0.4) +
          facet_wrap(~ common_name, scales = "free_y", ncol = 6,
                     labeller = as_labeller(function(x) stringr::str_wrap(gsub("/", "/ ", x), 16))) +
          scale_x_date(date_labels = "%b", date_breaks = "1 month") +
          labs(x = "Date", y = "Predicted detections per day (average site)",
               color = "Year", fill = "Year") +
          theme_minimal(base_size = 8) +
          theme(strip.text = element_text(size = 6.5),
                legend.position = "bottom"))
}
dev.off()

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


# snow and frost: no-data = 0 ---------------------------------------------------
# 17 Nov 2026. In the CHELSA-BIOCLIM+ snow cover days (scd), snow water
# equivalent (swe) and frost change frequency (fcf) layers, the no-data cells are
# the warm, low-elevation areas (no snow / no frost-change events), and the files
# contain no zeros. This looks like a mistake - NA is therefore treated as 0. These three layers are
# re-extracted with no-data set to 0 before averaging, so a buffer that is partly
# in a no-snow cell gets a correctly weighted mean. (Runs whether the loop above
# was run or the raw file was loaded.)

zeroFiles <- c("CA_Climate_CHELSA_scd_v2.tif", "CA_Climate_CHELSA_swe_v2.tif",
               "CA_Climate_CHELSA_fcf_v2.tif")
for (f in zeroFiles) {
  r <- terra::rast(file.path(path_CA, f))[[1]]
  r <- terra::ifel(is.na(r), 0, r)
  for (b in names(vectPolys)) {
    BioCube_varsBuffers[[b]][[f]] <- terra::extract(r, vectPolys[[b]], fun = mean,
                                                    exact = TRUE, na.rm = TRUE)[, 2]
  }
}
rm(r)


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


# exclude gsp ------------------------------------------------------------------
# 10 Nov 2026. Growing season precipitation (CHELSA-BIOCLIM+ gsp) is excluded: the
# source file has unflagged no-data cells (4294967295), giving extreme outliers
# (e.g. 3.3 million at C4833_U2, 150 m buffer).

BioCube_varsBuffers <- lapply(BioCube_varsBuffers, function(x) x[setdiff(names(x), "gsp")])

# exclude winter habitat indices and phenology amplitude -------------------------
# 17 Nov 2026. WHI_pFWOS, WHI_SCV, WHI_SSL (no-data at ~15 sites; these layers
# have true zeros, so no-data is not 0) and AMP (no-data at ~8 sites). Dropping
# these variables keeps those sites in the analysis.

BioCube_varsBuffers <- lapply(BioCube_varsBuffers, function(x)
  x[setdiff(names(x), c("WHI_pFWOS", "WHI_SCV", "WHI_SSL", "AMP"))])

# BioCube_varsBuffers has 122 vars
# merge with detections and remove sites with no foliar trait data -------------
# Nitrogen is NA where the 1 km cell(s) under a buffer have no valid 30 m trait pixels

siteDetections_foliarTraits_BioCubeBuffers <- lapply(names(bufferPolys), function(b) {
  x <- merge(siteDetectionsBuffers[[b]],
             BioCube_varsBuffers[[b]][setdiff(names(BioCube_varsBuffers[[b]]), c("long", "lat"))],
             by = "cell_unit")
  x <- x[!is.na(x$Nitrogen), ]
  # remove any site with a remaining missing habitat value (1 site with missing
  # SoilGrids bdod, clay, nitrogen and sand in the 50 and 150 m versions)
  envVars <- setdiff(names(BioCube_varsBuffers[[b]]), c("cell_unit", "long", "lat"))
  x[complete.cases(st_drop_geometry(x)[envVars]), ]
})
names(siteDetections_foliarTraits_BioCubeBuffers) <- names(bufferPolys)
sapply(siteDetections_foliarTraits_BioCubeBuffers, nrow) # sites for 50, 150, 300 m (30 m traits were: 584, 578, 570)(now 599, 592, 580)
sapply(names(siteDetections_foliarTraits_BioCubeBuffers), function(b)       # NA habitat values, should be 0. Yup. 
  sum(is.na(st_drop_geometry(siteDetections_foliarTraits_BioCubeBuffers[[b]])[
    setdiff(names(BioCube_varsBuffers[[b]]), c("cell_unit", "long", "lat"))])))

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
# 599 obs 50m
# 592 obs 150m
# 580 obs 300m

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
# 595 obs 50m
# 588 obs 150m
# 576 obs 300m

# how many unique cells? --------------------------------------------------------
sapply(siteDetections_foliarTraits_BioCubeBuffers,
       function(x) length(unique(str_extract(x$cell_unit, ".*(?=_)")))) # 285 before (150 m)
# 293 obs 50m
# 291 obs 150m
# 282 obs 300m

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
# 949m 150m
# 945m 300m

################################################################################
# 9) remove species with insufficient data
################################################################################
# 19 Nov 2026. Automated, run separately for each buffer version.
# Keep species detected at more than 30 sites and more than 300 total detections
# (as the curve fitting in section 3, except that section 3 fits species with
# 300 or more detections; a species with exactly 300 gets a curve there but is removed here)
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
  out$keep <- out$totalSites > 30 & out$totalDetections > 300
  out[order(out$totalSites), ]
}

speciesCounts <- lapply(siteDetections_foliarTraits_BioCubeBuffers, count_species)
lapply(speciesCounts, function(s) s[!s$keep, ]) # species removed from each version
sapply(speciesCounts, function(s) sum(s$keep)) # species kept, 50, 150, 300 m
# all buffer sizes removed: California Thrasher, Tree Swallow, Song Sparrow, 
# Wild Turkey, Brown-headed Cowbird, Vaux's Swift, Yellow-breasted Chat, 
# Brewer's Sparrow, Lawrence's Goldfinch, Sooty Grouse
# all buffer sizes kept 92 species

# drop those species from the main dataset
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

################################################################################
# 10) make a nice map
################################################################################

# 150 m version, used for the site map in section 10
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCubeBuffers[["150"]]
speciesSites <- speciesCounts[["150"]][c("common_name", "totalSites")]
speciesDetections <- speciesCounts[["150"]][c("common_name", "totalDetections")]

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
length(spatVarsAll) # 119 = 89 BioCube layers (94 minus gsp, WHI_pFWOS, WHI_SCV, WHI_SSL, AMP) + 14 trait means + 14 trait SDs + 2 canopy height
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
                "thermInd", "bio19", "Starch", "ngd0", "Sulfur", "DTM_std", 
                "bio18", "Nitrogen", "Sugar_std")
check_cor(removeVars) # rerun after each addition to removeVars
removeVars[!removeVars %in% spatVarsAll] # should be empty (catches typos in names)

siteDetections_foliarTraits_BioCubeBuffers <- lapply(siteDetections_foliarTraits_BioCubeBuffers,
                                                     function(x) x[setdiff(names(x), removeVars)])
spatVars <- setdiff(spatVarsAll, removeVars)
length(spatVars) # variables used in the PCA
# length 87 # now 91

# check: no remaining pairs above the threshold in any buffer version
lapply(siteDetections_foliarTraits_BioCubeBuffers, find_high_cor, vars = spatVars) # should all be empty
# 300m has 3 more vars above the 0.95 threshold
# LMA and Sulfur -0.952, remove sulfur
# DTM_std and Slope, 0.952, remove DTM_std
# bio17 and bio18 0.951, remove bio18
# 300m Sugar_std and NSC_std, 0.9508

# save updated version
for (b in names(siteDetections_foliarTraits_BioCubeBuffers)) {
  write_rds(siteDetections_foliarTraits_BioCubeBuffers[[b]],
            paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261119.rds"))
  x <- siteDetections_foliarTraits_BioCubeBuffers[[b]]
  x$geometry <- NULL
  write_csv(x, paste0("siteDetections_foliarTraits_BioCube_buffer", b, "m_20261119.csv"))
}

# 150 m version (used for the landcover labels in the PCA plots, section 14)
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCubeBuffers[["150"]]


# check for NA values ----------------------------------------------------------
lapply(siteDetections_foliarTraits_BioCubeBuffers, function(x) {
  n <- colSums(is.na(x[spatVars])); n[n > 0]
})
# no NAs. good to go. 

################################################################################
# 12) Niche Space PCA
################################################################################
# 19 Nov 2026. Run separately for each buffer version. Species and habitat
# variables are selected by name (spatVars from section 11), not by position.
# PC signs are arbitrary, so PC1/PC2 may be mirrored between buffer sizes; this
# does not affect areas or overlaps.

# species columns (species kept in section 9) -----------------------------------
speciesCols <- intersect(names(siteDetections_foliarTraits_BioCubeBuffers[["150"]]),
                         unique(siteCounts$common_name))
length(speciesCols) # 92 species

# site data, PCA, add PCs ------------------------------------------------------
# 20 Nov 2026: detection rates are no longer rounded to counts; section 15 uses
# them directly as kernel density weights. The PCA uses only the habitat
# variables, so it is unchanged.

siteDataBuffers <- lapply(siteDetections_foliarTraits_BioCubeBuffers, function(x) {
  x$geometry <- NULL
  as.data.frame(x)
})

traitsPCABuffers <- lapply(siteDataBuffers, function(x)
  prcomp(x[, spatVars], center = TRUE, scale = TRUE))

sitePCABuffers <- lapply(names(siteDataBuffers), function(b)
  cbind(siteDataBuffers[[b]], traitsPCABuffers[[b]]$x))
names(sitePCABuffers) <- names(siteDataBuffers)

# variance explained by PC1-PC4, for each buffer size (%)
sapply(traitsPCABuffers, function(p) round(100 * summary(p)$importance["Proportion of Variance", 1:min(10, ncol(p$x))], 1))
# PC1 23, PC2 18, PC3 8, PC4 7, PC5 4, PC6 4, PC7 3, PC8 3, PC9, 2, PC10 2

# save -------------------------------------------------------------------------
for (b in names(traitsPCABuffers)) {
  write_rds(traitsPCABuffers[[b]], paste0("traitsPCA_buffer", b, "m_20261120.rds"))
  write_rds(sitePCABuffers[[b]], paste0("sitePCA_buffer", b, "m_20261120.rds"))
}


################################################################################
# 13) PCA loadings (supp figures)
################################################################################
# 19 Nov 2026. Run for each buffer version; figure file names include the buffer size.

# variable categories (one row per variable: tidyNames.csv lists the trait and
# canopy height names twice, for the 30 m and 1 km files)
tidyCat <- read_csv("tidyNames.csv") %>% distinct(VariableName, Category)

catColors <- c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")

plot_loading <- function(loadings, pc) {
  ggplot(loadings, aes(x = reorder(variable, abs(.data[[pc]])), y = abs(.data[[pc]]), fill = Category)) +
    scale_fill_manual(values = catColors) +
    geom_col() + coord_flip() +
    xlab("Habitat Variable") +
    ylab(paste(pc, "loading (absolute value)")) +
    theme_minimal()
}

loadingsBuffers <- list()
for (b in names(traitsPCABuffers)) {
  p <- traitsPCABuffers[[b]]

  # what are the PCA loadings? -------------------------------------------------
  loadings <- as.data.frame(p$rotation[, 1:4])
  loadings$variable <- rownames(loadings)
  rownames(loadings) <- NULL
  loadings <- merge(loadings, tidyCat, by.x = "variable", by.y = "VariableName", all.x = TRUE)
  stopifnot(nrow(loadings) == length(spatVars), !anyNA(loadings$Category))
  loadingsBuffers[[b]] <- loadings

  # loadings plots PC1-PC4 - supplemental figures 1-4
  for (pc in paste0("PC", 1:4)) {
    print(plot_loading(loadings, pc))
    ggsave(paste0("Loading_", pc, "_buffer", b, "m_20261119.pdf"), height = 10, width = 7)
  }

  # visualize eigenvectors - supplemental figures 5 and 6
  print(fviz_pca_var(p, axes = c(1, 2)))
  ggsave(paste0("Eigenvectors_PC12_buffer", b, "m_20261119.pdf"), width = 10, height = 10)
  print(fviz_pca_var(p, axes = c(3, 4)))
  ggsave(paste0("Eigenvectors_PC34_buffer", b, "m_20261119.pdf"), width = 10, height = 10)
}

# contributions of variables to PC1 and PC2 (150 m)
fviz_contrib(traitsPCABuffers[["150"]], choice="var", axes = 1, sort.val ="asc") + 
  coord_flip() + labs(title = "Contribution of variables to PC1 (%), 150 m") +
  ylab("Habitat")

fviz_contrib(traitsPCABuffers[["150"]], choice="var", axes = 2, sort.val ="asc") + 
  coord_flip() + labs(title = "Contribution of variables to PC2 (%), 150 m") +
  ylab("Habitat")

# scree plot with broken-stick expectation - supplemental figure ---------------
# 20 Nov 2026. Supports keeping PC1-PC4 for the acoustic-use space: % variance
# levels off after PC4 (each later PC explains < 4.5%). Dashed line = broken-stick
# model (variance expected per PC if total variance were split at random).
screeDF <- map_dfr(names(traitsPCABuffers), function(b) {
  ev <- traitsPCABuffers[[b]]$sdev^2
  k  <- length(ev)
  tibble(buffer = paste(b, "m buffer"), PC = 1:k, pctVar = 100 * ev / sum(ev),
         cumVar = cumsum(pctVar), brokenStick = 100 * rev(cumsum(1 / (k:1))) / k)
}) %>%
  mutate(buffer = factor(buffer, levels = paste(names(traitsPCABuffers), "m buffer")))

screeDF %>% filter(PC <= 8) %>% mutate(across(where(is.numeric), ~ round(., 1))) %>% print(n = Inf)

ggplot(filter(screeDF, PC <= 20), aes(PC, pctVar)) +
  geom_vline(xintercept = 4.5, linetype = "dotted", colour = "grey50") +
  geom_line(aes(y = brokenStick), linetype = "dashed", colour = "grey40") +
  geom_line() + geom_point(size = 1.5) +
  facet_wrap(~ buffer) +
  scale_x_continuous(breaks = c(1, 4, 8, 12, 16, 20)) +
  labs(x = "Principal component", y = "Variance explained (%)") +
  theme_bw()
ggsave("ScreePlot_allBuffers_20261120.pdf", width = 9, height = 3.5)


################################################################################
# 14) render landcover labels
################################################################################
# 19 Nov 2026. Run for each buffer version. Sites are colored by GAP landcover
# type, and three on-figure labels (Oak, Pine, Mixed Conifer) are placed at the
# mean PC1/PC2 position of the sites in those landcover groups, so the labels
# follow the PCA of each buffer size (PC signs can differ between buffers - 300m buffer is mirrored to 50 and 150m).

# landcover group for the on-figure labels
lc_label <- function(lc) dplyr::case_when(
  grepl("Oak", lc)           ~ "Oak",   # checked first: "Blue Oak-Foothill Pine" counts as oak
  grepl("Pine", lc)          ~ "Pine",
  grepl("Mixed Conifer", lc) ~ "Mixed Conifer",
  TRUE                       ~ NA_character_)

landcoverLabelsBuffers <- list()
for (b in names(traitsPCABuffers)) {
  p  <- traitsPCABuffers[[b]]
  lc <- siteDataBuffers[[b]]$landcover_type # same rows, same order as the PCA

  # label positions: mean PC1/PC2 of the sites in each group
  lcLabels <- data.frame(p$x[, 1:2], label = lc_label(lc)) %>%
    filter(!is.na(label)) %>%
    group_by(label) %>%
    dplyr::summarise(x = mean(PC1), y = mean(PC2))
  landcoverLabelsBuffers[[b]] <- lcLabels

  # the plot
  landcover_pca_plot <- fviz_pca_ind(p,
                                     label = "none",
                                     habillage = lc) +
    guides(color = guide_legend(ncol = 1),
           fill = guide_legend(ncol = 1),
           shape = guide_legend(ncol = 1)) +
    geom_text(data = lcLabels, aes(x = x, y = y, label = label), inherit.aes = FALSE) +
    ggtitle(NULL)
  print(landcover_pca_plot)
  ggsave(paste0("landcover_pca_buffer", b, "m_20261119.pdf"), plot = landcover_pca_plot, height = 7, width = 12)

  # PC3 and PC4 panels - supplemental figure
  # 20 Nov 2026: shared landcover legend in the empty sixth panel, as in the
  # previous draft (Figure S11)
  make_pca_plot <- function(ax) {
    fviz_pca_ind(p,
                 axes = ax,
                 label = "none",
                 habillage = lc) +
      ggtitle(paste0("PC", ax[1], " vs PC", ax[2])) +
      guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2),
             shape = guide_legend(ncol = 2))
  }
  print(patchwork::wrap_plots(make_pca_plot(c(4,1)), make_pca_plot(c(3,1)),
                              make_pca_plot(c(4,2)), make_pca_plot(c(3,2)),
                              make_pca_plot(c(4,3)), patchwork::guide_area(), ncol = 2) +
          patchwork::plot_layout(guides = "collect") &
          theme(legend.text = element_text(size = 5), legend.key.size = unit(0.3, "cm")))
  ggsave(paste0("landcover_pca_PC34panels_buffer", b, "m_20261120.pdf"), height = 10, width = 11.5)
}


# buffer size sensitivity - supplemental figure ---------------------------------
# 20 Nov 2026. 50 m and 300 m buffers, PC1 x PC2 (with Oak/Pine/Mixed Conifer
# labels) and PC3 x PC4, with one shared landcover legend
buffer_panel <- function(b, ax) {
  # (all three buffers have the same 24 landcover types, so the legends merge)
  g <- fviz_pca_ind(traitsPCABuffers[[b]], axes = ax, label = "none",
                    habillage = siteDataBuffers[[b]]$landcover_type) +
    guides(color = guide_legend(ncol = 3), fill = guide_legend(ncol = 3),
           shape = guide_legend(ncol = 3)) +
    ggtitle(paste0(b, " m buffer: PC", ax[1], " vs PC", ax[2]))
  if (all(ax == c(1, 2))) {
    lab <- landcoverLabelsBuffers[[b]]  # Oak/Pine/Mixed Conifer positions from the loop above
    stopifnot(is.data.frame(lab), all(c("x", "y", "label") %in% names(lab)))
    g <- g + annotate("text", x = lab$x, y = lab$y, label = lab$label)
  }
  g
}
print(patchwork::wrap_plots(buffer_panel("50", c(1, 2)), buffer_panel("50", c(3, 4)),
                            buffer_panel("300", c(1, 2)), buffer_panel("300", c(3, 4)), ncol = 2) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(tag_levels = "A") &
        theme(legend.position = "bottom", legend.text = element_text(size = 6),
              legend.key.size = unit(0.3, "cm"), legend.title = element_blank()))
ggsave("landcover_pca_bufferSensitivity_20261120.pdf", height = 11, width = 10)

# 150 m label positions, in the format used by the overlap plots in section 16
landcover_centroids <- dplyr::rename(landcoverLabelsBuffers[["150"]], Groups = label)


################################################################################
# 15) Species acoustic-use space: weighted 4-D kernel density
################################################################################
# 20 Nov 2026. Replaces the rounded, jittered 2-D kde2d method (old sections
# 15 and 16). Run for each buffer size. For each species:
#   - space     = site scores on PC1-PC4 (raw PC scores, not rescaled)
#   - weights   = the species' detection rate at each site (section 3 rate)
#   - bandwidth = one diagonal plug-in bandwidth per buffer (ks::Hpi.diag),
#                 estimated from all sites and shared by all species, so
#                 volumes and overlaps are comparable among species
#   - grid      = one common grid per buffer; cell width on each PC is a fixed
#                 fraction of that PC's bandwidth (finer 4-D grids run out of
#                 memory on a laptop)
#   - acoustic-use space = smallest set of grid cells holding 90% of the density
#   - volume    = number of cells x cell volume (PC units^4)
#   - overlap   = Jaccard: cells shared by both species / cells in either
# 2-D regions (PC1 x PC2 and PC3 x PC4) are estimated the same way from their
# own 2-D densities and saved as sf polygons (holes kept as holes) for the
# figures and the Shiny app.
# Runtime: roughly 15 s per species per buffer (~25 min per buffer).

nPC        <- 4
prob       <- 0.90  # share of the density inside the acoustic-use space
cellFrac4D <- 0.75  # 4-D cell width, as a fraction of each PC's bandwidth
cellFrac2D <- 0.25  # 2-D cell width (2-D grids are cheap, so finer)

# functions --------------------------------------------------------------------
# grid points per axis = padded range / (cellFrac x bandwidth)
grid_points <- function(H, xmin, xmax, cellFrac)
  ceiling((xmax - xmin) / (cellFrac * sqrt(diag(H)))) + 1

# weighted KDE on a fixed grid, and the cells holding `prob` of the density
kde_region <- function(X, w, H, xmin, xmax, gridsize, prob = 0.90) {
  w  <- w * nrow(X) / sum(w)  # ks needs weights summing to n
  fh <- ks::kde(X, H = H, w = w, xmin = xmin, xmax = xmax, gridsize = gridsize,
                binned = TRUE, compute.cont = FALSE)
  est <- fh$estimate
  est[est < 0] <- 0           # binned estimates can have tiny negative values
  cellVol <- prod(sapply(fh$eval.points, function(g) diff(g[1:2])))
  z   <- sort(as.vector(est), decreasing = TRUE)
  thr <- z[which(cumsum(z) / sum(z) >= prob)[1]]
  inside <- est >= thr
  list(eval.points = fh$eval.points, estimate = est, threshold = thr,
       inside = inside, nCells = sum(inside), volume = sum(inside) * cellVol,
       massOnGrid = sum(est) * cellVol)  # ~1; lower means density off the grid
}

# 2-D region as an sf MULTIPOLYGON; isobands returns holes as interior rings
kde_polygon <- function(fit) {
  iso <- isoband::isobands(x = fit$eval.points[[1]], y = fit$eval.points[[2]],
                           z = t(fit$estimate), levels_low = fit$threshold,
                           levels_high = max(fit$estimate) + 1)
  st_cast(st_make_valid(st_sfc(isoband::iso_to_sfg(iso)[[1]])), "MULTIPOLYGON")
}

# 4-D regions are stored as packed bits (1 bit per cell) to save memory
pack_cells  <- function(v) packBits(c(v, rep(FALSE, (8 - length(v) %% 8) %% 8)))
bitsPerByte <- sapply(0:255, function(i) sum(as.integer(intToBits(i))))
count_cells <- function(r) sum(bitsPerByte[as.integer(r) + 1L])

# run for each buffer size -----------------------------------------------------
volumesBuffers <- list(); overlapsBuffers <- list()
polygonsPC12Buffers <- list(); polygonsPC34Buffers <- list()

for (b in names(traitsPCABuffers)) {
  pcs   <- traitsPCABuffers[[b]]$x[, 1:nPC]
  rates <- siteDataBuffers[[b]][, speciesCols]
  stopifnot(nrow(pcs) == nrow(rates))  # same sites, same order as the PCA

  # common bandwidths and grids for this buffer
  H4  <- ks::Hpi.diag(pcs)
  H12 <- ks::Hpi.diag(pcs[, 1:2])
  H34 <- ks::Hpi.diag(pcs[, 3:4])
  pad  <- 3 * sqrt(diag(H4))
  xmin <- apply(pcs, 2, min) - pad
  xmax <- apply(pcs, 2, max) + pad
  grid4D <- grid_points(H4,  xmin,      xmax,      cellFrac4D)
  grid12 <- grid_points(H12, xmin[1:2], xmax[1:2], cellFrac2D)
  grid34 <- grid_points(H34, xmin[3:4], xmax[3:4], cellFrac2D)
  cellVol4D <- prod((xmax - xmin) / (grid4D - 1))
  message(b, " m buffer: 4-D grid ", paste(grid4D, collapse = " x "),
          " (", format(prod(grid4D), big.mark = ","), " cells)")

  cells <- list(); vols <- list(); poly12 <- list(); poly34 <- list()
  for (sp in speciesCols) {
    w   <- rates[[sp]]
    f4  <- kde_region(pcs,         w, H4,  xmin,      xmax,      grid4D, prob)
    f12 <- kde_region(pcs[, 1:2], w, H12, xmin[1:2], xmax[1:2], grid12, prob)
    f34 <- kde_region(pcs[, 3:4], w, H34, xmin[3:4], xmax[3:4], grid34, prob)
    cells[[sp]]  <- pack_cells(as.vector(f4$inside))
    poly12[[sp]] <- kde_polygon(f12)
    poly34[[sp]] <- kde_polygon(f34)
    vols[[sp]] <- tibble(species = sp, nSites = sum(w > 0),
                         volume4D = f4$volume, nCells4D = f4$nCells,
                         area_PC12 = as.numeric(st_area(poly12[[sp]])),
                         area_PC34 = as.numeric(st_area(poly34[[sp]])),
                         mass4D = f4$massOnGrid)
    message(b, " m: ", sp)
  }
  volumes <- bind_rows(vols)
  stopifnot(all(volumes$mass4D > 0.99))  # density not cut off by the grid edge

  # pairwise 4-D overlaps (Jaccard) -------------------------------------------
  nCells <- setNames(volumes$nCells4D, volumes$species)
  overlaps <- as.data.frame(t(combn(sort(speciesCols), 2)))
  names(overlaps) <- c("species1", "species2")
  overlaps$overlap_cells <- mapply(function(s1, s2) count_cells(cells[[s1]] & cells[[s2]]),
                                   overlaps$species1, overlaps$species2)
  overlaps <- overlaps %>%
    mutate(volume1 = nCells[species1] * cellVol4D,
           volume2 = nCells[species2] * cellVol4D,
           overlap_volume = overlap_cells * cellVol4D,
           percentOverlap = overlap_cells / (nCells[species1] + nCells[species2] - overlap_cells)) %>%
    dplyr::select(-overlap_cells)

  # 2-D polygons as sf -----------------------------------------------------------
  polys12 <- st_sf(species = names(poly12), geometry = do.call(c, poly12))
  polys34 <- st_sf(species = names(poly34), geometry = do.call(c, poly34))

  volumesBuffers[[b]] <- volumes; overlapsBuffers[[b]] <- overlaps
  polygonsPC12Buffers[[b]] <- polys12; polygonsPC34Buffers[[b]] <- polys34

  # save -----------------------------------------------------------------------
  write_csv(volumes,  paste0("acousticSpace_volumes_buffer", b, "m_20261120.csv"))
  write_csv(overlaps, paste0("acousticSpace_4Doverlaps_buffer", b, "m_20261120.csv"))
  write_rds(polys12,  paste0("speciesPolygons_PC12_buffer", b, "m_20261120.rds"))
  write_rds(polys34,  paste0("speciesPolygons_PC34_buffer", b, "m_20261120.rds"))
  write_rds(list(cells = cells, grid4D = grid4D, xmin = xmin, xmax = xmax, H4 = H4,
                 cellVol4D = cellVol4D),
            paste0("acousticSpace_4Dcells_buffer", b, "m_20261120.rds"))
}

# checks -----------------------------------------------------------------------
lapply(volumesBuffers, function(v) summary(v$volume4D))
lapply(overlapsBuffers, function(o) summary(o$percentOverlap))
nrow(overlapsBuffers[["150"]]) # 4,186 pairs for 92 species

# check figure: 2-D polygons for four species (150 m), sites sized by rate
checkSp <- c("Hermit Thrush", "Swainson's Thrush", "Yellow-rumped Warbler", "Oak Titmouse")
checkSites <- map_dfr(checkSp, ~ tibble(species = .x, PC1 = traitsPCABuffers[["150"]]$x[, 1],
                                        PC2 = traitsPCABuffers[["150"]]$x[, 2],
                                        rate = siteDataBuffers[["150"]][[.x]]))
ggplot() +
  geom_sf(data = filter(polygonsPC12Buffers[["150"]], species %in% checkSp),
          fill = "#27548A", alpha = 0.3, colour = "#27548A") +
  geom_point(data = filter(checkSites, rate == 0), aes(PC1, PC2), colour = "grey70", size = 0.3) +
  geom_point(data = filter(checkSites, rate > 0), aes(PC1, PC2, size = rate), shape = 21) +
  scale_size_area(max_size = 3) +
  facet_wrap(~ species) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal()


################################################################################
# 16) plot 2 species acoustic-use space overlap (test plot)
################################################################################
# 23 Nov 2026. Checks the 2-D polygons from section 15 for any pair of species.
# Panels show PC1 x PC2 and PC3 x PC4. Axes use equal units (1 PC unit is the
# same length on both axes), which is how areas and volumes were calculated;
# PC1 looks longer than PC2 because its scores are more spread out.
# Colored text under each panel gives the areas and Jaccard % overlap; the
# block below gives the 4-D volumes and overlap.

# choose buffer and species ----------------------------------------------------
b        <- "150"
Species1 <- "White-breasted Nuthatch"
Species2 <- "Red-breasted Nuthatch"
sort(speciesCols)

# colors as in the previous draft
col1 <- "#25b9be"; col2 <- "#e41a1c"; colOverlap <- "#8849a5"

# polygons, overlaps and numbers -----------------------------------------------
pct <- round(100 * summary(traitsPCABuffers[[b]])$importance["Proportion of Variance", 1:4], 1)
sites <- as.data.frame(traitsPCABuffers[[b]]$x[, 1:4])

overlap_2D <- function(polys) {
  p1 <- st_make_valid(st_geometry(filter(polys, species == Species1)))
  p2 <- st_make_valid(st_geometry(filter(polys, species == Species2)))
  shared <- st_intersection(p1, p2)
  a1 <- sum(st_area(p1)); a2 <- sum(st_area(p2))
  aShared <- if (length(shared) == 0) 0 else sum(st_area(shared))
  list(p1 = p1, p2 = p2, shared = shared, area1 = a1, area2 = a2,
       jaccard = aShared / (a1 + a2 - aShared))
}
ov12 <- overlap_2D(polygonsPC12Buffers[[b]])
ov34 <- overlap_2D(polygonsPC34Buffers[[b]])

vol <- volumesBuffers[[b]]
ov4D <- overlapsBuffers[[b]] %>%
  filter((species1 == Species1 & species2 == Species2) |
           (species1 == Species2 & species2 == Species1))
stopifnot(nrow(ov4D) == 1)

# plotting function ------------------------------------------------------------
overlap_panel <- function(ov, ax) {
  g <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(data = sites, aes(x = .data[[paste0("PC", ax[1])]], y = .data[[paste0("PC", ax[2])]]),
               colour = "grey50", alpha = 0.2) +
    geom_sf(data = st_sf(geometry = ov$p1), fill = col1, alpha = 0.6, colour = NA) +
    geom_sf(data = st_sf(geometry = ov$p2), fill = col2, alpha = 0.6, colour = NA)
  if (length(ov$shared) > 0)
    g <- g + geom_sf(data = st_sf(geometry = ov$shared), fill = colOverlap, alpha = 0.6, colour = NA)
  if (all(ax == c(1, 2))) {
    lab <- landcoverLabelsBuffers[[b]]
    g <- g + annotate("text", x = lab$x, y = lab$y, label = lab$label)
  }
  g + coord_sf(expand = TRUE) +
    labs(x = paste0("PC", ax[1], " (", pct[ax[1]], "%)"),
         y = paste0("PC", ax[2], " (", pct[ax[2]], "%)")) +
    theme_minimal()
}

# colored text under each panel, as in the range maps (section 19) ---------------
text_block <- function(header, lines, cols) {
  n <- length(lines)
  g <- ggplot() + annotate("text", x = 0, y = n + 1, hjust = 0, size = 4, label = header)
  for (i in seq_len(n))
    g <- g + annotate("text", x = 0, y = n + 1 - i, hjust = 0, size = 4, colour = cols[i], label = lines[i])
  g + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0.5, n + 1.5)) + theme_void()
}
fmt <- function(x, d = 1) format(round(as.numeric(x), d), nsmall = d, big.mark = ",")
cols3 <- c(col1, col2, colOverlap)

v1 <- vol$volume4D[vol$species == Species1]
v2 <- vol$volume4D[vol$species == Species2]
txt12 <- text_block("PC1 x PC2 (area in PC units2)",
                    c(paste0(Species1, " area: ", fmt(ov12$area1)),
                      paste0(Species2, " area: ", fmt(ov12$area2)),
                      paste0("Overlap: ", fmt(100 * ov12$jaccard), "%")), cols3)
txt34 <- text_block("PC3 x PC4 (area in PC units2)",
                    c(paste0(Species1, " area: ", fmt(ov34$area1)),
                      paste0(Species2, " area: ", fmt(ov34$area2)),
                      paste0("Overlap: ", fmt(100 * ov34$jaccard), "%")), cols3)
txt4D <- text_block(paste0("4-D, PC1 to PC4 (volume in PC units4), ", b, " m buffer"),
                    c(paste0(Species1, " volume: ", fmt(v1, 0)),
                      paste0(Species2, " volume: ", fmt(v2, 0)),
                      paste0("Overlap: ", fmt(100 * ov4D$percentOverlap), "%")), cols3)

print(patchwork::wrap_plots(list(A = overlap_panel(ov12, c(1, 2)), B = overlap_panel(ov34, c(3, 4)),
                                 C = txt12, D = txt34, E = txt4D),
                            design = "AB\nCD\nE#", heights = c(5, 1.3, 1.3)))
ggsave(paste0("OverlapPoly_", gsub("[^A-Za-z]", "", Species1), "_", gsub("[^A-Za-z]", "", Species2),
              "_buffer", b, "m_20261124.pdf"), height = 7.5, width = 11)


################################################################################
# 17) Species niche areas and pairwise overlaps: PC1 x PC2, PC3 x PC4 and 4-D
################################################################################
# 23 Nov 2026. Section 15 already calculated, for each buffer, every species'
# PC1 x PC2 area, PC3 x PC4 area and 4-D volume, and the pairwise 4-D overlaps.
# This section adds the pairwise 2-D overlaps (from the polygons, holes
# excluded) and saves one species table and one pairwise table per buffer.
# % overlap = Jaccard: shared area (or volume) / area covered by either species,
# x 100. Units: areas in PC units^2, volumes in PC units^4.
# Runtime: a few minutes per buffer.

# pairwise 2-D overlap areas from a set of species polygons
pair_overlap_2D <- function(polys, pairs) {
  geoms <- st_make_valid(st_geometry(polys))
  i1 <- match(pairs$species1, polys$species)
  i2 <- match(pairs$species2, polys$species)
  stopifnot(!anyNA(i1), !anyNA(i2))
  mapply(function(a, b) {
    shared <- st_intersection(geoms[a], geoms[b])
    if (length(shared) == 0) 0 else sum(as.numeric(st_area(shared)))
  }, i1, i2)
}

speciesAreasBuffers <- list(); pairwiseOverlapsBuffers <- list()
for (b in names(volumesBuffers)) {
  # species table --------------------------------------------------------------
  areas <- volumesBuffers[[b]] %>%
    dplyr::select(species, nSites, area_PC12, area_PC34, volume4D)
  speciesAreasBuffers[[b]] <- areas

  # pairwise table: start from the 4-D overlaps (all species pairs) -----------
  pairs <- overlapsBuffers[[b]] %>%
    transmute(species1, species2,
              volume1_4D = volume1, volume2_4D = volume2,
              overlapVolume_4D = overlap_volume,
              percentOverlap_4D = 100 * percentOverlap)

  # add 2-D areas and overlaps
  for (plane in c("PC12", "PC34")) {
    polys <- if (plane == "PC12") polygonsPC12Buffers[[b]] else polygonsPC34Buffers[[b]]
    a <- setNames(areas[[paste0("area_", plane)]], areas$species)
    shared <- pair_overlap_2D(polys, pairs)
    pairs[[paste0("area1_", plane)]] <- a[pairs$species1]
    pairs[[paste0("area2_", plane)]] <- a[pairs$species2]
    pairs[[paste0("overlapArea_", plane)]] <- shared
    pairs[[paste0("percentOverlap_", plane)]] <- 100 * shared /
      (a[pairs$species1] + a[pairs$species2] - shared)
    message(b, " m: ", plane, " done")
  }
  pairs <- pairs %>%
    dplyr::select(species1, species2,
                  area1_PC12, area2_PC12, overlapArea_PC12, percentOverlap_PC12,
                  area1_PC34, area2_PC34, overlapArea_PC34, percentOverlap_PC34,
                  volume1_4D, volume2_4D, overlapVolume_4D, percentOverlap_4D)
  pairwiseOverlapsBuffers[[b]] <- pairs

  # save -----------------------------------------------------------------------
  write_csv(areas, paste0("speciesAreas_buffer", b, "m_20261123.csv"))
  write_csv(pairs, paste0("pairwiseOverlaps_buffer", b, "m_20261123.csv"))
}

# checks -----------------------------------------------------------------------
sapply(pairwiseOverlapsBuffers, nrow) # 4,186 pairs per buffer for 92 species
lapply(pairwiseOverlapsBuffers, function(p)
  summary(p[, c("percentOverlap_PC12", "percentOverlap_PC34", "percentOverlap_4D")]))
# how closely do 2-D (PC1 x PC2) and 4-D overlaps agree?
sapply(pairwiseOverlapsBuffers, function(p) cor(p$percentOverlap_PC12, p$percentOverlap_4D)) # more than 90% similar.


################################################################################
# 18) Range areas and pairwise range overlaps
################################################################################
# 23 Nov 2026. eBird Status and Trends range maps (RangeMaps/<species>/*_range_<year>.gpkg).
# - Breeding range only: "breeding" polygons for migratory species, "resident"
#   polygons for sedentary species.
# - Species names are taken from the folder names, which match speciesCols
#   (common names with spaces, as everywhere else in this script). The
#   sapsucker complex uses the union of the Red-breasted and Red-naped ranges.
# - Areas are calculated in EPSG 6933 (WGS 84 / EASE-Grid 2.0, equal-area,
#   global). Some breeding ranges extend outside North America and across the
#   180th meridian (e.g. Common Raven, Pacific Wren). CONUS Albers (EPSG 5070)
#   is equal-area in principle, but these polygons are not projected correctly
#   in it: it overstates Common Raven's range by ~16%; EPSG 6933 is within 0.3%
#   of the geodesic area.
# - Extent sensitivity (reviewer request): ranges are clipped to
#     full    : no clipping (whole breeding range)
#     NorthAm : North America north of 15 N (lon -170 to -50, lat 15 to 84)
#     WestNA  : western North America (as NorthAm, west of 100 W)
#     study   : square around the sites (side = the larger of the sites'
#               east-west and north-south extent)
# - % overlap = Jaccard x 100, as in section 17. Areas in km^2.

# extents ----------------------------------------------------------------------
sf_use_s2(FALSE) # planar geometry in EPSG 6933
equalArea <- 6933

lonlat_box <- function(xmin, xmax, ymin, ymax) {
  st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = 4326)) %>%
    st_transform(equalArea) %>%
    st_segmentize(10000) # straight in EPSG 6933; vertices every 10 km so the edges curve correctly in the map projection
}

sites <- st_as_sf(siteDataBuffers[["150"]], coords = c("long", "lat"), crs = 4326)
sb <- st_bbox(st_transform(sites, 3310))   # California Albers, meters
side <- max(sb["xmax"] - sb["xmin"], sb["ymax"] - sb["ymin"])
cx <- mean(sb[c("xmin", "xmax")]); cy <- mean(sb[c("ymin", "ymax")])
studySquare <- st_as_sfc(st_bbox(c(xmin = cx - side / 2, xmax = cx + side / 2,
                                   ymin = cy - side / 2, ymax = cy + side / 2), crs = 3310)) %>%
  st_segmentize(10000) %>% st_transform(equalArea)

extents <- list(full    = NULL,
                NorthAm = lonlat_box(-170, -50, 15, 84),
                WestNA  = lonlat_box(-170, -100, 15, 84),
                study   = studySquare)
round(side / 1000) # side of the study square (km)

# breeding ranges, one (multi)polygon per species --------------------------------
rangeFiles <- list.files(file.path(base_dir, "RangeMaps"), pattern = "_range_20(22|23)\\.gpkg$",
                         full.names = TRUE, recursive = TRUE)
ranges <- map_dfr(rangeFiles, function(f)
  st_read(f, quiet = TRUE) %>%
    filter(season %in% c("breeding", "resident")) %>%
    mutate(species = basename(dirname(f))) %>%
    dplyr::select(species, scientific_name, season))
ranges$species[ranges$species == "Eurasian Collared Dove"] <- "Eurasian Collared-Dove"

ranges <- ranges %>%
  st_transform(equalArea) %>%
  st_make_valid() %>%
  group_by(species) %>%
  dplyr::summarise(scientific_name = first(scientific_name), .groups = "drop")

# species complexes: union of the member ranges (complexMembers, section 2)
for (cx_name in unique(complexMembers$common_name)) {
  members <- ranges %>% filter(scientific_name %in% complexMembers$member[complexMembers$common_name == cx_name])
  stopifnot(nrow(members) == sum(complexMembers$common_name == cx_name))
  # summarise() unions the member geometries and keeps the same geometry column
  ranges <- bind_rows(ranges, members %>%
                        dplyr::summarise(species = cx_name, scientific_name = NA_character_))
}

setdiff(speciesCols, ranges$species) # should be empty
stopifnot(!any(st_is_empty(ranges)), length(attr(ranges, "sf_column")) == 1,
          sum(sapply(ranges, inherits, "sfc")) == 1) # one geometry column, no empty ranges
ranges <- ranges %>% filter(species %in% speciesCols) %>% arrange(species)
write_rds(ranges, "speciesRanges_breeding_EPSG6933_20261123.rds")

# areas and pairwise overlaps for each extent ------------------------------------
pairs <- as.data.frame(t(combn(sort(speciesCols), 2)))
names(pairs) <- c("species1", "species2")

rangeAreas <- tibble(species = ranges$species)
rangeOverlapsList <- list()
for (e in names(extents)) {
  g <- st_geometry(ranges)
  if (!is.null(extents[[e]])) g <- st_intersection(g, extents[[e]], dimensions = "polygon")
  # st_intersection drops empty results; put every species back in order
  if (length(g) != nrow(ranges)) {
    keep <- st_intersects(st_geometry(ranges), extents[[e]], sparse = FALSE)[, 1]
    gFull <- rep(st_sfc(st_multipolygon(), crs = equalArea), nrow(ranges))
    gFull[keep] <- g
    g <- gFull
  }
  areas <- setNames(as.numeric(st_area(g)) / 1e6, ranges$species)
  rangeAreas[[paste0("area_km2_", e)]] <- areas

  i1 <- match(pairs$species1, ranges$species); i2 <- match(pairs$species2, ranges$species)
  touching <- st_intersects(g)
  shared <- mapply(function(a, b) {
    if (!(b %in% touching[[a]])) return(0)
    as.numeric(sum(st_area(st_intersection(g[a], g[b])))) / 1e6
  }, i1, i2)

  rangeOverlapsList[[e]] <- pairs %>%
    mutate(extent = e, area1_km2 = areas[species1], area2_km2 = areas[species2],
           overlapArea_km2 = shared,
           percentOverlap = ifelse(area1_km2 + area2_km2 - shared > 0,
                                   100 * shared / (area1_km2 + area2_km2 - shared), NA))
  message(e, " done")
}
rangeOverlaps <- bind_rows(rangeOverlapsList)

write_csv(rangeAreas, "rangeAreas_20261123.csv")
write_csv(rangeOverlaps, "rangeOverlaps_20261123.csv") # long format: one row per pair and extent

# checks -----------------------------------------------------------------------
rangeAreas %>% filter(area_km2_study == 0)  # species with no breeding range in the study square
rangeOverlaps %>% group_by(extent) %>%
  dplyr::summarise(pairs = n(), meanOverlap = mean(percentOverlap, na.rm = TRUE),
                   noOverlap = sum(overlapArea_km2 == 0))
# how much does the extent change the % overlaps? (correlation with full ranges)
# 1 species pair never overlaps (ChBaChickadee + WHScrubJay)
# 10 species pairs don't overlap inside the study region (Tahoe, Yosemite) 4 of those are w/ Red-naped Sapsucker which should be merged.
rangeOverlaps %>% dplyr::select(species1, species2, extent, percentOverlap) %>%
  pivot_wider(names_from = extent, values_from = percentOverlap) %>%
  dplyr::select(-species1, -species2) %>% cor(use = "pairwise.complete.obs") %>% round(2)
# Full and NorthAm correlate 96%. Just the study region correlates only 46-58% with other range bboxes.
# using studyregion bbox might give me very different results.

# supplemental figure: the extents ------------------------------------------------
naAlbers <- "ESRI:102008" # North America Albers Equal Area, for display
countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  st_crop(c(xmin = -180, xmax = -20, ymin = 0, ymax = 90)) # keeps the projection well-behaved
extentLines <- st_sf(extent = c("North America", "Western North America", "Study area"),
                     geometry = c(extents$NorthAm, extents$WestNA, extents$study)) %>%
  st_transform(naAlbers)
extentLines$extent <- factor(extentLines$extent, levels = extentLines$extent)
view <- st_bbox(st_transform(extents$NorthAm, naAlbers))

ggplot() +
  geom_sf(data = st_transform(countries, naAlbers), fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = extentLines, aes(colour = extent, linetype = extent), fill = NA, linewidth = 0.8) +
  geom_sf(data = st_transform(sites, naAlbers), size = 0.3, colour = "black") +
  scale_colour_manual(values = c("#27548A", "#e41a1c", "#F5C45E")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid")) +
  coord_sf(xlim = view[c("xmin", "xmax")], ylim = view[c("ymin", "ymax")]) +
  labs(colour = "Range extent", linetype = "Range extent") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("FigureS_rangeExtents_20261123.pdf", height = 7, width = 7)


# are range overlaps zero-inflated? (reviewer comment) ---------------------------
# 24 Nov 2026. Distribution of pairwise % range overlap for each extent, for the
# 89 analysis species (excluded species as in section 20). All zeros are kept:
#   - ranges that do not touch (shared area 0)
#   - pairs where neither range reaches the extent (study extent only; overlap
#     is 0/0 = NA in rangeOverlaps), counted here as 0% overlap
# Zeros are shown as a separate colour in the first bar of each histogram.
histSpecies <- setdiff(speciesCols, c("Woodhouse's Scrub-Jay", "Red-breasted Sapsucker",
                                      "Red-naped Sapsucker")) # = analysisSpecies, section 20
histData <- rangeOverlaps %>%
  filter(species1 %in% histSpecies, species2 %in% histSpecies) %>%
  mutate(extent = factor(extent, levels = names(extents)),
         noRangeInExtent = is.na(percentOverlap),
         overlap = ifelse(noRangeInExtent, 0, percentOverlap),
         zero = overlapArea_km2 == 0)

zeroSummary <- histData %>%
  group_by(extent) %>%
  dplyr::summarise(pairs = n(),
                   zeroOverlap = sum(zero),                # all pairs with no shared range
                   pctZero = round(100 * zeroOverlap / pairs, 1),
                   of_which_noRangeInExtent = sum(noRangeInExtent),
                   exact100 = sum(overlap >= 99.999),
                   meanOverlap = round(mean(overlap), 1),
                   medianOverlap = round(median(overlap), 1),
                   .groups = "drop")
zeroSummary

ggplot(histData, aes(x = overlap, fill = zero)) +
  geom_histogram(binwidth = 2, boundary = 0, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = c(`FALSE` = "grey40", `TRUE` = "#e41a1c"),
                    labels = c(`FALSE` = "Overlap > 0", `TRUE` = "No shared range (0%)"), name = NULL) +
  geom_text(data = zeroSummary, aes(x = 98, y = Inf,
                label = paste0("Pairs with 0% overlap: ", zeroOverlap, 
                               "\nPairs with 100% overlap: ", exact100,
                               "\nMean overlap: ", meanOverlap, "%")),
            hjust = 1, vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  facet_wrap(~ extent, scales = "free_y") +
  labs(x = "Breeding range overlap (%)", y = "Species pairs") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("FigureS_rangeOverlapHistograms_20261124.pdf", height = 6, width = 8)


################################################################################
# 19) Range overlap maps (two species)
################################################################################
# 23 Nov 2026. Breeding ranges of any two species from section 18, as in the
# Shiny app: left, North America (NorthAm extent); right, the study square.
# Colors match the acoustic-use space plots (section 16), semi-transparent;
# purple = overlap. North America is drawn on flat grey land clipped to the
# extent; the study square on the ESRI physical basemap used for the site map
# (section 10).
# Under each map: both range areas and the % overlap (Jaccard) for that extent.
# Uses `ranges`, `extents`, `sites`, `countries`, `rangeAreas` and
# `rangeOverlaps` from section 18.

# choose species ---------------------------------------------------------------
sort(ranges$species)
Species1 <- "Cassin's Vireo"
Species2 <- "Plumbeous Vireo"

col1 <- "#25b9be"; col2 <- "#e41a1c"; colOverlap <- "#8849a5"
naAlbers <- "ESRI:102008" # North America Albers Equal Area
caAlbers <- "EPSG:3310"   # California Albers, for the study square
rangeAlpha <- 0.5

# study-area basemap (downloaded once; reused when you change species) ----------
# ESRI world physical map tiles (Web Mercator), reprojected to the map projection
get_basemap <- function(box, crs) {
  view <- st_as_sfc(st_bbox(st_transform(box, crs))) %>% st_segmentize(20000)
  bm <- basemaps::basemap_terra(ext = st_transform(view, 3857), map_service = "esri",
                                map_type = "world_physical_map", verbose = FALSE)
  bm <- terra::project(bm[[1:3]], crs, method = "bilinear")
  df <- as.data.frame(bm, xy = TRUE, na.rm = TRUE)
  names(df)[3:5] <- c("r", "g", "b")
  df$fill <- rgb(pmin(pmax(df$r, 0), 255), pmin(pmax(df$g, 0), 255), pmin(pmax(df$b, 0), 255),
                 maxColorValue = 255)
  df[c("x", "y", "fill")]
}
if (!exists("basemapStudy")) basemapStudy <- get_basemap(extents$study, caAlbers)

# functions --------------------------------------------------------------------
# clip to the extent in EPSG 6933 (where the areas were calculated), then project
clip_to <- function(g, box, crs) st_transform(st_intersection(st_geometry(g), box), crs)

range_map <- function(box, crs, basemap, showSites) {
  r1 <- clip_to(filter(ranges, species == Species1), box, crs)
  r2 <- clip_to(filter(ranges, species == Species2), box, crs)
  shared <- st_intersection(r1, r2)
  view <- st_bbox(st_transform(box, crs))
  outside <- st_difference(st_as_sfc(view), st_transform(box, crs)) # grey out beyond the extent
  g <- ggplot()
  if (is.null(basemap)) {
    # flat grey land, clipped to the extent
    land <- clip_to(st_transform(countries, equalArea), box, crs)
    g <- g + geom_sf(data = land, fill = "grey40", colour = "grey40")
  } else {
    g <- g + geom_raster(data = basemap, aes(x = x, y = y, fill = fill)) +
      scale_fill_identity() +
      geom_sf(data = outside, fill = "white", colour = NA, alpha = 0.7)
  }
  g <- g +
    geom_sf(data = r1, fill = col1, colour = NA, alpha = rangeAlpha) +
    geom_sf(data = r2, fill = col2, colour = NA, alpha = rangeAlpha)
  if (length(shared) > 0) g <- g + geom_sf(data = shared, fill = colOverlap, colour = NA, alpha = rangeAlpha)
  if (showSites) g <- g + geom_sf(data = st_transform(sites, crs), size = 0.3, colour = "black")
  g + coord_sf(crs = crs, xlim = view[c("xmin", "xmax")], ylim = view[c("ymin", "ymax")], expand = FALSE) +
    theme_void()
}

format_area <- function(km2) ifelse(km2 >= 1e6, paste(format(round(km2 / 1e6, 3), nsmall = 3), "million km2"),
                                    paste(format(round(km2), big.mark = ","), "km2"))

# extent name (black) and three colored lines of text for one extent
range_text <- function(e, extentLabel) {
  a1 <- rangeAreas[[paste0("area_km2_", e)]][rangeAreas$species == Species1]
  a2 <- rangeAreas[[paste0("area_km2_", e)]][rangeAreas$species == Species2]
  p <- rangeOverlaps %>%
    filter(extent == e, (species1 == Species1 & species2 == Species2) |
             (species1 == Species2 & species2 == Species1)) %>% pull(percentOverlap)
  stopifnot(length(a1) == 1, length(a2) == 1, length(p) == 1)
  ggplot() +
    annotate("text", x = 0, y = 4, hjust = 0, size = 4.5, colour = "black",
             label = paste0("Extent: ", extentLabel)) +
    annotate("text", x = 0, y = 3, hjust = 0, size = 4.5, colour = col1,
             label = paste0(Species1, " range area: ", format_area(a1))) +
    annotate("text", x = 0, y = 2, hjust = 0, size = 4.5, colour = col2,
             label = paste0(Species2, " range area: ", format_area(a2))) +
    annotate("text", x = 0, y = 1, hjust = 0, size = 4.5, colour = colOverlap,
             label = paste0("Range overlap: ", round(p, 1), "%")) +
    scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0.5, 4.5)) +
    theme_void()
}

# figure -----------------------------------------------------------------------
mapNA    <- range_map(extents$NorthAm, naAlbers, NULL,         showSites = FALSE) # flat grey land
mapStudy <- range_map(extents$study,   caAlbers, basemapStudy, showSites = TRUE)

print(patchwork::wrap_plots(mapNA, mapStudy, range_text("NorthAm", "North America"),
                            range_text("study", "Study area (Tahoe and Yosemite)"),
                            ncol = 2, widths = c(1.6, 1), heights = c(5, 1.3)))
ggsave(paste0("RangeOverlap_", gsub("[^A-Za-z]", "", Species1), "_", gsub("[^A-Za-z]", "", Species2),
              "_20261123.pdf"), height = 7, width = 12)


################################################################################
# 20) Combine range overlaps and niche overlaps
################################################################################
# 24 Nov 2026. One table for every combination of buffer size (section 17,
# pairwiseOverlapsBuffers) and range extent (section 18, rangeOverlapsList),
# with one row per species pair. Before joining, each pair is put in
# alphabetical order (species1 < species2) in both tables, and the
# species-specific columns (area1/area2, volume1/volume2) are swapped with the
# names, so a pair listed as (A, B) in one table and (B, A) in the other still
# matches correctly. Congeneric pairs and the remaining sister-species pairs are
# marked for sections 21 onward. All % overlaps are Jaccard x 100.

# species excluded from the analyses -------------------------------------------
# 23 Nov 2026. Woodhouse's Scrub-Jay: detected at 321 sites, none inside its
# eBird breeding range; 21% of its detections share a BirdNET segment with
# California Scrub-Jay (crossIdentification_species_20261124.csv), and away from
# its range it is detected at ~3-6% of the California Scrub-Jay rate, so most
# detections are likely misidentified California Scrub-Jays. California
# Scrub-Jay is kept (0.8% of its detections shared). Removing a species does not
# change any other species' results (PCA, bandwidths and grids use only the
# habitat variables), so earlier sections were not re-run; its rows are dropped
# here from the overlap tables.
# Red-breasted and Red-naped Sapsucker are analysed as the complex
# "Red-breasted/Red-naped Sapsucker" (section 2), so the two members are excluded.
excludedSpecies <- c("Woodhouse's Scrub-Jay", "Red-breasted Sapsucker", "Red-naped Sapsucker")
analysisSpecies <- setdiff(speciesCols, excludedSpecies)
length(analysisSpecies) # 89 species
nPairs <- choose(length(analysisSpecies), 2) # 3,916 pairs

# genus of each species (from the scientific names; the sapsucker complex is
# "Sphyrapicus ruber/nuchalis", so its genus is Sphyrapicus)
sciNames <- distinct(commonNames, common_name, scientific_name)
genus <- setNames(word(sciNames$scientific_name, 1), sciNames$common_name)
stopifnot(all(analysisSpecies %in% names(genus)))

# sister-species pairs among the analysis species (alphabetical order)
sisterPairs <- c("Cassin's Finch|Purple Finch", "Cassin's Vireo|Plumbeous Vireo")

# put each pair in alphabetical order, swapping the paired columns with it ------
order_pair <- function(df, pairedCols) {
  df <- as.data.frame(df)
  swap <- df$species1 > df$species2
  df[swap, c("species1", "species2")] <- df[swap, c("species2", "species1")]
  for (pc in pairedCols) df[swap, pc] <- df[swap, rev(pc)]
  df
}
keep_analysis <- function(df) filter(df, species1 %in% analysisSpecies, species2 %in% analysisSpecies)

nichePaired <- list(c("area1_PC12", "area2_PC12"), c("area1_PC34", "area2_PC34"),
                    c("volume1_4D", "volume2_4D"))
rangePaired <- list(c("range_area1_km2", "range_area2_km2"))

# combine ----------------------------------------------------------------------
nicheRangeOverlaps <- list()
for (b in names(pairwiseOverlapsBuffers)) {
  niche <- pairwiseOverlapsBuffers[[b]] %>% keep_analysis() %>% order_pair(nichePaired)
  for (e in names(rangeOverlapsList)) {
    range <- rangeOverlapsList[[e]] %>%
      keep_analysis() %>%
      transmute(species1, species2,
                range_area1_km2 = area1_km2, range_area2_km2 = area2_km2,
                range_overlapArea_km2 = overlapArea_km2,
                range_percentOverlap = percentOverlap) %>%
      order_pair(rangePaired)
    stopifnot(nrow(niche) == nPairs, nrow(range) == nPairs,
              !anyDuplicated(niche[c("species1", "species2")]),
              !anyDuplicated(range[c("species1", "species2")]))

    combined <- inner_join(niche, range, by = c("species1", "species2")) %>%
      mutate(buffer_m = as.numeric(b), extent = e,
             genus1 = unname(genus[species1]), genus2 = unname(genus[species2]),
             congeneric = if_else(genus1 == genus2, "Same Genus", "Different Genus"),
             sisterSpecies = paste(species1, species2, sep = "|") %in% sisterPairs) %>%
      relocate(buffer_m, extent, genus1, genus2, congeneric, sisterSpecies, .after = species2)
    stopifnot(nrow(combined) == nPairs) # every pair matched once

    nicheRangeOverlaps[[paste0("buffer", b, "m_", e)]] <- combined
  }
}

# checks -----------------------------------------------------------------------
names(nicheRangeOverlaps)  # 3 buffers x 4 extents = 12 tables
sapply(nicheRangeOverlaps, nrow)
table(nicheRangeOverlaps[[1]]$congeneric)    # 37 Same Genus, 3,879 Different Genus
sum(nicheRangeOverlaps[[1]]$sisterSpecies)   # 2
# pairs whose ranges do not reach the extent at all have NA range overlap - none of these
sapply(nicheRangeOverlaps, function(x) sum(is.na(x$range_percentOverlap)))

# save -------------------------------------------------------------------------
write_rds(nicheRangeOverlaps, "nicheRangeOverlaps_allBuffersExtents_20261124.rds")
write_csv(bind_rows(nicheRangeOverlaps), "nicheRangeOverlaps_allBuffersExtents_20261124.csv")

# main analysis: 150 m buffer, full breeding ranges (used in sections 21 onward)
NicheRangeOverlap <- nicheRangeOverlaps[["buffer150m_full"]]


################################################################################
# 21) Mixed effects models
################################################################################
# 24 Nov 2026. One Bayesian mixed model for each of the 12 datasets from
# section 20 (3 buffer sizes x 4 range extents), as before:
#   range overlap ~ niche overlap * congeneric,  random = ~ mm(species1 + species2)
# - Niche overlap is the 4-D acoustic-use space overlap; range overlap is the
#   breeding range overlap for that extent. Both are Jaccard overlaps, used as
#   proportions (0-1), the scale the priors were set for and that earlier
#   results were reported on.
# - No species pairs are excluded any more: the cross-identified taxa were
#   handled in sections 2 and 20, so all 37 congeneric pairs are included.
# - Pairs whose ranges both fall outside the extent (range overlap NA; study
#   extent only) are dropped.
# - The multiple-membership random effect accounts for each species appearing
#   in many pairs. The sapsucker complex is one species level.
# - 4 chains per model, run in parallel (forked processes, macOS/Linux).
#   Chains are saved per dataset; if a saved file exists it is read instead of
#   refitting, so the loop can be restarted where it stopped.

# settings and helpers ---------------------------------------------------------
nChains <- 4
mcmcSettings <- list(nitt = 260000, burnin = 60000, thin = 200)

prior <- list(
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))
)

# model data: proportions, species factors with the same levels in both columns
prep_model_data <- function(df) {
  d <- df %>%
    transmute(species1, species2, sisterSpecies,
              congeneric = factor(congeneric, levels = c("Different Genus", "Same Genus")),
              niche_overlap_percent = percentOverlap_4D / 100,
              range_overlap_percent = range_percentOverlap / 100) %>%
    filter(!is.na(range_overlap_percent))
  sp <- sort(unique(c(d$species1, d$species2)))
  d$species1 <- factor(d$species1, levels = sp)
  d$species2 <- factor(d$species2, levels = sp)
  d
}

fit_chains <- function(d) {
  parallel::mclapply(seq_len(nChains), function(i) {
    set.seed(i)
    MCMCglmm(range_overlap_percent ~ niche_overlap_percent * congeneric,
             random = ~ mm(species1 + species2),
             family = "gaussian", data = d, prior = prior,
             nitt = mcmcSettings$nitt, burnin = mcmcSettings$burnin,
             thin = mcmcSettings$thin, verbose = FALSE)
  }, mc.cores = nChains)
}

# pMCMC from posterior draws (same definition MCMCglmm uses)
p_from_draws <- function(x) 2 * max(0.5 / length(x), min(mean(x > 0), mean(x < 0)))
ci <- function(x) sprintf("%.2f [%.2f, %.2f]", mean(x), quantile(x, 0.025), quantile(x, 0.975))

# fit all 12 models ----------------------------------------------------(slooow)
modelChains <- list(); modelResults <- list()
for (k in names(nicheRangeOverlaps)) {
  d <- prep_model_data(nicheRangeOverlaps[[k]])
  f <- paste0("MixedEffectsModel_", k, "_4chains_20261124.rds")
  if (file.exists(f)) {
    chains <- read_rds(f)
  } else {
    chains <- fit_chains(d)
    write_rds(chains, f)
  }
  modelChains[[k]] <- chains

  # convergence and pooled posterior
  sol_list <- mcmc.list(lapply(chains, function(m) m$Sol))
  sol <- do.call(rbind, lapply(chains, function(m) m$Sol))
  rhat <- gelman.diag(sol_list)$psrf[, "Point est."]
  ess  <- effectiveSize(sol_list)
  slopeDiff <- sol[, "niche_overlap_percent:congenericSame Genus"]
  slopeHet  <- sol[, "niche_overlap_percent"]
  slopeCon  <- slopeHet + slopeDiff

  modelResults[[k]] <- tibble(
    dataset = k, buffer_m = nicheRangeOverlaps[[k]]$buffer_m[1], extent = nicheRangeOverlaps[[k]]$extent[1],
    nPairs = nrow(d), nCongeneric = sum(d$congeneric == "Same Genus"),
    slope_differentGenus = ci(slopeHet), slope_sameGenus = ci(slopeCon),
    slope_difference = ci(slopeDiff), pMCMC = p_from_draws(slopeDiff),
    max_Rhat = round(max(rhat), 3), min_ESS = round(min(ess)))
  message(k, " done")
}


modelSummary <- bind_rows(modelResults)
modelSummary
write_csv(modelSummary, "MixedEffectsModels_summary_20261124.csv") # supplementary table

# main model: convergence plots and figures --------------------------------------
figKey  <- "buffer150m_full"   # 150 m buffer, full breeding ranges
figData <- prep_model_data(nicheRangeOverlaps[[figKey]])
chains  <- modelChains[[figKey]]

sol_list <- mcmc.list(lapply(chains, function(m) m$Sol))   # fixed effects
gelman.diag(sol_list)      # R-hat: point est + upper CI (want ~1.00)
gelman.plot(sol_list)      # shows R-hat settling to 1 as iterations accrue
plot(sol_list)             # trace + density; want all 4 chains overlapping in one band
effectiveSize(sol_list)    # pooled ESS across chains (~4,000 if each was ~1,000)
autocorr.diag(sol_list)    # autocorrelation should drop off fast

# pooled posterior
sol_all <- do.call(rbind, lapply(chains, function(m) m$Sol))
vcv_all <- do.call(rbind, lapply(chains, function(m) m$VCV))
apply(sol_all, 2, mean)                         # posterior means
t(apply(sol_all, 2, quantile, c(.025, .975)))   # 95% credible intervals

# ---- Build fitted lines from the POOLED posterior (all 4 chains) ----
sol  <- sol_all
lev  <- levels(figData$congeneric)          # lev[1] = reference (heterogeneric)
xseq <- seq(min(figData$niche_overlap_percent),
            max(figData$niche_overlap_percent), length.out = 100)

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

# slope labels (posterior mean + 95% credible interval), placed near the top
yTop <- max(pred_df$upper)
lab_df <- do.call(rbind, lapply(lev, function(g) {
  b1 <- coefs(g)$b1
  data.frame(congeneric = factor(g, levels = lev),
             label = sprintf("slope = %.2f [%.2f, %.2f]",
                             mean(b1), quantile(b1, 0.025), quantile(b1, 0.975)))
}))
lab_df$x <- min(xseq)
lab_df$y <- yTop * c(1.15, 1.08)[seq_along(lev)]
eff_df <- data.frame(x = min(xseq), y = yTop * 1.01, label = eff_label)

# ---- Plot it ----
# Figure 8
ggplot(figData,
       aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_ribbon(data = pred_df,
              aes(y = fit, ymin = lower, ymax = upper, fill = congeneric),
              alpha = 0.2, color = NA) +
  geom_line(data = pred_df, aes(y = fit), linewidth = 1) +
  geom_point(data = subset(figData, sisterSpecies == TRUE), color = "black") +
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
ggsave(paste0("Figure8_NicheRangeOverlap_", figKey, "_20261124.pdf"), height = 5, width = 7)

# supp figure: raw pairs with the fitted lines, one panel per group
ggplot(figData,
       aes(x = niche_overlap_percent, y = range_overlap_percent, color = congeneric)) +
  geom_ribbon(data = pred_df,
              aes(y = fit, ymin = lower, ymax = upper, fill = congeneric),
              alpha = 0.2, color = NA) +
  geom_point(alpha = 0.2) +
  geom_line(data = pred_df, aes(y = fit), linewidth = 1) +
  geom_text(data = lab_df, aes(x = x, y = max(figData$range_overlap_percent) * 1.05, label = label),
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
ggsave(paste0("FigureS_NicheRangeOverlap_points_", figKey, "_20261124.pdf"), height = 5, width = 7)

# congeneric pairs only, sister species in black
ggplot(subset(figData, congeneric == "Same Genus"),
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
# 22) Leave-one-genus-out models (main dataset)
################################################################################
# 24 Nov 2026. Reviewer request: is the congeneric vs heterogeneric difference
# driven by a single genus? Main dataset only (150 m buffer, full breeding
# ranges). For each genus with at least one congeneric pair, all of its species
# and every pair involving them are removed and the section 21 model is refitted.
# Uses prep_model_data(), fit_chains(), prior, mcmcSettings, p_from_draws() and
# ci() from section 21. Chains are saved per genus and re-read if present, so
# the loop can be restarted where it stopped.

logoKey  <- "buffer150m_full"
logoData <- nicheRangeOverlaps[[logoKey]]
congenericGenera <- sort(unique(logoData$genus1[logoData$congeneric == "Same Genus"]))
length(congenericGenera) # 16 genera

# slope summaries from a list of chains (numbers for the figure, text for the table)
summarise_chains <- function(chains) {
  sol_list <- mcmc.list(lapply(chains, function(m) m$Sol))
  sol <- do.call(rbind, lapply(chains, function(m) m$Sol))
  slopeDiff <- sol[, "niche_overlap_percent:congenericSame Genus"]
  slopeHet  <- sol[, "niche_overlap_percent"]
  slopeCon  <- slopeHet + slopeDiff
  tibble(slope_differentGenus = ci(slopeHet), slope_sameGenus = ci(slopeCon),
         slope_difference = ci(slopeDiff), pMCMC = p_from_draws(slopeDiff),
         diff_mean = mean(slopeDiff), diff_lower = quantile(slopeDiff, 0.025),
         diff_upper = quantile(slopeDiff, 0.975),
         max_Rhat = round(max(gelman.diag(sol_list)$psrf[, "Point est."]), 3),
         min_ESS = round(min(effectiveSize(sol_list))))
}

logoResults <- list()
for (g in congenericGenera) {
  d <- logoData %>% filter(genus1 != g, genus2 != g) %>% prep_model_data()
  f <- paste0("MixedEffectsModel_LOGO_", g, "_", logoKey, "_4chains_20261124.rds")
  if (file.exists(f)) {
    chains <- read_rds(f)
  } else {
    chains <- fit_chains(d)
    write_rds(chains, f)
  }
  logoResults[[g]] <- bind_cols(
    tibble(genusRemoved = g, speciesRemoved = sum(genus[analysisSpecies] == g),
           nPairs = nrow(d), nCongeneric = sum(d$congeneric == "Same Genus")),
    summarise_chains(chains))
  message(g, " done")
}

# full model (section 21) for comparison
fullChains <- if (exists("modelChains") && !is.null(modelChains[[logoKey]])) modelChains[[logoKey]] else
  read_rds(paste0("MixedEffectsModel_", logoKey, "_4chains_20261124.rds"))
fullModel <- bind_cols(
  tibble(genusRemoved = "none (full model)", speciesRemoved = 0,
         nPairs = nrow(prep_model_data(logoData)),
         nCongeneric = sum(logoData$congeneric == "Same Genus")),
  summarise_chains(fullChains))

logoSummary <- bind_rows(fullModel, bind_rows(logoResults))
logoSummary %>% dplyr::select(genusRemoved, nCongeneric, slope_difference, pMCMC, max_Rhat)
write_csv(logoSummary, "MixedEffectsModels_leaveOneGenusOut_20261124.csv") # supplementary table

# supplemental figure: slope difference with each genus removed -------------------
logoPlot <- bind_rows(logoResults) %>%
  mutate(genusRemoved = factor(genusRemoved, levels = rev(congenericGenera)))
ggplot(logoPlot, aes(x = diff_mean, y = genusRemoved)) +
  annotate("rect", xmin = fullModel$diff_lower, xmax = fullModel$diff_upper,
           ymin = -Inf, ymax = Inf, fill = "grey85") +
  geom_vline(xintercept = fullModel$diff_mean, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_pointrange(aes(xmin = diff_lower, xmax = diff_upper), colour = "darkgreen") +
  labs(x = "Difference in slopes, congeneric minus heterogeneric (95% CI)",
       y = "Genus removed") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))
ggsave("FigureS_leaveOneGenusOut_20261124.pdf", height = 6, width = 7)


################################################################################
# 23) permutation test
################################################################################
# 24 Nov 2026. Could the congeneric vs heterogeneric slope difference arise by
# chance from how species are grouped into genera? Genus labels are shuffled
# across species (keeping the number of genera and their sizes), congeneric
# status is recalculated for every pair, and the interaction (slope difference)
# is refitted with a linear model. Main dataset (150 m buffer, full breeding
# ranges), same variables as section 21: 4-D niche overlap and range overlap as
# proportions. Uses nicheRangeOverlaps, analysisSpecies and genus from section 20.

permKey  <- "buffer150m_full"
permData <- nicheRangeOverlaps[[permKey]] %>%
  transmute(species1, species2,
            niche_overlap_percent = percentOverlap_4D / 100,
            range_overlap_percent = range_percentOverlap / 100) %>%
  filter(!is.na(range_overlap_percent))

# ---- 1. Species -> genus map (section 20) ----
genus_of <- genus[analysisSpecies]
stopifnot(!anyNA(genus_of))

# ---- 2. Statistic: interaction (slope difference) for a given species -> genus map ----
interaction_coef <- function(gmap) {
  cong_status <- gmap[permData$species1] == gmap[permData$species2]
  fit <- lm(range_overlap_percent ~ niche_overlap_percent * cong_status, data = permData)
  unname(coef(fit)["niche_overlap_percent:cong_statusTRUE"])
}

obs <- interaction_coef(genus_of)                     # observed slope difference
sum(genus_of[permData$species1] == genus_of[permData$species2]) # 37 congeneric pairs

# ---- 3. Null: shuffle genus labels across species (preserves genus-size distribution) ----
set.seed(1)
nperm   <- 5000
null <- replicate(nperm, {
  permuted <- setNames(sample(unname(genus_of)), names(genus_of))
  interaction_coef(permuted)
})

# ---- 4. Two-sided permutation p-value ----
p_perm <- (sum(abs(null) >= abs(obs)) + 1) / (nperm + 1)
c(observed_slope_diff = obs, p_perm = p_perm)

write_rds(list(dataset = permKey, observed = obs, null = null, p = p_perm, nperm = nperm),
          paste0("permutationTest_", permKey, "_20261124.rds"))

# ---- 5. Supplemental figure ----
ggplot(data.frame(null = null), aes(x = null)) +
  geom_histogram(bins = 40, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = obs, colour = "red", linewidth = 1) +
  labs(title = "Permutation null: slope difference",
       subtitle = paste0("Observed = ", round(obs, 2), ", p = ", signif(p_perm, 2),
                         " (", format(nperm, big.mark = ","), " permutations)"),
       x = "Interaction coefficient (difference in slopes)", y = "Frequency") +
  theme_minimal()
ggsave(paste0("FigureS_permutationTest_", permKey, "_20261124.pdf"), height = 5, width = 6)


################################################################################
# 24) Figure 4 acoustic-use space volumes
################################################################################
# 24 Nov 2026. 4-D acoustic-use space volume (PC1-PC4, section 15) of each of
# the 89 analysis species, sorted by size. Main dataset (150 m buffer). Uses
# speciesAreasBuffers (section 17) and analysisSpecies (section 20).

b <- "150"
volumes <- speciesAreasBuffers[[b]] %>%
  filter(species %in% analysisSpecies)
stopifnot(nrow(volumes) == length(analysisSpecies))

ggplot(volumes, aes(x = reorder(species, volume4D), y = volume4D)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip the coordinates to make it a horizontal bar chart
  labs(
    x = "Species",
    y = "Acoustic-use space volume (PC units^4)",
    title = "Species acoustic-use space volumes sorted by size"
  ) +
  theme_minimal()

ggsave(paste0("Figure4_acousticUseSpaceVolumes_buffer", b, "m_20261124.pdf"), height = 12, width = 7)


################################################################################
# 25) Figure 5 - matrix
################################################################################
# 24 Nov 2026. Hierarchically clustered matrix of pairwise 4-D acoustic-use
# space overlap (Jaccard, %, section 15) for the 89 analysis species, with a
# bar plot of each species' 4-D volume. Main dataset (150 m buffer). Uses
# pairwiseOverlapsBuffers and speciesAreasBuffers (section 17), siteDataBuffers
# and traitsPCABuffers (section 12), analysisSpecies (section 20).
# Bars are coloured by each species' detection-weighted mean PC1 score (the
# centre of its acoustic-use space on PC1), replacing the old 2-D polygon
# centroid. Colour legends are drawn under the figure (previously added by
# hand). The land cover types on the PC1 legend sit at the mean PC1 score of
# the sites of that type (previously placed by hand from Figure 3).

b <- "150"

# detection-weighted mean PC1 per species --------------------------------------
pc1   <- traitsPCABuffers[[b]]$x[, 1]
rates <- siteDataBuffers[[b]][, analysisSpecies]
stopifnot(length(pc1) == nrow(rates))
speciesPC1 <- tibble(species = analysisSpecies,
                     PC1 = sapply(analysisSpecies, function(sp) weighted.mean(pc1, rates[[sp]])))

# land cover positions on PC1 (mean PC1 of the sites of each type) -------------
legendCover <- c("Mediterranean California Red Fir Forest",
                 "California Montane Jeffrey Pine-(Ponderosa Pine) Woodland",
                 "Mediterranean California Mesic Mixed Conifer Forest and Woodland",
                 "Mediterranean California Mixed Oak Woodland",
                 "California Central Valley Mixed Oak Savanna")
coverPC1 <- tibble(landcover = siteDataBuffers[[b]]$landcover_type, PC1 = pc1) %>%
  filter(landcover %in% legendCover) %>%
  group_by(landcover) %>%
  dplyr::summarise(nSites = n(), PC1 = mean(PC1))
stopifnot(nrow(coverPC1) == length(legendCover))
coverPC1

# one colour range for the bars and the legend, wide enough for both
pc1Lims <- range(c(speciesPC1$PC1, coverPC1$PC1))

# function: vertical colour bar legend, optional text labels to the right ------
colour_legend <- function(lims, option, title, breaks, labels = NULL, textSize = 3.5) {
  y <- seq(lims[1], lims[2], length.out = 200)
  g <- ggplot() +
    geom_tile(data = tibble(y = y), aes(x = 0.5, y = y, fill = y),
              width = 1, height = diff(y[1:2]) * 1.02) +
    scale_fill_viridis_c(option = option, limits = lims, guide = "none") +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, if (is.null(labels)) 3 else 10), expand = c(0, 0)) +
    labs(title = title) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 10),
          axis.ticks.y = element_line(), plot.title = element_text(size = 12))
  if (!is.null(labels))
    g <- g + ggrepel::geom_text_repel(data = labels, aes(x = 1, y = PC1, label = landcover),
                                      hjust = 0, direction = "y", nudge_x = 0.4,
                                      xlim = c(1.2, NA), size = textSize,
                                      segment.size = 0.2, min.segment.length = 0)
  g
}

# function: heatmap + bar plot for a set of species ----------------------------
# Clustering is redone for each set of species.
overlap_matrix_plot <- function(spp, b = "150", barWidth = 1, legendHeight = 0.2) {
  pairs <- pairwiseOverlapsBuffers[[b]] %>%
    filter(species1 %in% spp, species2 %in% spp) %>%
    dplyr::select(species1, species2, overlap = percentOverlap_4D)
  stopifnot(nrow(pairs) == choose(length(spp), 2))

  # both directions plus the diagonal (100% with itself)
  allOverlaps <- bind_rows(pairs,
                           rename(pairs, species1 = species2, species2 = species1),
                           tibble(species1 = spp, species2 = spp, overlap = 100))
  overlapMatrix <- acast(allOverlaps, species1 ~ species2, value.var = "overlap")

  # hierarchical clustering on 1 - overlap
  hc <- hclust(as.dist(1 - overlapMatrix / 100), method = "complete")
  speciesOrder <- hc$labels[hc$order]
  allOverlaps <- allOverlaps %>%
    mutate(species1 = factor(species1, levels = speciesOrder),
           species2 = factor(species2, levels = speciesOrder))

  vols <- speciesAreasBuffers[[b]] %>%
    filter(species %in% spp) %>%
    left_join(speciesPC1, by = "species") %>%
    mutate(species = factor(species, levels = speciesOrder))

  heatmapPlot <- ggplot(allOverlaps, aes(x = species1, y = species2, fill = overlap)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(speciesOrder)) +
    scale_fill_viridis_c(name = "Overlap (%)", limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
          legend.position = "none") +
    xlab("") + ylab("")

  barPlot <- ggplot(vols, aes(y = species, x = volume4D, fill = PC1)) +
    geom_col() +
    scale_y_discrete(limits = rev(speciesOrder)) +
    xlab("Acoustic-use space\n(4D volume)") +
    scale_fill_viridis_c(option = "A", limits = pc1Lims) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          legend.position = "none")

  # legends under the figure
  overlapLegend <- colour_legend(c(0, 100), "D", "Percent overlap of\nacoustic-use space",
                                 breaks = seq(0, 100, 25))
  pc1Legend <- colour_legend(pc1Lims, "A", "Species location along PC1\n(detection-weighted mean)",
                             breaks = pretty(pc1Lims), labels = coverPC1)
  legends <- plot_grid(overlapLegend, pc1Legend, NULL, nrow = 1, rel_widths = c(1.6, 4.5, 2))

  top <- plot_grid(heatmapPlot, barPlot, align = "h", axis = "tblr",
                   rel_widths = c(8, barWidth))
  list(heatmap = heatmapPlot, bar = barPlot, order = speciesOrder, legends = legends,
       combined = plot_grid(top, legends, ncol = 1, rel_heights = c(1, legendHeight)))
}

# all 89 species ---------------------------------------------------------------
matrixAll <- overlap_matrix_plot(analysisSpecies, b)
matrixAll$combined
ggsave(paste0("Figure5_matrix_4D_buffer", b, "m_20261124.pdf"), matrixAll$combined,
       height = 18, width = 15)

# a reviewer asked for a simplified version with fewer species -----------------
# Same species as before, minus Woodhouse's Scrub-Jay (excluded, section 20);
# the two sapsuckers are now the complex.
spSubset <- c("Pacific Wren", "Brown Creeper", "Red-breasted Nuthatch",
              "Red-breasted/Red-naped Sapsucker", "White-headed Woodpecker",
              "Townsend's Solitaire", "Yellow-rumped Warbler", "Hammond's Flycatcher",
              "Hermit Warbler", "Hermit Thrush", "Golden-crowned Kinglet",
              "Mountain Chickadee", "Olive-sided Flycatcher", "Cassin's Finch",
              "Red Crossbill", "Green-tailed Towhee", "Dusky Flycatcher",
              "California Quail", "California Towhee", "Oak Titmouse",
              "Blue-gray Gnatcatcher", "Bewick's Wren", "Wrentit", "Nuttall's Woodpecker",
              "Mourning Dove", "California Scrub-Jay", "Black Phoebe", "Acorn Woodpecker",
              "Bushtit", "Black-headed Grosbeak", "Plumbeous Vireo", "Purple Finch",
              "Cassin's Vireo", "Nashville Warbler", "Common Raven", "American Robin",
              "Hairy Woodpecker", "Northern Flicker", "Steller's Jay", "Dark-eyed Junco",
              "Western Tanager", "Western Wood-Pewee")
missingSubset <- setdiff(spSubset, analysisSpecies)
if (length(missingSubset) > 0) stop("Not in analysisSpecies: ", paste(missingSubset, collapse = ", "))

matrixSub <- overlap_matrix_plot(spSubset, b, barWidth = 2, legendHeight = 0.35)
matrixSub$combined
ggsave(paste0("Figure5_matrix_4D_subset_buffer", b, "m_20261124.pdf"), matrixSub$combined,
       height = 11, width = 11)


################################################################################
# 26) Species details supp table
################################################################################
# 24 Nov 2026. One row per species (the 92 species kept in section 9, 150 m
# buffer), in the format of Table S1 (S1_speciesDetails in the Draft 7
# Supplemental Tables.xlsx), with updated values:
#   - totalDetections, totalSites: speciesCounts[["150"]] (section 9)
#   - area_PC12, area_PC34: 2-D acoustic-use space areas, 150 m buffer
#   - volume4D_50m, volume4D_150m, volume4D_300m: 4-D acoustic-use space
#     volume for each buffer (speciesAreasBuffers, section 17)
#   - genus (from the scientific name), family (AVONET eBird taxonomy,
#     AnimalTraits/AVONET2_eBird.csv, matched by genus)
#   - nCongeners: number of other analysis species in the same genus
#   - rangeArea_km2_*: breeding range area in each extent (rangeAreas, section 18)
#   - crossID_pctShared: % of the species' BirdNET detections shared with a
#     congener in the same 3 s window (crossID, section 2); NA = no congener detected
#   - inAnalysis: FALSE for Woodhouse's Scrub-Jay and the two sapsuckers
#     (excludedSpecies, section 20)
# BirdNET validation columns (n_samples to cutoff90.l_recall, added by Connor
# 06 July 2026) are copied from the Draft 7 table. That table used dot names
# (Bewick.s.Wren), so names are converted only to match; the output keeps the
# normal names. Downy Woodpecker (dropped before because the old 2-D KDE failed
# for it) and the sapsucker complex were not in that table, so they are NA
# (the complex members are validated separately).

b <- "150"
S1old <- readxl::read_excel(file.path(base_dir, "Draft 7_OikosRevision_20260629",
                                      "Supplemental Tables.xlsx"),
                            sheet = "S1_speciesDetails")
validationCols <- c("n_samples", "intercept.l", "beta.l", "cutoff90.l",
                    "cutoff90.l_precision", "cutoff90.l_recall")
dot_name <- function(x) gsub("[^A-Za-z]", ".", x)  # old table used dot names
validation <- S1old[c("species", validationCols)] %>% rename(dotName = species)

# family by genus (AVONET, eBird taxonomy)
avonet <- read_csv(file.path(base_dir, "AnimalTraits", "AVONET2_eBird.csv"),
                   locale = locale(encoding = "latin1"), show_col_types = FALSE)
families <- avonet %>%
  transmute(genus = word(Species2, 1), family = Family2) %>%
  distinct()
# (a genus listed under two families would duplicate rows; the check below catches it)

# 4-D volume for each buffer, wide
volumesWide <- bind_rows(lapply(names(speciesAreasBuffers), function(bb)
  speciesAreasBuffers[[bb]] %>% transmute(species, buffer = paste0("volume4D_", bb, "m"), volume4D))) %>%
  pivot_wider(names_from = buffer, values_from = volume4D)

speciesDetails <- speciesCounts[[b]] %>%
  filter(keep) %>%
  transmute(species = common_name, totalDetections, totalSites) %>%
  left_join(dplyr::select(commonNames, species = common_name, scientific_name),
            by = "species") %>%
  mutate(inAnalysis = !species %in% excludedSpecies,
         genus = word(scientific_name, 1)) %>%
  left_join(families, by = "genus") %>%
  group_by(genus) %>%
  mutate(nCongeners = sapply(species, function(s) sum(species[inAnalysis] != s))) %>%
  ungroup() %>%
  left_join(dplyr::select(speciesAreasBuffers[[b]], species, area_PC12, area_PC34),
            by = "species") %>%
  left_join(volumesWide, by = "species") %>%
  left_join(rename_with(rangeAreas, ~ sub("^area_km2_", "rangeArea_km2_", .x)),
            by = "species") %>%
  left_join(dplyr::select(crossID, species = common_name, crossID_pctShared = pct_shared),
            by = "species") %>%
  mutate(dotName = dot_name(species)) %>%
  left_join(validation, by = "dotName") %>%
  dplyr::select(species, scientific_name, genus, family, nCongeners, inAnalysis,
                area_PC12, area_PC34, starts_with("volume4D_"),
                totalDetections, totalSites, starts_with("rangeArea_km2_"),
                crossID_pctShared, all_of(validationCols)) %>%
  arrange(species)

# checks
stopifnot(nrow(speciesDetails) == n_distinct(speciesDetails$species),
          sum(speciesDetails$inAnalysis) == length(analysisSpecies),
          !anyNA(dplyr::select(speciesDetails, family, starts_with("volume4D_"))))
setdiff(validation$dotName, dot_name(speciesDetails$species)) # old rows not matched (should be none)
speciesDetails$species[is.na(speciesDetails$n_samples)]       # no validation data

write_rds(speciesDetails, "speciesDetails_20261124.rds")
write_csv(speciesDetails, "speciesDetails_20261124.csv")

# RUNNING LEFT OFF HERE ######
