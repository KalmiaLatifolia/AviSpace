
# resampling.R
# Resample the fine-resolution layers (foliar traits and canopy height, 30 m) to
# the BioCube grid (30 arc-seconds, ~1 km, WGS84), so that every PCA variable in
# fullProtocol.R has the same grain ("equal-grain" version of the analysis).
# Laura Berman
# created 09 Nov 2026

# Outputs, written to California_biocube_1km_layers/ (read by fullProtocol section 7):
#   CA_Physiology_<trait>_mean_1km.tif  mean of the 30 m pixels in each 1 km cell
#   CA_Physiology_<trait>_std_1km.tif   standard deviation of those 30 m pixels
#   Favrichon_Sierra_<layer>_1km.tif    mean canopy height (and RH50)
# Quality layers, written to resampling_QC/ (NOT read by fullProtocol):
#   <trait>_validShare_1km.tif  share of 30 m pixels in each 1 km cell with a trait value

library(terra)

base_dir <- "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics"
in_dir  <- file.path(base_dir, "native_resolution_layers")
out_dir <- file.path(base_dir, "California_biocube_1km_layers")
qc_dir  <- file.path(base_dir, "resampling_QC")
dir.create(qc_dir, showWarnings = FALSE)

# template: an existing BioCube layer, so new layers match its grid, extent and CRS exactly
template <- rast(file.path(out_dir, "CA_Climate_CHELSA_bio01_v2.tif"))
template


# function: 30 m layer(s) -> 1 km mean, sd and share of valid pixels -----------
# How the SD is calculated: for each 1 km cell, we take the average of the 30 m
# values (mean) and the average of the squared 30 m values. The variance of the
# 30 m pixels in that cell is then (average of squares) - (mean squared), and the
# SD is its square root. This is the usual SD of all valid 30 m pixels in the
# cell (population SD), computed with terra's "average" resampling because terra
# has no direct "sd" resampling option.

to_1km <- function(files, zero_is_na = TRUE) {
  parts <- lapply(files, function(f) {
    r <- rast(f)[[1]] # first band only. The trait mosaics have 12 bands
                      # (s6v3_mean, s6v3_sd, QC flags, flight_id); band 1 is the
                      # trait value, the same band raster() reads in fullProtocol
    if (zero_is_na) r <- classify(r, cbind(0, NA)) # 0 = no data in the trait maps
    valid <- !is.na(r)                             # 1 = has a value, 0 = missing
    s <- c(r, r^2, valid)
    names(s) <- c("x", "x2", "valid")
    # reproject from UTM 11N and average all 30 m pixels falling in each 1 km cell
    # (missing pixels are ignored in the averages of x and x2)
    project(s, template, method = "average")
  })
  out <- if (length(parts) == 1) parts[[1]] else do.call(mosaic, c(parts, fun = "mean"))
  sd <- sqrt(max(out$x2 - out$x^2, 0)) # max(..., 0) guards against tiny negative rounding errors
  names(sd) <- "sd"
  list(mean = out$x, sd = sd, validShare = out$valid)
}


# foliar traits (Tahoe and Yosemite mosaics, 30 m) ------------------------------(slow)

traitNames <- c("Nitrogen", "LMA", "Cellulose", "Calcium", "Chlorophylls",
                "Fiber", "Lignin", "NSC", "Phenolics", "Phosphorus",
                "Potassium", "Starch", "Sugar", "Sulfur")

for (trait in traitNames) {
  files <- file.path(in_dir, c(paste0("tahoe_20180621_", trait, "_30m_polished.tif"),
                               paste0("yosemite_20180622_", trait, "_30m_polished.tif")))
  res <- to_1km(files, zero_is_na = TRUE)
  writeRaster(res$mean, file.path(out_dir, paste0("CA_Physiology_", trait, "_mean_1km.tif")), overwrite = TRUE)
  writeRaster(res$sd, file.path(out_dir, paste0("CA_Physiology_", trait, "_std_1km.tif")), overwrite = TRUE)
  writeRaster(res$validShare, file.path(qc_dir, paste0(trait, "_validShare_1km.tif")), overwrite = TRUE)
  message(trait, " done")
}

# check
x <- rast(file.path(out_dir, "CA_Physiology_Nitrogen_mean_1km.tif"))
compareGeom(x, template) # TRUE = same grid, extent and CRS as the BioCube layers
plot(c(x, rast(file.path(out_dir, "CA_Physiology_Nitrogen_std_1km.tif")),
       rast(file.path(qc_dir, "Nitrogen_validShare_1km.tif"))))


# canopy height (Favrichon, 30 m) ---------------------------------------------
# 0 m is a real value here (no canopy), so zeros are kept; -9999 is the nodata flag

for (layer in c("CanopyHeight", "RH50")) {
  res <- to_1km(file.path(in_dir, paste0("Favrichon_Sierra_", layer, "_30m_20250508.tif")),
                zero_is_na = FALSE)
  writeRaster(res$mean, file.path(out_dir, paste0("Favrichon_Sierra_", layer, "_1km.tif")), overwrite = TRUE)
  writeRaster(res$validShare, file.path(qc_dir, paste0(layer, "_validShare_1km.tif")), overwrite = TRUE)
  message(layer, " done")
}


# add the new file names to tidyNames.csv (used by fullProtocol section 7) ------
# only rows that are not already there are added

tidy_file <- file.path(base_dir, "Draft 1", "tidyNames.csv")
tidy <- read.csv(tidy_file, fileEncoding = "UTF-8-BOM")
new <- rbind(
  data.frame(VariableName = traitNames,
             FileName = paste0("CA_Physiology_", traitNames, "_mean_1km.tif"),
             Category = "Physiology"),
  data.frame(VariableName = paste0(traitNames, "_std"),
             FileName = paste0("CA_Physiology_", traitNames, "_std_1km.tif"),
             Category = "Physiology"),
  data.frame(VariableName = c("CanopyHeight", "MidCanopyHeight"),
             FileName = c("Favrichon_Sierra_CanopyHeight_1km.tif", "Favrichon_Sierra_RH50_1km.tif"),
             Category = "Function"))
new <- new[!new$FileName %in% trimws(tidy$FileName), ]
if (nrow(new) > 0) {
  cat(paste0("\r\n", paste(new$VariableName, new$FileName, new$Category, sep = ","), collapse = ""),
      file = tidy_file, append = TRUE)
}
nrow(new) # rows added (30 the first time, 0 after that)
