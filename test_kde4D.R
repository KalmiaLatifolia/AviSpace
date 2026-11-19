################################################################################
# TEST: weighted 4-D kernel density (ks::kde) acoustic-use space
################################################################################
# 19 Nov 2026. Stand-alone test, not part of fullProtocol.R.
# Run AFTER fullProtocol.R sections 1-14 so these objects are in memory:
#   siteDetections_foliarTraits_BioCubeBuffers, traitsPCABuffers, speciesCols
# Writes test figures/tables (prefix "test_kde4D_") to the working directory.
#
# Method (matches the 2-D method in sections 16-19 as closely as possible):
#   - space   = site scores on PC1-PC4 (raw scores, not rescaled)
#   - weights = each species' unrounded detection rate at each site
#               (no rounding, no uncount, no jitter needed)
#   - bandwidth = one diagonal plug-in bandwidth, estimated from all sites and
#               shared by every species, so volumes and overlaps are comparable
#   - grid    = one common grid for every species; cell width on each PC is a
#               fixed fraction of that PC's bandwidth
#   - space use = smallest set of grid cells holding 90% of the density
#               (same 90% threshold as section 16)
#   - overlap = Jaccard: shared cells / union of cells (same as section 19)
# The 2-D (PC1 x PC2) version is computed the same way for comparison.

library(ks)
library(tidyverse)

b         <- "150"   # buffer version to test
nPC       <- 4
prob      <- 0.90    # share of density inside the space-use region
cellFrac4D <- 0.75   # 4-D grid cell width, as a fraction of that PC's bandwidth
cellFrac2D <- 0.25   # 2-D grid cell width (2-D grids are cheap, so finer)
# Grid spacing is set per axis from the bandwidth, so every PC is resolved at
# the same scale relative to how much it is smoothed. PCs with more spread
# (PC1) get wider bandwidths and more grid points. Cells narrower than ~0.7 x
# bandwidth in 4-D run out of memory (>8 GB) on a laptop.
testSpecies <- c("Hermit Thrush", "Swainson's Thrush", "Yellow-rumped Warbler",
                 "Red-breasted Nuthatch", "Mountain Chickadee", "Spotted Towhee")

# data -------------------------------------------------------------------------
dat <- siteDetections_foliarTraits_BioCubeBuffers[[b]]
dat$geometry <- NULL
dat <- as.data.frame(dat)
pcs <- traitsPCABuffers[[b]]$x[, 1:nPC]
stopifnot(nrow(pcs) == nrow(dat))  # PCA rows are the same sites, same order

missingSp <- setdiff(testSpecies, speciesCols)
if (length(missingSp) > 0) stop("Not in speciesCols: ", paste(missingSp, collapse = ", "))
rates <- dat[, testSpecies]
colSums(rates > 0) # sites with detections per test species

# variance explained by the PCs used
round(100 * summary(traitsPCABuffers[[b]])$importance["Proportion of Variance", 1:nPC], 1)

# common bandwidth and grid ----------------------------------------------------
H4 <- Hpi.diag(pcs)          # diagonal: one bandwidth per PC, scaled to its spread
H2  <- Hpi.diag(pcs[, 1:2])
H34 <- Hpi.diag(pcs[, 3:4])  # only for the 2-D outline on the PC3 x PC4 figure
sqrt(diag(H4))               # bandwidth (SD units of each PC)

pad  <- 3 * sqrt(diag(H4))
xmin <- apply(pcs, 2, min) - pad
xmax <- apply(pcs, 2, max) + pad

# grid points per axis = padded range / (cellFrac x bandwidth)
grid_points <- function(H, xmin, xmax, cellFrac) ceiling((xmax - xmin) / (cellFrac * sqrt(diag(H)))) + 1
grid4D  <- grid_points(H4,  xmin,      xmax,      cellFrac4D)
grid2D  <- grid_points(H2,  xmin[1:2], xmax[1:2], cellFrac2D)
grid34  <- grid_points(H34, xmin[3:4], xmax[3:4], cellFrac2D)
grid4D; prod(grid4D)         # points per PC; total 4-D cells
round((xmax - xmin) / (grid4D - 1), 2)  # cell width per PC (PC units)

# function: weighted KDE + 90% region on the common grid -----------------------
kde_region <- function(X, w, H, xmin, xmax, gridsize, prob = 0.90) {
  w  <- w * nrow(X) / sum(w)  # ks needs weights summing to n
  fh <- kde(X, H = H, w = w, xmin = xmin, xmax = xmax,
            gridsize = gridsize, binned = TRUE, compute.cont = FALSE)
  est <- fh$estimate
  est[est < 0] <- 0           # binned estimates can have tiny negative values
  cellVol <- prod(sapply(fh$eval.points, function(g) diff(g[1:2])))
  z   <- sort(as.vector(est), decreasing = TRUE)
  thr <- z[which(cumsum(z) / sum(z) >= prob)[1]]
  inside <- est >= thr
  list(eval.points = fh$eval.points, estimate = est, inside = inside,
       volume = sum(inside) * cellVol,
       massOnGrid = sum(est) * cellVol)  # should be ~1; <1 means density off the grid
}

# run 4-D and 2-D for each test species ----------------------------------------
fits4D <- list(); fits2D <- list(); fits34 <- list(); timing <- c()
for (sp in testSpecies) {
  t0 <- Sys.time()
  fits4D[[sp]] <- kde_region(pcs, rates[[sp]], H4, xmin, xmax, grid4D, prob)
  timing[sp] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fits2D[[sp]] <- kde_region(pcs[, 1:2], rates[[sp]], H2, xmin[1:2], xmax[1:2], grid2D, prob)
  fits34[[sp]] <- kde_region(pcs[, 3:4], rates[[sp]], H34, xmin[3:4], xmax[3:4], grid34, prob)
}
round(timing, 1) # seconds per species in 4-D; x 92 species for the full run

volumes <- tibble(species   = testSpecies,
                  volume4D  = sapply(fits4D, `[[`, "volume"),
                  area2D    = sapply(fits2D, `[[`, "volume"),
                  mass4D    = sapply(fits4D, `[[`, "massOnGrid"),
                  mass2D    = sapply(fits2D, `[[`, "massOnGrid"))
volumes

# pairwise Jaccard overlap -----------------------------------------------------
jaccard <- function(a, b) sum(a & b) / sum(a | b)
overlaps <- combn(testSpecies, 2, simplify = FALSE) %>%
  map_dfr(~ tibble(species1 = .x[1], species2 = .x[2],
                   overlap4D = jaccard(fits4D[[.x[1]]]$inside, fits4D[[.x[2]]]$inside),
                   overlap2D = jaccard(fits2D[[.x[1]]]$inside, fits2D[[.x[2]]]$inside))) %>%
  arrange(desc(overlap4D))
overlaps

write_csv(volumes,  "test_kde4D_volumes.csv")
write_csv(overlaps, "test_kde4D_overlaps.csv")

# figure 1: 4-D region projected onto PC1 x PC2, with 2-D region outline -------
# fill = share of the PC3 x PC4 grid that is inside the 4-D region at that
# PC1/PC2 location. Red line = the 2-D 90% region (current method, unrounded).
project_region <- function(fit, keep) {
  m <- apply(fit$inside, keep, mean)
  expand.grid(x = fit$eval.points[[keep[1]]], y = fit$eval.points[[keep[2]]]) %>%
    mutate(share = as.vector(m))
}

proj12 <- map_dfr(testSpecies, ~ project_region(fits4D[[.x]], c(1, 2)) %>% mutate(species = .x))
outline_df <- function(fits) map_dfr(names(fits), ~ expand.grid(x = fits[[.x]]$eval.points[[1]],
                                                                y = fits[[.x]]$eval.points[[2]]) %>%
                       mutate(inside = as.numeric(as.vector(fits[[.x]]$inside)), species = .x))
outline2D  <- outline_df(fits2D)
outline34  <- outline_df(fits34)
# sites: grey dots = species never detected there; black circles = detected,
# area proportional to detection rate
sitePts <- map_dfr(testSpecies, ~ tibble(PC1 = pcs[, 1], PC2 = pcs[, 2], PC3 = pcs[, 3],
                                         PC4 = pcs[, 4], rate = rates[[.x]], species = .x))
siteZero <- filter(sitePts, rate == 0)
siteDet  <- filter(sitePts, rate > 0)

p1 <- ggplot() +
  geom_raster(data = filter(proj12, share > 0), aes(x, y, fill = share)) +
  geom_contour(data = outline2D, aes(x, y, z = inside), breaks = 0.5,
               colour = "red", linewidth = 0.5) +
  geom_point(data = siteZero, aes(PC1, PC2), colour = "grey60", size = 0.3) +
  geom_point(data = siteDet, aes(PC1, PC2, size = rate), shape = 21,
             colour = "black", fill = NA, stroke = 0.25) +
  scale_fill_viridis_c(name = "Share of PC3 x PC4\ninside 4-D region") +
  scale_size_area(max_size = 3, name = "Detection rate") +
  facet_wrap(~ species) + coord_equal() +
  labs(x = "PC1", y = "PC2",
       title = paste0("4-D acoustic-use space projected on PC1 x PC2 (", b, " m buffer)"),
       subtitle = "Red outline = 2-D 90% region from PC1 and PC2 only. Grey dots = sites with no detections") +
  theme_bw()
ggsave("test_kde4D_fig1_projection_PC12.pdf", p1, width = 12, height = 8)

# figure 2: same 4-D regions projected onto PC3 x PC4 --------------------------
proj34 <- map_dfr(testSpecies, ~ project_region(fits4D[[.x]], c(3, 4)) %>% mutate(species = .x))
p2 <- ggplot() +
  geom_raster(data = filter(proj34, share > 0), aes(x, y, fill = share)) +
  geom_contour(data = outline34, aes(x, y, z = inside), breaks = 0.5,
               colour = "red", linewidth = 0.5) +
  geom_point(data = siteZero, aes(PC3, PC4), colour = "grey60", size = 0.3) +
  geom_point(data = siteDet, aes(PC3, PC4, size = rate), shape = 21,
             colour = "black", fill = NA, stroke = 0.25) +
  scale_fill_viridis_c(name = "Share of PC1 x PC2\ninside 4-D region") +
  scale_size_area(max_size = 3, name = "Detection rate") +
  facet_wrap(~ species) + coord_equal() +
  labs(x = "PC3", y = "PC4",
       title = "Same 4-D regions projected on PC3 x PC4",
       subtitle = "Red outline = 2-D 90% region from PC3 and PC4 only. Grey dots = sites with no detections") +
  theme_bw()
ggsave("test_kde4D_fig2_projection_PC34.pdf", p2, width = 12, height = 8)

# figure 3: pairwise overlap, 2-D vs 4-D ---------------------------------------
p3 <- ggplot(overlaps, aes(overlap2D, overlap4D)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(size = 2) +
  geom_text(aes(label = paste(species1, "/", species2)), size = 2.5,
            hjust = 0, nudge_x = 0.01, check_overlap = TRUE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Jaccard overlap, PC1 x PC2 (2-D)", y = "Jaccard overlap, PC1-PC4 (4-D)",
       title = "Pairwise overlap: 2-D vs 4-D") +
  theme_bw()
ggsave("test_kde4D_fig3_overlap_2D_vs_4D.pdf", p3, width = 7, height = 7)

# check: does the grid resolution change the 4-D volume? ------------------------
# (run for one species; if volumes differ by more than a few %, use a finer grid)
sp <- testSpecies[1]
cellFracs <- c(1, 0.85, 0.75)
gridCheck <- sapply(cellFracs, function(f)
  kde_region(pcs, rates[[sp]], H4, xmin, xmax, grid_points(H4, xmin, xmax, f), prob)$volume)
names(gridCheck) <- paste0("cell_", cellFracs, "xBandwidth")
round(gridCheck, 1)

