# =============================================================================
# Grey Seal Distribution Modelling - Data Preparation Pipeline
# =============================================================================
# PURPOSE: Prepare used/available RSF (Resource Selection Function) dataset
# for Grey Seal habitat modelling in the Baltic Sea.
#
# WORKFLOW:
#   1. Build reference raster grid (EPSG:3035, 1km resolution)
#   2. Filter tracking data (remove breeding/moulting months: Feb-Jun)
#   3. Prepare bathymetry raster (project, resample, mask)
#   4. Generate availability points for MULTIPLE used:available ratios
#   5. Extract environmental covariates from raster stack (per ratio)
#   6. Add distances, stratification indices and weightings (per ratio)
#   7. Save one dataset per ratio
#
# KEY INPUTS:
#   - trackingData_sampled.rda    : GPS tracking data (used locations)
#   - stack_env_cov.grd           : Raster stack of environmental covariates
#   - GreySeal_capture_sites.xlsx : Capture/haulout site coordinates
#   - EUSeaMap sediment raster    : Used to define sea-only availability mask
# =============================================================================

library(raster)
library(sp)
library(sf)
library(dplyr)
library(lubridate)
library(readxl)
library(mefa)

# ---- File paths (edit once here, used throughout) ---------------------------
BASE       <- "~/Documents/GitHub/SDM_seals_IP21"
SEAL_DIR   <- file.path(BASE, "Seal data")
COV_STACK  <- file.path(SEAL_DIR, "rasterStack/stack_env_cov.grd")
HAULOUT_XL <- file.path(SEAL_DIR, "haulout data/GreySeal_capture_sites.xlsx")

# ---- Ratios to generate -----------------------------------------------------
# Add or remove ratios here — everything downstream runs automatically.
# For each ratio R, availability points = R per used fix.
RATIOS <- c(1, 5, 10)

# =============================================================================
# 1. REFERENCE RASTER
# =============================================================================

build_reference_raster <- function() {
  r <- raster(
    xmn = 4255000, xmx = 5650000,
    ymn = 3420000, ymx = 4800000,
    res = c(1000, 1000),
    crs = "+init=epsg:3035"
  )
  r[] <- seq_len(ncell(r))
  r
}

refRaster <- build_reference_raster()

# =============================================================================
# 2. LOAD & FILTER TRACKING DATA
# =============================================================================
# Removes Feb–Jun (months 2–6) to exclude breeding and moulting seasons.

load_tracking_data <- function(path = file.path(SEAL_DIR, "Rdata/trackingData_sampled.rda")) {
  load(path, envir = environment())
  final_dat %>%
    mutate(Month = month(date)) %>%
    filter(Month < 2 | Month > 6)
}

final_dat <- load_tracking_data()

# =============================================================================
# 3. BATHYMETRY — project, resample, mask
# =============================================================================
# Defines valid sea-only cells for availability sampling.

prepare_bathymetry <- function(refRaster) {
  bathy_raw <- raster(file.path(SEAL_DIR,
    "sediment data/sediment_large_reduced_classes_entire_Baltic.grd"))
  bathy_proj   <- projectRaster(bathy_raw, crs = "+init=epsg:3035")
  bathy_resamp <- resample(bathy_proj, refRaster, method = "bilinear")
  mask_r       <- raster(file.path(SEAL_DIR, "rasterStack/maskMESS.grd"))
  mask(crop(bathy_resamp, mask_r), mask_r)
}

bathy <- prepare_bathymetry(refRaster)

# Remove tracking points on land
final_dat <- {
  sp_dat <- final_dat
  coordinates(sp_dat) <- ~x + y
  final_dat$bathymetry <- extract(bathy, sp_dat)
  final_dat[!is.na(final_dat$bathymetry), ]
}

# =============================================================================
# 4. GENERATE AVAILABILITY POINTS — multiple ratios
# =============================================================================
# For each fix, availability points are sampled once at the MAXIMUM ratio
# requested. Smaller ratios are created by subsampling from this master set,
# avoiding redundant raster operations.
#
# For each fix:
#   - Buffer 700 km (99th percentile movement distance)
#   - Clip bathymetry raster to buffer (sea cells only)
#   - Sample MAX_RATIO points from valid sea cells
#   - Subsampling to smaller ratios done in combine_used_available()

MAX_RATIO <- max(RATIOS)

generate_availability_for_fix <- function(row1, bathy, n = MAX_RATIO) {
  fix_sf   <- st_as_sf(row1, coords = c("lon", "lat"), crs = 4326) %>%
               st_transform(crs(bathy))
  fix_buff <- st_buffer(fix_sf, dist = 700000)
  bathy_clip <- crop(bathy, fix_buff, mask = TRUE)

  avail <- tryCatch(
    as.data.frame(sampleRandom(bathy_clip, size = n, xy = TRUE, na.rm = TRUE)),
    error = function(e) NULL
  )
  if (is.null(avail) || nrow(avail) < n) return(NULL)

  # Add a within-fix availability index (1:n) — used for subsampling later
  fix_rep        <- as.data.frame(rep(as.data.frame(row1), times = n))
  avail_out      <- cbind(fix_rep, avail)
  avail_out$used <- 0
  avail_out$avail_idx <- seq_len(n)   # index 1:MAX_RATIO per fix
  avail_out
}

generate_availability <- function(final_dat, bathy, max_ratio = MAX_RATIO) {
  final_dat2 <- final_dat[, 1:11]

  availDat_master <- purrr::map_dfr(
    unique(final_dat2$ref),
    function(id) {
      indv <- subset(final_dat2, ref == id)
      rownames(indv) <- NULL
      indv$row <- seq_len(nrow(indv))
      message("Processing individual: ", id, " (", nrow(indv), " fixes)")

      purrr::map_dfr(
        seq_len(nrow(indv)),
        function(j) generate_availability_for_fix(indv[j, ], bathy, n = max_ratio)
      )
    }
  )
  availDat_master
}

# Cache the master availability set (MAX_RATIO points per fix).
# All ratio-specific datasets are derived from this one file.
AVAIL_MASTER <- file.path(SEAL_DIR,
                  paste0("Rdata/AvailDat_master_1_", MAX_RATIO, ".rda"))

if (file.exists(AVAIL_MASTER)) {
  message("Loading cached master availability data from: ", AVAIL_MASTER)
  load(AVAIL_MASTER)
} else {
  availDat_master <- generate_availability(final_dat, bathy)
  save(availDat_master, file = AVAIL_MASTER)
  message("Master availability points saved to: ", AVAIL_MASTER)
}

# =============================================================================
# 5–10. BUILD ONE COMPLETE DATASET PER RATIO
# =============================================================================
# For each ratio R:
#   5. Subsample availability to R points per fix (avail_idx <= R)
#   6. Reproject availability coordinates and combine with used locations
#   7. Extract all environmental covariates
#   8. Add distance to haulout + distance from capture site
#   9. Compute stratification indices
#  10. Assign case weights (1 for used; 1/R for available)
#  11. Save

# ---- Helper functions (ratio-independent) -----------------------------------

combine_used_available <- function(final_dat, avail_sub) {
  avail_sp <- SpatialPointsDataFrame(
    coords      = avail_sub[, c(13, 14)],
    data        = avail_sub,
    proj4string = CRS("+init=epsg:3035")
  ) %>%
    spTransform(CRS("+proj=longlat +datum=WGS84")) %>%
    as.data.frame()

  avail_clean <- avail_sp[, c("ref","ref_burst","date","lon","lat",
                               "sex","age_class","device","country",
                               "x...13","y...14")]
  names(avail_clean)[10:11] <- c("x", "y")
  avail_clean$used <- 0

  used_clean <- final_dat[, c("ref","ref_burst","date","lon","lat",
                               "sex","age_class","device","country","x","y")]
  used_clean$used <- 1

  rbind(used_clean, avail_clean)
}

extract_stack_covariates <- function(dat, stack_path = COV_STACK) {
  UniStack <- stack(stack_path)

  layer_map <- list(
    sediment          =  1,
    bathymetry        =  2,
    rugosity          =  3,
    sss               =  4,
    sbs               =  5,
    sst               =  6,
    sbt               =  7,
    ssoxy             =  8,
    sboxy             =  9,
    ssh               = 10,
    curr_velocity     = 11,
    IceThickness      = 15,
    IceFraction       = 16,
    MixLayerThickness = 17,
    WaveExpIndex      = 18
  )

  sp_dat <- dat
  coordinates(sp_dat) <- ~x + y
  for (varname in names(layer_map)) {
    dat[[varname]] <- extract(UniStack[[ layer_map[[varname]] ]], sp_dat)
  }
  dat
}

add_dist_to_haulout <- function(dat, refRaster) {
  haulouts_sp <- read_excel(HAULOUT_XL) %>%
    distinct(location, .keep_all = TRUE) %>%
    { SpatialPointsDataFrame(.[, c(2, 3)], .,
        proj4string = CRS("+proj=longlat +datum=WGS84")) } %>%
    spTransform(CRS("+init=epsg:3035"))

  dist_rast <- distanceFromPoints(refRaster, haulouts_sp)

  coordinates(dat) <- ~x + y
  dat$dist2haulout <- extract(dist_rast, dat)
  as.data.frame(dat)
}

add_dist_from_capture <- function(dat) {
  capture <- read_excel(HAULOUT_XL) %>% select(1:3)
  dat     <- merge(dat, capture, by = "ref", all.x = TRUE)
  seal    <- cbind(dat$lon,         dat$lat)
  capt    <- cbind(dat$lon.capture, dat$lat.capture)
  dat$dist_home <- pointDistance(seal, capt, lonlat = TRUE)
  dat
}

# ---- Distance raster is the same for all ratios — compute once --------------
haulouts_sp <- read_excel(HAULOUT_XL) %>%
  distinct(location, .keep_all = TRUE) %>%
  { SpatialPointsDataFrame(.[, c(2, 3)], .,
      proj4string = CRS("+proj=longlat +datum=WGS84")) } %>%
  spTransform(CRS("+init=epsg:3035"))

dist_rast_haulout <- distanceFromPoints(refRaster, haulouts_sp)

# ---- Main loop: one complete SDM dataset per ratio --------------------------
for (R in RATIOS) {

  message("\n========== Processing ratio 1:", R, " ==========")

  # -- Step 5: subsample master availability to R points per fix --------------
  avail_sub <- availDat_master %>% filter(avail_idx <= R)

  # -- Step 6: combine with used locations ------------------------------------
  rsf_dat <- combine_used_available(final_dat, avail_sub)

  # -- Step 7: extract environmental covariates --------------------------------
  rsf_dat <- extract_stack_covariates(rsf_dat)

  # -- Step 8a: distance to haulout (reuse pre-computed raster) ---------------
  coordinates(rsf_dat) <- ~x + y
  rsf_dat$dist2haulout <- extract(dist_rast_haulout, rsf_dat)
  rsf_dat <- as.data.frame(rsf_dat)

  # Sanity check
  message("Mean dist2haulout — used: ",
    round(mean(rsf_dat$dist2haulout[rsf_dat$used==1], na.rm=TRUE)),
    " | available: ",
    round(mean(rsf_dat$dist2haulout[rsf_dat$used==0], na.rm=TRUE)))

  # -- Step 8b: distance from capture site ------------------------------------
  rsf_dat <- add_dist_from_capture(rsf_dat)

  # -- Step 9: stratification indices -----------------------------------------
  rsf_dat <- rsf_dat %>%
    mutate(
      strat_temp = sst   - sbt,
      strat_salt = sss   - sbs,
      strat_oxy  = ssoxy - sboxy
    )

  # -- Step 10: case weights --------------------------------------------------
  rsf_dat$weighting <- ifelse(rsf_dat$used == 1, 1, 1/R)

  # -- Step 11: save ----------------------------------------------------------
  # Dynamically name the object and file for each ratio
  obj_name <- paste0("SDM_data_1_", R)
  assign(obj_name, rsf_dat)

  out_path <- file.path(SEAL_DIR, "model data", paste0(obj_name, ".rda"))
  save(list = obj_name, file = out_path)

  message(
    "Saved: ", obj_name, "\n",
    "  Total rows : ", nrow(rsf_dat), "\n",
    "  Used       : ", sum(rsf_dat$used == 1), "\n",
    "  Available  : ", sum(rsf_dat$used == 0), "\n",
    "  File       : ", out_path
  )
}

message("\nAll ratios complete: ", paste(paste0("1:", RATIOS), collapse=", "))