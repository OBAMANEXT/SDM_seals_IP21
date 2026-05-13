# =============================================================================
# Grey Seal SDM — Prediction Pipeline
# =============================================================================
# PURPOSE: Generate relative density and habitat suitability surfaces for
#          Grey Seals across the Baltic Sea using a fitted BAM/GAM model.
#
# TWO OUTPUT PRODUCTS:
#   A) Relative population density — predictions weighted by mean haulout
#      counts and summed across all sites, then normalised to % per km²
#   B) Habitat suitability — unweighted probability surface (0–1) from
#      a single distance-to-nearest-haulout raster across all sites
#
# UNCERTAINTY: Both products include mean + 95% Bayesian credible intervals
#   via posterior simulation (mvrnorm over GAM coefficient covariance matrix)
#
# WORKFLOW:
#   1. Load model, build reference raster, prepare covariate stack
#   2. [Density]      Loop over haulouts, predict + weight, accumulate
#   3. [Density]      MESS-mask, normalise, save, plot
#   4. [Hab.suit.]    Single prediction over all-haulout distance raster
#   5. [Hab.suit.]    MESS-mask, save, plot
#   6. Combine mean + CI panels into final figures
# =============================================================================

library(raster)
library(mgcv)
library(MASS)
library(sf)
library(ggplot2)
library(terra)
library(tidyterra)
library(grid)
library(plyr)

# ---- Paths (edit once) ------------------------------------------------------
BASE       <- "~/Documents/GitHub/SDM_seals_IP21"
SEAL     <- file.path(BASE, "Seal data")
GIS      <- file.path(SEAL, "GIS")
FIG_DIR  <- file.path(SEAL, "Figures")
LAY_DIR  <- file.path(SEAL, "final layers")
STACK_PATH  <- file.path(SEAL, "rasterStack/stack_env_cov.grd")
MESS_PATH   <- file.path(SEAL, "rasterStack/maskMESS.grd")
MODEL_PATH  <- file.path(SEAL, "Rdata/optimal_UNIFIED_model_1_10.RData")
HAULOUT_CSV <- file.path(SEAL, "haulout data/meanCounts_for_prediction.csv")

# Country recode: used to match factor levels in fitted model
COUNTRY_RECODE <- c("Germany"="1","Denmark"="1","Sweden"="2",
                    "Estonia"="3","Poland"="3","Finland"="4")

N_SIM      <- 25    # posterior simulations for uncertainty; increase for final outputs

# =============================================================================
# 1. SHARED SETUP — model, reference raster, covariate stack, haulouts
# =============================================================================

# ---- 1a. Load fitted model --------------------------------------------------
load(MODEL_PATH)
gam_final <- gam_full
rm(gam_full)

# ---- 1b. Reference raster (EPSG:3035, 1 km) ---------------------------------
# Resampled from the Copernicus bathymetry NetCDF to define the prediction grid.
# All prediction layers must align to this grid.

build_reference_raster <- function() {
  ref <- raster(xmn=4255000, xmx=5650000, ymn=3420000, ymx=4800000,
                res=c(1000,1000), crs="+init=epsg:3035")
  ref[] <- seq_len(ncell(ref))

  bathy_nc <- stack(file.path(SEAL,
    "rasterStack/bathy_cmems_mod_bal_phy_my_static_1705672302686.nc"))
  bathy_proj <- projectRaster(bathy_nc, crs="+init=epsg:3035")
  resample(bathy_proj, ref, method="bilinear")
}

refRaster <- build_reference_raster()
names(refRaster) <- "reference"

# ---- 1c. Covariate stack ----------------------------------------------------
# Adds log-transformed WaveExpIndex as an extra layer (required by model formula)

prepare_covariate_stack <- function(stack_path) {
  s <- stack(stack_path)
  names(s)[1]  <- "sediment"
  names(s)[18] <- "WaveExpIndex"
  log_wei <- log(s[[18]])
  names(log_wei) <- "log_WaveExpIndex"
  stack(s, log_wei)
}

UniStack <- prepare_covariate_stack(STACK_PATH)

# ---- 1d. Haulout data -------------------------------------------------------
load_haulouts <- function(csv_path) {
  h <- read.csv(csv_path)
  h$country <- revalue(h$Country, COUNTRY_RECODE)
  h
}

haulouts <- load_haulouts(HAULOUT_CSV)

# ---- 1e. MESS mask ----------------------------------------------------------
# Masks out grid cells with environmental conditions dissimilar to training
# data (MESS < 0 = extrapolation zone → set to NA)

load_mess <- function(mess_path, ref_rast) {
  m <- raster(mess_path)
  m[m < 0] <- NA
  crop(m, ref_rast)
}

MESS <- load_mess(MESS_PATH, refRaster)
MESS
# =============================================================================
# 2. BAYESIAN PREDICTION HELPER
# =============================================================================
# Given a prepared newdata frame and the fitted model, draws N_SIM samples
# from the posterior and returns mean + 95% CI on the response (probability) scale.
#
# Steps:
#   - draw beta vectors from multivariate normal (mean=coef, Sigma=vcov)
#   - multiply linear predictor matrix Xp by each beta draw
#   - apply inverse logit to convert to probability scale
#   - summarise across simulations

inv_logit <- function(x) exp(x) / (1 + exp(x))

bayesian_predict <- function(model, newdata, n_sim = N_SIM) {
  # Track which rows have complete data for all predictors in the model
  pred_vars <- all.vars(formula(model))
  pred_vars <- pred_vars[pred_vars %in% names(newdata)]
  complete_rows <- complete.cases(newdata[, pred_vars])
  
  # Initialise output vectors with NA for the full grid
  pred_mean  <- rep(NA, nrow(newdata))
  pred_lower <- rep(NA, nrow(newdata))
  pred_upper <- rep(NA, nrow(newdata))
  
  # Only predict for complete rows
  newdata_complete <- newdata[complete_rows, ]
  
  if (nrow(newdata_complete) == 0) {
    return(list(mean = pred_mean, lower = pred_lower, upper = pred_upper))
  }
  
  Vp       <- vcov(model, unconditional = TRUE)
  beta_hat <- coef(model)
  beta_sim <- mvrnorm(n_sim, mu = beta_hat, Sigma = Vp)

  Xp          <- predict(model, newdata = newdata_complete, type = "lpmatrix")
  linpred_mat <- Xp %*% t(beta_sim)
  prob_mat    <- apply(linpred_mat, 2, inv_logit)

  # Insert predictions back into the full-grid vectors
  pred_mean[complete_rows]  <- rowMeans(prob_mat)
  pred_lower[complete_rows] <- apply(prob_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  pred_upper[complete_rows] <- apply(prob_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

  list(mean = pred_mean, lower = pred_lower, upper = pred_upper)
}

# ---- Helper: build distance + Year layers for a given coordinate matrix -----
make_dist_year_layers <- function(coord_matrix, ref_rast, 
                                  training_dat = optimal_dat) {
  
  # Use the most common year in training data for prediction
  PRED_YEAR <- as.numeric(as.character(
    names(sort(table(training_dat$Year), decreasing = TRUE))[1]
  ))
  message("Predicting for year: ", PRED_YEAR)

  dist_rast     <- distanceFromPoints(ref_rast, coord_matrix)
  log_dist_rast <- log(dist_rast)
  names(log_dist_rast) <- "log_dist2haulout"

  year_rast <- setValues(dist_rast, PRED_YEAR)
  names(year_rast) <- "Year"

  full_stack <- stack(UniStack, log_dist_rast, year_rast)
  newdata    <- as.data.frame(full_stack)
  
  # Match factor levels from training data
  newdata$Year     <- factor(newdata$Year,     levels = levels(training_dat$Year))
  newdata$sediment <- factor(newdata$sediment, levels = levels(training_dat$sediment))

  list(stack = full_stack, newdata = newdata)
}


# ---- Helper: convert prediction vectors to rasters --------------------------
preds_to_rasters <- function(preds, template_rast) {
  list(
    mean  = setValues(template_rast, preds$mean),
    lower = setValues(template_rast, preds$lower),
    upper = setValues(template_rast, preds$upper)
  )
}

# =============================================================================
# 3. PRODUCT A — RELATIVE POPULATION DENSITY
# =============================================================================
# For each haulout site:
#   - Build a distance raster from that single site
#   - Predict habitat use probability across the grid (Bayesian)
#   - Weight predictions by mean seal count at that haulout
#   - Accumulate weighted predictions across all sites
#
# Using purrr::map to replace the explicit for-loop + rbind accumulation.
# Accumulation is done with Reduce("+") rather than repeated raster addition.

message("--- Product A: Relative density ---")

density_rasters <- purrr::map(unique(haulouts$X)[1:25], function(i) {
  ho1 <- subset(haulouts, X == i)
  message("  Haulout ", i, " of ", max(haulouts$X))

  layers  <- make_dist_year_layers(matrix(c(ho1$x, ho1$y), ncol=2), refRaster)
  preds   <- bayesian_predict(gam_final, layers$newdata)
  rasts   <- preds_to_rasters(preds, layers$stack[[1]])

  # Weight by mean count — seals contribute proportionally to local abundance
  list(
    mean  = rasts$mean  * ho1$MeanCount,
    lower = rasts$lower * ho1$MeanCount,
    upper = rasts$upper * ho1$MeanCount
  )
})

# Sum across all haulout sites
cumulative_mean  <- Reduce("+", purrr::map(density_rasters, "mean"))
cumulative_lower <- Reduce("+", purrr::map(density_rasters, "lower"))
cumulative_upper <- Reduce("+", purrr::map(density_rasters, "upper"))

# ---- MESS mask & normalise --------------------------------------------------
# Upper and lower CIs are normalised by the MEAN total (not their own sum)
# to preserve the interpretation of uncertainty relative to central estimate

mask_and_normalise <- function(rast, mess, norm_total) {
  masked <- mask(rast, mess)
  norm   <- masked / norm_total * 100
  norm
}

MESS_resampled <- resample(MESS, cumulative_mean, method = "ngb")

masked_mean  <- mask(cumulative_mean,  MESS_resampled)
masked_lower <- mask(cumulative_lower, MESS_resampled)
masked_upper <- mask(cumulative_upper, MESS_resampled)

total_mean   <- sum(masked_mean[], na.rm = TRUE)

norm_mean  <- mask_and_normalise(cumulative_mean,  MESS_resampled, total_mean)
norm_upper <- mask_and_normalise(cumulative_upper, MESS_resampled, total_mean)
norm_lower <- mask_and_normalise(cumulative_lower, MESS_resampled, total_mean)

# ---- Save raw cumulative rasters --------------------------------------------
save_rasters <- function(rast_list, names_vec, paths) {
  purrr::walk2(seq_along(rast_list), names_vec, function(i, nm) {
    r <- rast_list[[i]]
    names(r) <- nm
    writeRaster(r, paths[i], format = "raster", overwrite = TRUE)
  })
}

save_rasters(
  list(cumulative_mean, cumulative_upper, cumulative_lower),
  c("meanRelDens", "95upperRelDens", "95lowerRelDens"),
  file.path(LAY_DIR, c(
    "GS_RelDensity_MEAN_BAYES_cummulative.grd",
    "GS_RelDensity_UPPER_95CI_BAYES_cummulative.grd",
    "GS_RelDensity_LOWER_95CI_BAYES_cummulative.grd"))
)

# ---- Save normalised density rasters ----------------------------------------
save_rasters(
  list(norm_mean, norm_upper, norm_lower),
  c("meanRelDens", "95upperRelDens", "95lowerRelDens"),
  file.path(LAY_DIR, c(
    "GS_RelDensity_MEAN_BAYES.grd",
    "GS_RelDensity_UPPER_95CI_BAYES.grd",
    "GS_RelDensity_LOWER_95CI_BAYES.grd"))
)

# =============================================================================
# 4. PRODUCT B — HABITAT SUITABILITY
# =============================================================================
# Single prediction using minimum distance to ANY haulout across the full grid.
# No count weighting — gives the raw model probability surface (0–1).
# Useful for identifying suitable habitat independent of current population size.

message("--- Product B: Habitat suitability ---")

all_coords  <- matrix(c(haulouts$x, haulouts$y), ncol = 2)
hab_layers  <- make_dist_year_layers(all_coords, refRaster)
hab_preds   <- bayesian_predict(gam_final, hab_layers$newdata)
hab_rasts   <- preds_to_rasters(hab_preds, hab_layers$stack[[1]])

# MESS mask
hab_mean_masked  <- mask(hab_rasts$mean,  MESS_resampled)
hab_upper_masked <- mask(hab_rasts$upper, MESS_resampled)
hab_lower_masked <- mask(hab_rasts$lower, MESS_resampled)

save_rasters(
  list(hab_mean_masked, hab_upper_masked, hab_lower_masked),
  c("meanHabSuit", "95upperHabSuit", "95lowerHabSuit"),
  file.path(LAY_DIR, c(
    "GS_HabSuitability_MEAN_BAYESIAN.grd",
    "GS_HabSuitability_UPPER_95CI_BAYESIAN.grd",
    "GS_HabSuitability_LOWER_95CI_BAYESIAN.grd"))
)

# =============================================================================
# 5. SHARED PLOTTING INFRASTRUCTURE
# =============================================================================
# Country labels and shapefiles are identical across all four maps.
# Load once; reuse in every plot call.

WGS84 <- "+proj=longlat +ellps=WGS84 +no_defs"

load_shapefiles <- function(gis_dir) {
  coast <- sf::st_read(file.path(gis_dir, "Europe_border_polygon.shp"),
                       layer = "Europe_border_polygon") %>%
    st_transform(st_crs(WGS84))
  borders <- sf::st_read(file.path(gis_dir,
                         "maritime_bordersEPSG-4326 - WGS 84.shp")) %>%
    st_transform(st_crs(WGS84))
  list(coast = coast, borders = borders)
}

shapes <- load_shapefiles(GIS)

# Country annotation positions — defined once, applied everywhere
COUNTRY_LABELS <- tibble::tribble(
  ~label,  ~x,    ~y,
  "NOR",   9.5,  59.5,
  "EST",  25.3,  58.8,
  "LAT",  24.0,  56.75,
  "LIT",  23.2,  55.65,
  "RUS",  21.5,  54.6,
  "FIN",  23.5,  61.0,
  "SWE",  13.7,  57.0,
  "DNK",   9.9,  56.35,
  "POL",  17.0,  54.0,
  "GER",  12.5,  54.0
)

# ---- Base map theme (shared across all panels) ------------------------------
theme_seal_map <- function(show_axes = TRUE) {
  base <- theme_classic() +
    theme(
      panel.border    = element_rect(colour="black", fill=NA, linewidth=1),
      axis.line       = element_line(linewidth=0.25),
      plot.margin     = unit(c(0,0,0,0), "pt"),
      plot.title      = element_text(face="bold", size=14,
                                     margin=margin(l=10, t=20, b=-40))
    )
  if (!show_axes) {
    base <- base + theme(
      legend.position = "none",
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      plot.title      = element_text(face="bold", size=14,
                                     margin=margin(l=10, t=20, b=-20))
    )
  }
  base
}

# ---- Core map builder -------------------------------------------------------
# Builds a ggplot panel for one raster layer.
# `fill_scale`  : a ggplot2 scale object (e.g. scale_fill_distiller(...))
# `is_main`     : TRUE = show axes + legend; FALSE = inset panel (no axes)
# `title`       : panel title string
# `fill_mult`   : optional multiplier applied to fill values (density scaling)

build_map <- function(rast_path, col_name, fill_scale, title,
                      is_main = FALSE, fill_mult = 1,
                      shapes, country_labels,
                      xlim = c(8.8, 30.2), ylim = c(53.8, 66.2)) {

  r    <- raster(rast_path) %>% projectRaster(crs = WGS84)
  dat  <- as.data.frame(r, xy = TRUE)
  names(dat)[3] <- "val"

  p <- ggplot() +
    geom_raster(data = dat, aes(x = x, y = y, fill = val * fill_mult)) +
    fill_scale +
    scale_x_continuous("Longitude", limits = xlim, expand = c(0,0)) +
    scale_y_continuous("Latitude",  limits = ylim, expand = c(0,0)) +
    geom_sf(data = shapes$coast,   fill = "lightgrey", colour = NA,
            linewidth = 0.25) +
    geom_sf(data = shapes$borders, fill = NA, linetype = "dashed",
            colour = "black", linewidth = 0.1) +
    coord_sf(expand = FALSE) +
    ggtitle(title) +
    labs(x = if (is_main) "Longitude" else "",
         y = if (is_main) "Latitude"  else "",
         fill = "") +
    theme_seal_map(show_axes = is_main)

  if (is_main) {
    p <- p +
      theme(legend.position      = c(0.15, 0.9),
            legend.direction     = "horizontal",
            legend.key.width     = unit(1.0, "cm"),
            legend.background    = element_rect(fill = NA),
            plot.margin          = unit(c(0,0,0,-5), "pt"))
    # Add country labels on main panel only
    p <- p + geom_text(data = country_labels,
                       aes(x = x, y = y, label = label),
                       colour = "black", size = 3)
  }
  p
}

# ---- Compose and save a 3-panel figure (main + 2 insets) --------------------
save_three_panel <- function(main_plot, upper_plot, lower_plot, out_path,
                              width = 12.55, height = 12) {
  png(out_path, height = height, width = width,
      units = "in", res = 600, type = "quartz")
  grid.newpage()
  print(main_plot,  vp = viewport(width=0.6,  height=2,     x=0.35, y=0.50))
  print(upper_plot, vp = viewport(width=0.5,  height=0.345, x=0.795, y=0.698))
  print(lower_plot, vp = viewport(width=0.5,  height=0.345, x=0.795, y=0.344))
  dev.off()
}

# =============================================================================
# 6. PLOT RELATIVE DENSITY
# =============================================================================
# Density values are multiplied by 9 for visual scaling (as in original).
# Shared colour limits ensure mean, upper and lower panels are comparable.

DENSITY_SCALE_MAIN <- scale_fill_distiller(
  palette = "Spectral", limits = c(0, 0.006383791),
  na.value = "white", breaks = c(0.001, 0.003, 0.005)
)
DENSITY_SCALE_CI <- scale_fill_distiller(
  palette = "Spectral", limits = c(0, 0.006383791), na.value = "white"
)

dens_main  <- build_map(
  file.path(LAY_DIR, "GS_RelDensity_MEAN_BAYES.grd"),
  col_name   = "meanRelDens",
  fill_scale = DENSITY_SCALE_MAIN,
  title      = expression(bold("Relative population \n density (% at-sea km"^-2*")")),
  is_main    = TRUE, fill_mult = 9,
  shapes = shapes, country_labels = COUNTRY_LABELS
)

dens_upper <- build_map(
  file.path(LAY_DIR, "GS_RelDensity_UPPER_95CI_BAYES.grd"),
  col_name   = "95upperRelDens",
  fill_scale = DENSITY_SCALE_CI,
  title = "Upper 95% CI", fill_mult = 9,
  shapes = shapes, country_labels = COUNTRY_LABELS
)

dens_lower <- build_map(
  file.path(LAY_DIR, "GS_RelDensity_LOWER_95CI_BAYES.grd"),
  col_name   = "95lowerRelDens",
  fill_scale = DENSITY_SCALE_CI,
  title = "Lower 95% CI", fill_mult = 9,
  shapes = shapes, country_labels = COUNTRY_LABELS
)

save_three_panel(dens_main, dens_upper, dens_lower,
  file.path(FIG_DIR, "GS_RelDensity_plot_UNIFIED_1_10_BAYESIAN.png"))

# =============================================================================
# 7. PLOT HABITAT SUITABILITY
# =============================================================================

HAB_SCALE_MAIN <- scale_fill_terrain_c(
  limits = c(0, 1), na.value = "white", breaks = c(0.25, 0.5, 0.75)
)
HAB_SCALE_CI <- scale_fill_terrain_c(
  limits = c(0, 1), na.value = "white", breaks = c(0.25, 0.5, 0.75)
)

hab_main  <- build_map(
  file.path(LAY_DIR, "GS_HabSuitability_MEAN_BAYESIAN.grd"),
  col_name   = "meanHabSuit",
  fill_scale = HAB_SCALE_MAIN,
  title      = "Habitat suitability \n (prob.)",
  is_main    = TRUE,
  shapes = shapes, country_labels = COUNTRY_LABELS
)

hab_upper <- build_map(
  file.path(LAY_DIR, "GS_HabSuitability_UPPER_95CI_BAYESIAN.grd"),
  col_name   = "95upperHabSuit",
  fill_scale = HAB_SCALE_CI,
  title = "Upper 95% CI",
  shapes = shapes, country_labels = COUNTRY_LABELS
)

hab_lower <- build_map(
  file.path(LAY_DIR, "GS_HabSuitability_LOWER_95CI_BAYESIAN.grd"),
  col_name   = "95lowerHabSuit",
  fill_scale = HAB_SCALE_CI,
  title = "Lower 95% CI",
  shapes = shapes, country_labels = COUNTRY_LABELS
)

save_three_panel(hab_main, hab_upper, hab_lower,
  file.path(FIG_DIR, "GS_HabSuit_plot_UNIFIED_1_10_BAYESIAN.png"),
  height = 12, width = 12.55)