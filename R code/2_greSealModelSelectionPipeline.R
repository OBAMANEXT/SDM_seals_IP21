# =============================================================================
# Grey Seal SDM — Model Selection Pipeline
# =============================================================================
# PURPOSE: Unified GAM/BAM model selection for grey seal habitat use
#
# WORKFLOW:
#   1. Load & clean all used:available ratio datasets
#   2. Assess collinearity of predictors (VIF / correlation)
#   3. Select optimal used:available ratio by comparing smooth stability (BIC)
#   4. Select optimal smoothing degree (k) per predictor (BIC)
#   5. Quantify deviance explained by each predictor (leave-one-out)
#
# MODEL STRUCTURE:
#   - Response  : used (1 = tracking fix, 0 = available point)
#   - Family    : binomial(logit)
#   - Method    : fREML with AR1 autocorrelation correction (rho per ratio)
#   - RE        : s(Year, bs="re") to account for interannual variation
#   - Weights   : 1 for used; 1/ratio for available (e.g. 0.1 for 1:10)
# =============================================================================

library(mgcv)
library(gratia)
library(itsadug)
library(dplyr)
library(ggplot2)
library(corrplot)
library(usdm)

# ---- Paths ------------------------------------------------------------------
BASE    <- "~/Documents/GitHub/SDM_seals_IP21"
SEAL    <- file.path(BASE, "Seal data")
FIG_DIR <- file.path(SEAL, "Figures")
MOD_DIR <- file.path(SEAL, "model output/optimal_smooth")

# =============================================================================
# 1. LOAD, CLEAN & PREPARE ALL RATIO DATASETS
# =============================================================================
# Each dataset represents a different used:available sampling ratio.
# For each:
#   - Drop rows with NA in key predictors
#   - Recode country to numeric factor (DNK=1, SWE=2, EST=3, FIN=4)
#   - Correct available-point weights (1/ratio for available, 1 for used)
#   - Extract Year as factor for random effect
#   - Mark temporal event starts for AR1 correction (itsadug::start_event)

COUNTRY_RECODE <- c("DNK" = "1", "SWE" = "2", "EST" = "3", "FIN" = "4")
NA_VARS        <- c("sediment", "WaveExpIndex", "rugosity", "sss")

# Config table: one row per ratio dataset.
# To add more ratios, extend the c() vectors here only — nothing else changes.
ratio_config <- tibble(
  file     = file.path(SEAL, "model data", paste0("SDM_data_1_", c(1, 5, 10), ".rda")),
  obj_name = paste0("SDM_data_1_", c(1, 5, 10)),
  label    = paste0("1:", c(1, 5, 10)),
  avail_wt = 1 / c(1, 5, 10)
)

prepare_ratio_data <- function(cfg_row) {
  env <- new.env()
  load(cfg_row$file, envir = env)
  dat <- get(cfg_row$obj_name, envir = env)

  dat <- dat %>%
    filter(if_all(all_of(NA_VARS), ~ !is.na(.x))) %>%
    mutate(
      ratio            = cfg_row$label,
      sediment         = as.factor(as.character(sediment)),
      weighting        = if_else(used == 0, cfg_row$avail_wt, weighting),
      country = as.factor(COUNTRY_RECODE[as.character(country)]),
      Year             = as.factor(format(date, "%Y")),
      log_dist2haulout = log(dist2haulout),
      log_WaveExpIndex = log(WaveExpIndex)
    )

  dat <- start_event(dat, column = "date", event = c("ref", "ref_burst"),
                     label.event = "Event")
  dat
}

# Load and prepare all ratio datasets into a named list
ratio_list <- purrr::map(
  split(ratio_config, seq_len(nrow(ratio_config))),
  prepare_ratio_data
) %>% setNames(ratio_config$label)

# Convenience aliases
SDM_data_1_1  <- ratio_list[["1:1"]]
SDM_data_1_5  <- ratio_list[["1:5"]]
SDM_data_1_10 <- ratio_list[["1:10"]]

# =============================================================================
# 2. COLLINEARITY ASSESSMENT
# =============================================================================
# Subset to used locations only (collinearity is a property of the predictors,
# not the availability design). Check:
#   - Pairwise Pearson correlations (corrplot)
#   - Variance Inflation Factors (vif / vifcor)
#
# Predictors with VIF > 10 or |r| > 0.7 flagged for removal.
# Final retained set: rugosity, dist2haulout, sboxy, sss, curr_velocity,
#                     WaveExpIndex (sediment retained as categorical)

assess_collinearity <- function(dat, fig_dir = FIG_DIR) {
  used_dat <- subset(dat, used == 1)

  all_preds <- c("bathymetry", "rugosity", "dist2haulout", "sss", "sbs", "sst",
                 "sbt", "ssoxy", "sboxy", "MixLayerThickness",
                 "ssh", "curr_velocity", "WaveExpIndex")

  preds_full <- used_dat[, all_preds]
  message("VIF — full predictor set:")
  print(vif(preds_full))

  pdf(file.path(fig_dir, "collinearityPlot.pdf"))
  corrplot(cor(preds_full), type = "upper")
  dev.off()

  message("vifcor threshold = 0.5:")
  print(vifcor(preds_full, th = 0.5))

  preds_final <- used_dat[, c("rugosity", "dist2haulout", "sboxy",
                               "sss", "curr_velocity", "WaveExpIndex")]
  message("VIF — reduced predictor set:")
  print(vif(preds_final))

  pdf(file.path(fig_dir, "collinearityPlot_final.pdf"))
  corrplot(cor(preds_final), type = "upper")
  dev.off()
}

assess_collinearity(SDM_data_1_10)

# =============================================================================
# 3. STEP 2 — OPTIMAL USED:AVAILABLE RATIO
# =============================================================================
# Fit the same BAM formula to each ratio dataset, with and without AR1.
# AR1 rho is estimated from residual ACF at lag = ratio value, which
# approximates autocorrelation induced by the sampling design.
# Compare models via BIC and smooth shape stability (gratia::compare_smooths).

BASE_FORMULA <- used ~ sediment +
  s(rugosity,          k = 3, bs = "cs") +
  s(log_dist2haulout,  k = 3, bs = "cs") +
  s(log_WaveExpIndex,  k = 3, bs = "cs") +
  s(curr_velocity,     k = 3, bs = "cs") +
  s(sss,               k = 3, bs = "cs") +
  s(sboxy,             k = 3, bs = "cs") +
  s(Year,              bs = "re")

fit_ratio_model <- function(dat, formula = BASE_FORMULA, rho = 0, ar_start = NULL) {
  do.call(bam, list(
    formula   = formula,
    family    = binomial(link = "logit"),
    data      = dat,
    weights   = dat[["weighting"]],
    method    = "fREML",
    discrete  = TRUE,
    na.action = na.fail,
    rho       = rho,
    AR.start  = ar_start
  ))
}


# ACF lag per ratio: lag = ratio value + 2 (empirical approach from original code)
# e.g. 1:1 → lag 3, 1:5 → lag 7, 1:10 → lag 12
acf_lags <- c("1:1" = 3, "1:5" = 7, "1:10" = 12)

# Two-pass fitting: estimate rho from pass-1 residuals, refit with AR1
ratio_models <- purrr::map(seq_along(ratio_list), function(i) {
  lbl <- names(ratio_list)[i]
  dat <- ratio_list[[i]]          # index by position, not name
  
  message("Fitting model for ratio: ", lbl)
  message("Rows in dat: ", nrow(dat))   # debug check
  
  mod0 <- fit_ratio_model(dat)
  rho  <- acf(resid(mod0), plot = FALSE)$acf[ acf_lags[lbl] ]
  message(lbl, " — estimated rho: ", round(rho, 4))
  fit_ratio_model(dat, rho = rho, ar_start = dat$start.event)
}) %>% setNames(names(ratio_list))

# Store rho values for reuse in steps 4 and 5
rho_vals <- purrr::map_dbl(seq_along(ratio_list), function(i) {
  lbl  <- names(ratio_list)[i]
  dat  <- ratio_list[[i]]
  mod0 <- fit_ratio_model(dat)
  acf(resid(mod0), plot = FALSE)$acf[ acf_lags[lbl] ]
}) %>% setNames(names(ratio_list))

# Convenience: rho for the chosen optimal ratio (update label after inspecting plots)
r1  <- rho_vals[["1:1"]]
r5  <- rho_vals[["1:5"]]
r10 <- rho_vals[["1:10"]]

# BIC comparison across ratios
bic_ratio <- data.frame(
  ratio = names(ratio_models),
  BIC   = purrr::map_dbl(ratio_models, BIC)
)
print(bic_ratio)

# Compare smooth shapes — stabilisation across ratios indicates optimal choice
cs <- compare_smooths(ratio_models[["1:1"]], 
                      ratio_models[["1:5"]], 
                      ratio_models[["1:10"]],
                      unconditional = TRUE)
cs_plot <- cs %>%
  mutate(.model = factor(.model, levels = names(ratio_list),
                         labels = names(ratio_list))) %>%
  filter(.smooth != "s(Year)")

# Add model labels manually (3 models × 6 smooths = 18 rows)
cs_plot$.model <- rep(names(ratio_models), times = length(unique(cs_plot$.smooth)))

# Unnest the data column
library(tidyr)
cs_unnested <- cs_plot %>%
  unnest(data)

# Check columns now available
names(cs_unnested)
head(cs_unnested)

# Pivot the predictor columns into a single x column
cs_long <- cs_unnested %>%
  tidyr::pivot_longer(
    cols      = c(curr_velocity, log_WaveExpIndex, log_dist2haulout,
                  rugosity, sboxy, sss),
    names_to  = "predictor",
    values_to = "x"
  ) %>%
  # Keep only rows where the smooth matches the predictor
  filter(paste0("s(", predictor, ")") == .smooth) %>%
  # Ensure ratio order is correct
  mutate(.model = factor(.model, levels = names(ratio_models)))

# Plot
ratio_plot <- ggplot(cs_long, aes(x = x, y = .estimate,
                                   colour = .model, fill = .model)) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se,
                  ymax = .estimate + 1.96 * .se), alpha = 0.2, colour = NA) +
  geom_line() +
  facet_wrap(~ .smooth, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(colour = "Ratio", fill = "Ratio", x = "", y = "Partial effect")

ggsave(file.path(FIG_DIR, "SelectRatio.png"), ratio_plot, width = 6, height = 12)


# ** Inspect BIC table and ratio_plot, then set OPTIMAL_RATIO below **
# Update this label and the corresponding data/rho objects to match your choice
OPTIMAL_RATIO   <- "1:10"           # <-- update after inspection
optimal_dat     <- ratio_list[[ OPTIMAL_RATIO ]]
optimal_rho     <- rho_vals[[ OPTIMAL_RATIO ]]
optimal_arstart <- optimal_dat$start.event

# Fit the final full model
gam_full <- fit_drop_model(NULL, optimal_dat, FULL_FORMULA, 
                            optimal_rho, optimal_arstart)

# Save to RData file
save(gam_full, file = file.path(SEAL, "Rdata/optimal_UNIFIED_model_1_10.RData"))

# To load it back later:
#load(file.path(SEAL, "Rdata/optimal_UNIFIED_model_1_10.RData"))

# =============================================================================
# 4. STEP 3 — OPTIMAL SMOOTHING DEGREE (k) PER PREDICTOR
# =============================================================================
# For each predictor, fit single-term BAMs with k = 3:8 on the optimal ratio
# dataset. Record BIC and deviance explained. Optimal k = lowest BIC.

smooth_vars <- list(
  log_dist2haulout = "Dist. to haulout",
  log_WaveExpIndex = "Wave exposure index",
  rugosity         = "Sea floor rugosity",
  sboxy            = "Sea bottom oxygen",
  sss              = "Sea surface salinity",
  curr_velocity    = "Current velocity"
)

run_k_selection <- function(var, label, dat, k_seq = 3:8, rho, ar_start) {
  purrr::map_dfr(k_seq, function(k) {
    formula <- as.formula(
      paste0("used ~ s(", var, ", k=", k, ", bs='cs') + s(Year, bs='re')")
    )
    mod <- bam(formula,
               family    = binomial(link = "logit"),
               data      = dat,
               weights   = dat$weighting,
               method    = "fREML",
               discrete  = TRUE,
               na.action = na.fail,
               rho       = rho,
               AR.start  = ar_start)
    tibble(
      var      = label,
      k        = k,
      BIC      = BIC(mod),
      dev.expl = round(summary(mod)$r.sq * 100, 1)
    )
  })
}

t1 <- Sys.time()
k_results <- purrr::map2_dfr(
  names(smooth_vars), smooth_vars,
  ~ run_k_selection(.x, .y, optimal_dat,
                    rho = optimal_rho, ar_start = optimal_arstart)
)
message("k selection runtime: ", round(difftime(Sys.time(), t1, units = "mins"), 1), " min")

k_results <- k_results %>%
  group_by(var) %>%
  mutate(selected = BIC == min(BIC)) %>%
  ungroup()

write.table(k_results,
  file.path(MOD_DIR, "optimal_k_unified_model_all_ratios.txt"),
  row.names = FALSE, sep = "\t")

k_plot <- ggplot(k_results, aes(x = k, y = BIC)) +
  geom_line(linetype = "dashed") +
  geom_point(aes(size = dev.expl, colour = selected)) +
  scale_size_continuous(range = c(1, 5), breaks = c(1, 15, 30)) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "steelblue")) +
  facet_wrap(~var, scales = "free_y", ncol = 2) +
  labs(x = "Smoothing factor (k)", y = "BIC",
       size = "Deviance explained (%)", colour = "Selected") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(file.path(FIG_DIR, "SelectSmooth.png"), k_plot, width = 6.5, height = 8)

# =============================================================================
# 5. STEP 4 — DEVIANCE EXPLAINED PER PREDICTOR (leave-one-out)
# =============================================================================
# Fit the full model, then refit dropping one predictor at a time.
# Change in deviance explained vs. full model = contribution of that predictor.
#
# ** Update k values below after inspecting k_results / SelectSmooth.png **

FULL_FORMULA <- used ~ sediment +
  s(rugosity,         k = 6, bs = "cs") +
  s(log_dist2haulout, k = 3, bs = "cs") +
  s(log_WaveExpIndex, k = 3, bs = "cs") +
  s(curr_velocity,    k = 3, bs = "cs") +
  s(sss,              k = 3, bs = "cs") +
  s(sboxy,            k = 6, bs = "cs") +
  s(Year,             bs = "re")

drop_terms <- list(
  "full model"          = NULL,
  "minus sediment"      = "sediment",
  "minus rugosity"      = "s(rugosity, k = 6, bs = \"cs\")",
  "minus dist2haulout"  = "s(log_dist2haulout, k = 3, bs = \"cs\")",
  "minus WaveExpIndex"  = "s(log_WaveExpIndex, k = 3, bs = \"cs\")",
  "minus curr_velocity" = "s(curr_velocity, k = 3, bs = \"cs\")",
  "minus sss"           = "s(sss, k = 3, bs = \"cs\")",
  "minus sboxy"         = "s(sboxy, k = 6, bs = \"cs\")",
  "minus Year"          = "s(Year, bs = \"re\")"
)

fit_drop_model <- function(drop_term, dat, base_formula, rho, ar_start) {
  f <- if (is.null(drop_term)) {
    base_formula
  } else {
    update(base_formula, paste0(". ~ . - ", drop_term))
  }
  do.call(bam, list(
    formula   = f,
    family    = binomial(link = "logit"),
    data      = dat,
    weights   = dat[["weighting"]],
    method    = "fREML",
    discrete  = TRUE,
    na.action = na.fail,
    rho       = rho,
    AR.start  = ar_start
  ))
}
dev_results <- purrr::imap_dfr(drop_terms, function(drop_term, model_label) {
  mod <- fit_drop_model(drop_term, optimal_dat, FULL_FORMULA,
                        optimal_rho, optimal_arstart)
  tibble(
    Model         = model_label,
    dev.explained = round(summary(mod)$dev.expl * 100, 1),
    BIC           = BIC(mod)
  )
})

print(dev_results)

write.table(dev_results,
  file.path(MOD_DIR, "deviance_explained.txt"),
  row.names = FALSE, sep = "\t")