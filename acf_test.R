library(readxl)
library(data.table)
library(dplyr)
library(zoo)
library(ggplot2)
library(segmented)
library(tseries) 

# ==============================================================================
# CONFIGURATION
# ==============================================================================

SITE_NAME  <- "Site1262"            
OUTPUT_DIR <- "Results/Diagnostics_Check"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ==============================================================================
# FUNCTIONS
# ==============================================================================

resample_timegrid <- function(age, d13C, dt_kyr = 1) {
  dt_Ma <- dt_kyr / 1000
  t_grid <- seq(min(age), max(age), by = dt_Ma)
  if(anyDuplicated(age)) {
    tmp <- aggregate(d13C, by=list(age), FUN=mean)
    age <- tmp[,1]; d13C <- tmp[,2]
  }
  d_interp <- approx(age, d13C, xout = t_grid, rule = 2)$y
  return(data.frame(age_Ma = t_grid, d13C = d_interp))
}

asymmetric_despike <- function(d13C, median_win_kyr=200, dt_kyr=1, mad_mult=3) {
  k <- round(median_win_kyr / dt_kyr)
  if(k %% 2 == 0) k <- k + 1
  
  # Extend median to edges for boundary coverage
  run_med <- zoo::rollmedian(d13C, k = k, fill = NA, align = "center")
  run_med <- zoo::na.fill(run_med, "extend")
  
  resid <- d13C - run_med
  mad_val <- mad(resid[!is.na(resid)], constant = 1)
  is_outlier <- (resid < (-mad_mult * mad_val)) & !is.na(resid)
  return(!is_outlier) 
}

fit_ramp <- function(age, d13C, psi_init = NULL) {
  ok <- !is.na(age) & !is.na(d13C)
  age <- age[ok]; d13C <- d13C[ok]
  lm0 <- lm(d13C ~ age)
  
  if(is.null(psi_init)) {
    r <- range(age)
    psi_init <- c(r[1] + 0.33*diff(r), r[1] + 0.66*diff(r))
  }
  
  seg <- tryCatch({
    segmented(lm0, seg.Z = ~age, npsi = 2, psi = psi_init, 
              control = seg.control(display = FALSE, n.boot = 0))
  }, error = function(e) return(NULL))
  
  if(is.null(seg)) return(NULL)
  
  fitted_y <- broken.line(seg)$fit
  return(list(fitted_df = data.frame(age=age, fitted=fitted_y)))
}

get_efolding <- function(residuals) {
  acf_vals <- acf(residuals, plot=FALSE, lag.max=100)$acf[-1]
  thresh <- 1/exp(1)
  min(which(acf_vals < thresh))
}

# ==============================================================================
# DATA LOADING
# ==============================================================================

d13c_data <- read_xlsx("Data/d13C_dd13C_benthic.xlsx") %>% as.data.table()
site_num <- gsub("Site", "", SITE_NAME) 
col_age  <- paste0("age_", site_num)
col_d13c <- paste0("d13c_", site_num)

df <- d13c_data[, .SD, .SDcols = c(col_age, col_d13c)]
colnames(df) <- c("age", "d13c")
df <- na.omit(df)

if(mean(df$age) > 1000) df$age <- df$age / 1000
df <- df[order(df$age),]

# Filter to common analysis interval (53.5 - 49.0 Ma)
df <- df[df$age <= 53.5 & df$age >= 49.0, ]

avg_res_kyr <- mean(diff(df$age)) * 1000
cat(sprintf("\nNative Resolution for %s: %.2f kyr\n", SITE_NAME, avg_res_kyr))

# ==============================================================================
# COMPARATIVE EXPERIMENTS
# ==============================================================================

# --- A: Native Resolution (Baseline) ---
res_native <- resample_timegrid(df$age, df$d13c, dt_kyr = round(avg_res_kyr))
fit_native <- fit_ramp(res_native$age_Ma, res_native$d13C, psi_init = c(51.5, 51.0))
resid_native <- res_native$d13C - fit_native$fitted_df$fitted
efold_native <- get_efolding(resid_native) * round(avg_res_kyr) 

# --- B: 1-kyr Resampling (Raw) ---
res_1kyr <- resample_timegrid(df$age, df$d13c, dt_kyr = 1)
fit_1kyr <- fit_ramp(res_1kyr$age_Ma, res_1kyr$d13C, psi_init = c(51.5, 51.0))
resid_1kyr <- res_1kyr$d13C - fit_1kyr$fitted_df$fitted
efold_1kyr <- get_efolding(resid_1kyr) * 1 

# --- C: 1-kyr Resampling + Despiking (Final Method) ---
mask <- asymmetric_despike(res_1kyr$d13C, median_win_kyr=200, dt_kyr=1)
clean_age <- res_1kyr$age_Ma[mask]
clean_val <- res_1kyr$d13C[mask]
fit_clean <- fit_ramp(clean_age, clean_val, psi_init = c(51.5, 51.0))
resid_clean <- clean_val - fit_clean$fitted_df$fitted
efold_clean <- get_efolding(resid_clean) * 1

# ==============================================================================
# RESULTS & INTERPRETATION
# ==============================================================================

cat("\n--- Autocorrelation Forensic Analysis ---\n")
cat(sprintf("1. Native (~%d kyr) e-folding: %d kyr\n", round(avg_res_kyr), efold_native))
cat(sprintf("2. 1-kyr Resampled   e-folding: %d kyr\n", efold_1kyr))
cat(sprintf("3. Despiked (Final)  e-folding: %d kyr\n", efold_clean))

cat("\n--- Interpretation ---\n")
diff_interp <- efold_1kyr - efold_native
diff_despike <- efold_clean - efold_1kyr

if(diff_interp > 5) {
  cat(sprintf(" -> Interpolation added %d kyr of artificial memory.\n", diff_interp))
} else {
  cat(" -> Interpolation had minimal effect on memory.\n")
}

if(diff_despike > 5) {
  cat(sprintf(" -> Despiking added %d kyr of smoothness (removing noise).\n", diff_despike))
} else {
  cat(" -> Despiking did not significantly change autocorrelation.\n")
}

# ==============================================================================
# VISUALIZATION (ACF COMPARISON)
# ==============================================================================

acf_A <- acf(resid_native, plot=FALSE, lag.max=100)
acf_B <- acf(resid_1kyr, plot=FALSE, lag.max=100)
acf_C <- acf(resid_clean, plot=FALSE, lag.max=100)

df_acf <- rbind(
  data.frame(Lag_kyr = acf_A$lag[-1] * round(avg_res_kyr), ACF = acf_A$acf[-1], Experiment = "A: Native Res"),
  data.frame(Lag_kyr = acf_B$lag[-1] * 1,                  ACF = acf_B$acf[-1], Experiment = "B: 1-kyr Resampled"),
  data.frame(Lag_kyr = acf_C$lag[-1] * 1,                  ACF = acf_C$acf[-1], Experiment = "C: Despiked (Final)")
)

p_comp <- ggplot(df_acf[df_acf$Lag_kyr <= 100, ], aes(x=Lag_kyr, y=ACF, color=Experiment)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept = 0, color="black") +
  geom_hline(yintercept = 1/exp(1), linetype="dashed", color="grey50", label="1/e") +
  scale_color_manual(values = c("black", "blue", "red")) +
  labs(title = paste("Source of Autocorrelation:", SITE_NAME),
       subtitle = "Comparing Original avg resolution vs. 1kyr-interpolation vs. Despiking",
       y = "Autocorrelation", x = "Lag (kyr)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_comp)
ggsave(file.path(OUTPUT_DIR, paste0(SITE_NAME, "_ACF_Comparison.pdf")), p_comp, width=8, height=6)