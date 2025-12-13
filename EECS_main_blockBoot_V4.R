library(readxl)
library(writexl)
library(data.table)
library(ggplot2)
library(dplyr)
library(zoo)
library(robustbase) 
library(segmented)
library(boot)
library(tseries) 

# ==============================================================================
# CONFIGURATION
# ==============================================================================

SITE_NAME  <- "Site1409"            
OUTPUT_DIR <- "Results/Main analysis"

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
  
  # Extend median to edges to catch boundary outliers
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
  
  seg <- segmented(lm0, seg.Z = ~age, npsi = 2, psi = psi_init, 
                   control = seg.control(display = FALSE, n.boot = 0))
  
  bps <- seg$psi[, "Est."]
  preds <- predict(seg, newdata = data.frame(age = bps))
  sorted_indices <- order(bps, decreasing = TRUE) 
  magnitude <- preds[sorted_indices[2]] - preds[sorted_indices[1]]
  fitted_y <- broken.line(seg)$fit
  
  return(list(model = seg, bps = bps, magnitude = magnitude, 
              fitted_df = data.frame(age=age, fitted=fitted_y)))
}

boot_ramp_residual <- function(age, d13C, psi_init, R=500, block_size=20) {
  fit_orig <- fit_ramp(age, d13C, psi_init)
  fitted_val <- fit_orig$fitted_df$fitted
  residuals  <- d13C - fitted_val
  n <- length(residuals)
  
  boot_results <- matrix(NA, nrow=R, ncol=3)
  colnames(boot_results) <- c("Onset", "End", "Magnitude")
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  for(i in 1:R) {
    res_star <- numeric(n)
    filled <- 0
    # Circular block bootstrap sampling
    while(filled < n) {
      start_idx <- sample(1:n, 1)
      idx_seq <- seq(start_idx, start_idx + block_size - 1)
      idx_seq <- ((idx_seq - 1) %% n) + 1 
      take   <- min(length(idx_seq), n - filled)
      res_star[(filled + 1):(filled + take)] <- residuals[idx_seq[1:take]]
      filled <- filled + take
    }
    d13C_star <- fitted_val + res_star
    
    try({
      fit_boot <- fit_ramp(age, d13C_star, psi_init = fit_orig$bps)
      if(!is.null(fit_boot)) {
        bps_sorted <- sort(fit_boot$bps, decreasing = TRUE)
        mag <- fit_boot$magnitude
        boot_results[i, ] <- c(bps_sorted[1], bps_sorted[2], mag)
      }
    }, silent = TRUE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ci_onset <- quantile(boot_results[,1], c(0.025, 0.5, 0.975), na.rm=TRUE)
  ci_end   <- quantile(boot_results[,2], c(0.025, 0.5, 0.975), na.rm=TRUE)
  ci_mag   <- quantile(boot_results[,3], c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  return(list(onset_ci = ci_onset, end_ci = ci_end, mag_ci = ci_mag))
}

# ==============================================================================
# DATA LOADING & PRE-PROCESSING
# ==============================================================================

cat("Loading Data...\n")
d13c_data <- read_xlsx("Data/d13C_dd13C_benthic.xlsx") %>% as.data.table()

cat("Processing", SITE_NAME, "...\n")
site_num <- gsub("Site", "", SITE_NAME) 
col_age  <- paste0("age_", site_num)
col_d13c <- paste0("d13c_", site_num)

df <- d13c_data[, .SD, .SDcols = c(col_age, col_d13c)]
colnames(df) <- c("age", "d13c")
df <- na.omit(df)

# Ensure Ma units
if(mean(df$age) > 1000) df$age <- df$age / 1000
df <- df[order(df$age),]

# Filter to common analysis interval (53.5 - 49.1 Ma)
df <- df[df$age <= 53.6 & df$age >= 49.1, ]

cat("1. Resampling to 1 kyr grid...\n")
res <- resample_timegrid(df$age, df$d13c, dt_kyr = 1)

cat("2. Removing asymmetric outliers (negative excursions)...\n")
mask_keep <- asymmetric_despike(res$d13C, median_win_kyr=200, dt_kyr=1, mad_mult=3)
clean_age <- res$age_Ma[mask_keep]
clean_val <- res$d13C[mask_keep]

cat("3. Fitting Segmented Ramp Model...\n")
fit_obj <- fit_ramp(clean_age, clean_val, psi_init = c(51.5, 51.0))

# Extract estimates
bps_sorted <- sort(fit_obj$bps, decreasing = TRUE)
onset_est <- bps_sorted[1]
end_est   <- bps_sorted[2]
mag_est   <- fit_obj$magnitude

cat(sprintf("   -> Est Onset:     %.3f Ma\n", onset_est))
cat(sprintf("   -> Est End:       %.3f Ma\n", end_est))
cat(sprintf("   -> Est Magnitude: %.3f permil\n", mag_est))

# Save main plot
p <- ggplot() +
  geom_line(data=res, aes(x=age_Ma, y=d13C), color="grey80") +
  geom_point(aes(x=clean_age, y=clean_val), color="steelblue", size=0.5, alpha=0.5) +
  geom_line(data=fit_obj$fitted_df, aes(x=age, y=fitted), color="red", linewidth=1.2) +
  geom_vline(xintercept = c(onset_est, end_est), linetype="dashed", color="red") +
  labs(title = paste("EECS Analysis:", SITE_NAME),
       subtitle = paste0("Onset: ", round(onset_est,3), " Ma | End: ", round(end_est,3), 
                         " Ma | Mag: +", round(mag_est,2), " permil"),
       x = "Age (Ma)", y = "d13C (permil)") +
  theme_minimal()

print(p)
plot_path <- file.path(OUTPUT_DIR, paste0(SITE_NAME, "_EECS_Plot.pdf"))
ggsave(plot_path, p, width=8, height=5, device = "pdf", useDingbats = FALSE)
cat("   -> Plot Saved:", plot_path, "\n")

# ==============================================================================
# RESIDUAL DIAGNOSTICS & BLOCK SIZE OPTIMIZATION
# ==============================================================================

cat("\n--- Running Residual Diagnostics ---\n")
resid <- clean_val - fit_obj$fitted_df$fitted

# Stationarity Tests
adf <- adf.test(resid)
kpss <- kpss.test(resid, null="Level")
lb <- Box.test(resid, lag=20, type="Ljung-Box")

cat(sprintf("   ADF p-value:  %.4f\n", adf$p.value))
cat(sprintf("   KPSS p-value: %.4f\n", kpss$p.value))
cat(sprintf("   Ljung-Box p:  %.4e\n", lb$p.value))

# Autocorrelation Structure
acf_obj <- acf(resid, lag.max=100, plot=T)
acf_vals <- acf_obj$acf[-1] 

# Calculate e-folding time (decay to 1/e)
efold_thresh <- 1 / exp(1) 
efold_lag <- min(which(acf_vals < efold_thresh))

# Calculate first local minimum (elbow)
diffs <- diff(acf_vals)
possible_mins <- which(diffs > 0)

if(length(possible_mins) > 0) {
  min_lag <- min(possible_mins)
} else {
  # Fallback for monotonic decay
  min_lag <- round(efold_lag * 1.5)
}

cat(sprintf("\n--- Suggested Block Sizes ---\n"))
cat(sprintf("   e-folding time:    %d kyr\n", efold_lag))
cat(sprintf("   First Local Min:   %d kyr\n", min_lag))

# Save Diagnostic Plots
diag_path <- file.path(OUTPUT_DIR, paste0(SITE_NAME, "_Diagnostics.pdf"))
pdf(diag_path, width=8, height=8)
par(mfrow=c(3,1))
plot(clean_age, resid, type="l", main=paste("Residuals:", SITE_NAME), 
     ylab="Residuals"); abline(h=0, col="red")
acf(resid, lag.max=100, main="ACF (Red=e-fold, Blue=Elbow)")
abline(v=efold_lag, col="red", lwd=2, lty=2)
abline(v=min_lag, col="blue", lwd=2, lty=2)
pacf(resid, lag.max=100, main="PACF")
par(mfrow=c(1,1))
dev.off()
cat("   -> Diagnostics Saved:", diag_path, "\n")

# ==============================================================================
# SENSITIVITY TEST (BLOCK BOOTSTRAP)
# ==============================================================================

# Define block sizes based on diagnostics
BLOCK_SIZES_TO_TEST <- c(efold_lag, min_lag, min_lag * 2)
R_ITERATIONS <- 1000 

cat(paste0("\n--- Running Block Bootstrap Sensitivity (R=", R_ITERATIONS, ") ---\n"))
cat("Testing Block Sizes:", paste(BLOCK_SIZES_TO_TEST, collapse=", "), "kyr\n")

sensitivity_results <- data.frame()

for(b_size in BLOCK_SIZES_TO_TEST) {
  cat(sprintf("\nRunning Block Size: %d kyr... ", b_size))
  
  cis <- boot_ramp_residual(clean_age, clean_val, psi_init = fit_obj$bps, 
                            R = R_ITERATIONS, block_size = b_size)
  
  row <- data.frame(
    Site = SITE_NAME,
    Method = paste0("Block Bootstrap (", b_size, "kyr)"),
    Block_Size_kyr = b_size,
    Onset_Median = cis$onset_ci[2],
    Onset_Lower_2.5 = cis$onset_ci[1],
    Onset_Upper_97.5 = cis$onset_ci[3],
    End_Median = cis$end_ci[2],
    End_Lower_2.5 = cis$end_ci[1],
    End_Upper_97.5 = cis$end_ci[3],
    Mag_Median = cis$mag_ci[2],
    Mag_Lower_2.5 = cis$mag_ci[1],
    Mag_Upper_97.5 = cis$mag_ci[3]
  )
  sensitivity_results <- rbind(sensitivity_results, row)
  cat("Done.")
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

point_est_row <- data.frame(
  Site = SITE_NAME,
  Method = "Point Estimate (Fit Ramp)",
  Block_Size_kyr = 0, 
  Onset_Median = onset_est,
  Onset_Lower_2.5 = NA,
  Onset_Upper_97.5 = NA,
  End_Median = end_est,
  End_Lower_2.5 = NA,
  End_Upper_97.5 = NA,
  Mag_Median = mag_est,
  Mag_Lower_2.5 = NA,
  Mag_Upper_97.5 = NA
)

final_table <- rbind(point_est_row, sensitivity_results)

cat("\nFinal Results:\n")
print(final_table)
cat("\nSaving to Excel...\n")

excel_path <- file.path(OUTPUT_DIR, paste0(SITE_NAME, "_EECS_Results.xlsx"))
write_xlsx(final_table, excel_path)
cat("   -> Saved:", excel_path, "\n")

# Generate Sensitivity Plot
source("utils/plot_block_sensitivity.R") 
plot_sensitivity_results(sensitivity_results, SITE_NAME, OUTPUT_DIR)

cat("\nDONE.\n")