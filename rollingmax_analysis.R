library(readxl)
library(data.table)
library(dplyr)
library(zoo)
library(ggplot2)
library(reshape2)

d13c_data = read_xlsx("Data/d13C_dd13C_benthic.xlsx") %>% as.data.table()

min_age <- 49.1   # e.g. 49 Ma
max_age <- 52   # e.g. 54 Ma

site_list <- list(
  "Site1258" = d13c_data[, .(age = age_1258, d13c = d13c_1258)] %>% na.omit(),
  "Site1262" = d13c_data[, .(age = age_1262, d13c = d13c_1262)] %>% na.omit(),
  "Site1409" = d13c_data[, .(age = age_1409, d13c = d13c_1409)] %>% na.omit(),
  "Site1209" = d13c_data[, .(age = age_1209, d13c = d13c_1209)] %>% na.omit()
)

for(name in names(site_list)) {
  df <- site_list[[name]]
  
  # Convert to Ma if necessary
  if(mean(df$age, na.rm = TRUE) > 1000) {
    df$age <- df$age / 1000
  }
  
  df <- df %>% dplyr::filter(age >= min_age, age <= max_age)
  
  site_list[[name]] <- df
}

# ---- Average resolution on filtered interval ----
for(name in names(site_list)) {
  ages <- sort(site_list[[name]]$age)
  avg_res <- mean(diff(ages))
  cat(sprintf("%s Average Resolution (filtered): %.3f kyr\n", 
              name, avg_res * 1000))
}


dt = 1
test_windows = c(20, 40, 60, 100, 400)
colors = c("blue", "black", "red", "green", "orange")
names(colors) <- paste0("Win: ", test_windows, "kyr")

resample_timegrid <- function(age, value, dt_kyr) {
  # This function expects age in Ma
  dt_Ma <- dt_kyr / 1000 
  tgrid <- seq(min(age), max(age), by = dt_Ma)
  
  # Handle duplicates by averaging
  if(anyDuplicated(age)) {
    tmp <- aggregate(value, by=list(age), FUN=mean)
    age <- tmp[,1]; value <- tmp[,2]
  }
  
  vinterp <- approx(age, value, xout = tgrid, rule = 2)$y
  data.frame(age_Ma = tgrid, value = vinterp)
}

for(site_name in names(site_list)) { 
  
  df = site_list[[site_name]]
  
  # If ages are > 1000, assume they are kyr and convert to Ma
  if(mean(df$age, na.rm=TRUE) > 1000) {
    cat("  -> Detected kyr units. Converting to Ma.\n")
    df$age <- df$age / 1000
  }
  
  ages = df$age
  ord = order(ages)
  ages = ages[ord]
  vals = df$d13c[ord]
  
  res = resample_timegrid(ages, vals, dt)
  
  #rollin max for each window
  for(w in test_windows) {
    pts = round(w/dt)
    col_name = paste0("Win_", w, "kyr")
    res[[col_name]] = zoo::rollapply(res$value, width = pts, FUN = max, 
                                     align = "center", fill = NA)
  }
  
  #reshape for plotting
  env_data = melt(res, id.vars = "age_Ma", 
                  measure.vars = paste0("Win_", test_windows, "kyr"))
  levels(env_data$variable) <- gsub("_", ": ", levels(env_data$variable))
  
  # Generate Plot
  p <- ggplot() +
    # Layer 1: Raw Data (Background - same for all facets)
    geom_line(data = res, aes(x = age_Ma, y = value), 
              color = "grey80", size = 0.5) +
    # Layer 2: Envelopes
    geom_line(data = env_data, aes(x = age_Ma, y = value, color = variable), size = 1) +
    # Faceting: One column, rows defined by window variable
    facet_grid(variable ~ .) +
    # Styling
    scale_color_manual(values = colors) +
    labs(title = paste("Sensitivity of EECS estimate -", site_name),
         subtitle = "Comparison of different rolling windows",
         y = expression(delta^{13}*C~"(‰)"),
         x = "Age (Ma)") +
    theme_minimal() +
    theme(legend.position = "none",  # Legend redundant with facet labels
          strip.text = element_text(size = 12, face = "bold"))
  
  # Save Plot
  filename <- paste0("Sensitivity_", site_name, ".png")
  ggsave(filename, p, width = 8, height = 8) # Increased height for stacked plots
  print(p)
  cat("  -> Saved:", filename, "\n")
}
  
