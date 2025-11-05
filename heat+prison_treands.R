###
# Heat Exposure Analysis - Part 1: Data Processing & Initial Mapping
###
# 
# Overview: This script processes ERA5 climate data (2000-2023) to analyze heat 
# exposure trends at California prison facilities. It loads daily temperature records,
# matches them with prison locations, calculates heat metrics (days ≥28°C), computes 
# 23-year trends, and generates an initial heat exposure map. This is the first of 
# three scripts in the analysis pipeline.

# Load required libraries using pacman
if (!require(pacman)) install.packages("pacman")
pacman::p_load(data.table, here, lubridate, ggplot2, maps, viridis, grid)

#### Set data paths ####
era5_path <- here("data/raw/ERA5_iWBGT_WBGT_2000_2023")
prison_csv <- here("data/raw/prison_boundaries_centroids_fl.csv")
output_csv <- here("data/processed/map_trends/daily_era5_heat_data_ca_prisons.csv")

#### Function to load and validate ERA5 data files ####
load_era5_data <- function(file_path) {
  filename <- basename(file_path)
  tryCatch({
    # Read only necessary columns and filter dates
    data <- fread(file_path, select = c("site", "temp_max", "year", "month", "day"))
    
    # Log column names and sample data for the first file
    if (!exists("first_file_logged", envir = .GlobalEnv)) {
      cat("Column names in", filename, ":\n")
      print(names(data))
      cat("First few rows:\n")
      print(head(data, 3))
      assign("first_file_logged", TRUE, envir = .GlobalEnv)
    }
    
    # Validate and create date column
    if (all(c("year", "month", "day") %in% names(data))) {
      data[, date := as.Date(sprintf("%d-%02d-%02d", year, month, day))]
    } else {
      warning("Missing year/month/day columns in", filename)
      data[, date := NA_Date_]
    }
    
    # Filter valid data and date range
    data <- data[!is.na(date) & !is.na(temp_max) & !is.na(site) &
                   between(date, as.Date("2000-01-01"), as.Date("2023-12-31"))]
    
    # Check for required columns
    if (!"temp_max" %in% names(data)) warning("No temp_max column in", filename)
    if (!"site" %in% names(data)) warning("No site column in", filename)
    
    return(data)
  }, error = function(e) {
    warning("Error reading", filename, ":", e$message)
    return(NULL)
  })
}

#### Load all ERA5 data files sequentially ####
cat("Loading ERA5 temperature data...\n")
heat_files <- list.files(era5_path, pattern = "\\.csv$", full.names = TRUE)
cat("Found", length(heat_files), "heat data files\n")

# Load files one by one using lapply
all_temp_data <- lapply(heat_files, load_era5_data)

# Filter out NULL results and combine
all_temp_data <- Filter(Negate(is.null), all_temp_data)

if (length(all_temp_data) == 0) {
  stop("No valid ERA5 CSV files found at ", era5_path)
}

cat("Successfully loaded", length(all_temp_data), "files\n")
temp_data <- rbindlist(all_temp_data, use.names = TRUE, fill = TRUE)
cat("Combined data has", nrow(temp_data), "rows\n")

# Log summary statistics
cat("Loaded temperature data for", length(unique(temp_data$site)), "unique sites\n")
cat("Date range:", format(min(temp_data$date, na.rm = TRUE)), "to", 
    format(max(temp_data$date, na.rm = TRUE)), "\n")
cat("Temperature range:", round(min(temp_data$temp_max, na.rm = TRUE), 1), "to",
    round(max(temp_data$temp_max, na.rm = TRUE), 1), "°C\n")
cat("Sample site IDs:", paste(head(unique(temp_data$site), 5), collapse = ", "), "\n")

if (nrow(temp_data) == 0) stop("No valid temperature data loaded")

# Clean up to free memory
rm(all_temp_data)
gc()

#### Load and filter prison facility data ####
prison_data <- fread(prison_csv, select = c("NAME", "ID_CROH", "Lat", "Lon", "STATE"))
ca_prisons <- prison_data[STATE == "CA"]
cat("\nLoaded", nrow(ca_prisons), "California prison facilities\n")

# Match temperature data with prison sites
temp_sites <- unique(temp_data$site)
prison_sites <- ca_prisons$ID_CROH
matched_sites <- intersect(temp_sites, prison_sites)
cat("Sites in temperature data:", length(temp_sites), "\n")
cat("Prison sites (CA):", length(prison_sites), "\n")
cat("Matched sites:", length(matched_sites), "\n")

# Prepare facility locations and join with temperature data
facility_locations <- ca_prisons[ID_CROH %in% matched_sites & !is.na(Lat) & !is.na(Lon),
                                 .(facility_name = NAME, site_id = ID_CROH, Lat, Lon)]
temp_data <- temp_data[facility_locations, on = .(site = site_id), nomatch = 0L][
  , facility_name := i.facility_name]

# Save processed temperature data
fwrite(temp_data, output_csv)

# Function to calculate heat exposure metrics
calculate_heat_metrics <- function(temp_data, facility_info) {
  # Add year and hot_day columns
  temp_data[, `:=`(year = year(date), hot_day = temp_max >= 28)]
  temp_data <- temp_data[between(year, 2000, 2023)]
  
  # Calculate annual heat metrics
  annual_heat <- temp_data[, .(total_days = .N,
                               hot_days = sum(hot_day, na.rm = TRUE),
                               heat_pct = (sum(hot_day, na.rm = TRUE) / .N) * 100),
                           by = .(facility_name, year)][
                             total_days >= 300]
  
  # Calculate current heat level (2019-2023)
  current_heat <- annual_heat[year >= 2019,
                              .(current_heat_level = mean(heat_pct, na.rm = TRUE)),
                              by = facility_name]
  
  # Calculate trends for facilities with sufficient data
  trend_data <- annual_heat[annual_heat[, .N, by = facility_name][N >= 20, facility_name],
                            .(trend_slope = lm(heat_pct ~ year)$coefficients[2],
                              trend_pvalue = summary(lm(heat_pct ~ year))$coefficients[2, 4],
                              years_analyzed = .N),
                            by = facility_name][
                              , total_trend_change := trend_slope * 23]
  
  # Combine results and merge with facility info
  results <- merge(current_heat, trend_data, by = "facility_name", all.x = FALSE)
  results <- merge(results, facility_info[, .(facility_name, Lat, Lon)],
                   by = "facility_name", all.x = FALSE)
  
  # Ensure valid ranges
  results[, `:=`(current_heat_level = pmax(0, pmin(100, current_heat_level)),
                 total_trend_change = pmax(-50, pmin(50, total_trend_change)))]
  
  return(results)
}

# Perform heat exposure analysis
cat("\nAnalyzing heat exposure...\n")
facility_robust_trend_data <- calculate_heat_metrics(temp_data, facility_locations)

cat("Analysis complete for", nrow(facility_robust_trend_data), "facilities\n")

if (nrow(facility_robust_trend_data) == 0) {
  stop("No facilities with sufficient data for visualization")
}

# Convert to data.frame for ggplot2 compatibility
facility_robust_trend_data <- as.data.frame(facility_robust_trend_data)

#### Create heat exposure map ####
cat("\nCreating heat exposure map...\n")
ca_map <- map_data("state") %>% filter(region == "california")

create_heat_map <- function(data, map_data) {
  ggplot() +
    geom_polygon(data = map_data, aes(x = long, y = lat, group = group),
                 fill = "white", color = "gray50", linewidth = 0.5) +
    geom_point(data = data, aes(x = Lon, y = Lat,
                                color = current_heat_level,
                                size = abs(total_trend_change)),
               alpha = 0.8) +
    scale_color_viridis_c(option = "plasma",
                          name = "Current (2019-2023)\n% Days ≥ 28°C",
                          labels = scales::percent_format(scale = 1, accuracy = 0.1),
                          trans = "sqrt", breaks = c(0, 5, 10, 20, 30, 40)) +
    scale_size_continuous(name = "23-Year Trend\nChange Magnitude\n(2000-2023)",
                          range = c(0.5, 5), breaks = c(0, 2, 5, 10, 15, 20),
                          labels = scales::percent_format(scale = 1),
                          limits = c(0, 25),
                          guide = guide_legend(override.aes = list(alpha = 0.8))) +
    coord_fixed(ratio = 1.3) +
    theme_minimal() +
    labs(title = "Heat Exposure at California Prison Facilities (2000-2023)",
         subtitle = "Current heat levels (2019-2023) and 23-year trend") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          plot.caption = element_text(size = 10, color = "gray60"),
          legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(1.0, "cm"),
          plot.margin = margin(t = 10, r = 10, b = 40, l = 10))
}

# Generate and display map
heat_map <- create_heat_map(facility_robust_trend_data, ca_map)
print(heat_map)

# Save outputs (uncomment to enable)
# ggsave(output_png, heat_map, width = 8, height = 10, dpi = 300)
# write.csv(facility_robust_trend_data, output_trends_csv, row.names = FALSE)

# Clean up
rm(temp_data, facility_locations, ca_prisons, prison_data)
gc()
