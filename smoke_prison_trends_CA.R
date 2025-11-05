
# Heat Exposure Analysis for California Prisons (2000–2023)
# Computes heat exposure metrics and missing years (2006–2017) for prison facilities


# Load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  dplyr,         # Data manipulation
  lubridate,     # Date handling
  ggplot2,       # Visualization
  maps,          # Map data
  viridis,       # Color scales
  grid,          # Plot layout
  data.table,    # Efficient data processing
  furrr,         # Parallel processing
  here,          # Robust file paths
  kableExtra     # Table formatting
)

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)

# Define file paths using here()
era5_path <- here("data/raw/era5_iwbgt_wbgt_2000_2023")
prison_csv <- here("data/raw/prison_boundaries_centroids_fl.csv")
output_csv <- here("data/processed/daily_era5_heat_data_ca_prisons.csv")
output_trends_csv <- here("data/processed/heat_prison_trends_ca.csv")


# Function to load and validate ERA5 CSV files
load_era5_data <- function(file_path) {
  filename <- basename(file_path)
  tryCatch({
    data <- read_csv(file_path, col_types = cols(.default = "?", year = "i", month = "i", day = "i"))
    
    # Create date column
    if (all(c("year", "month", "day") %in% names(data))) {
      data$date <- as.Date(sprintf("%d-%02d-%02d", data$year, data$month, data$day))
    } else {
      warning(sprintf("Missing year/month/day columns in %s", filename))
      data$date <- NA
    }
    
    # Validate required columns
    required_cols <- c("site", "date", "temp_max")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      warning(sprintf("Missing columns %s in %s", paste(missing_cols, collapse = ", "), filename))
      return(NULL)
    }
    
    # Log column names and sample data for the first file
    if (!exists("first_file_logged", .GlobalEnv)) {
      cat(sprintf("Column names in %s:\n", filename))
      print(names(data))
      cat("First few rows:\n")
      print(head(data, 3))
      assign("first_file_logged", TRUE, .GlobalEnv)
    }
    
    return(data)
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", filename, e$message))
    return(NULL)
  })
}

# Load ERA5 data with parallel processing
cat("Loading ERA5 temperature data from", era5_dir, "...\n")
heat_files <- list.files(era5_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(heat_files) == 0) {
  stop(sprintf("No CSV files found in %s", era5_dir))
}
all_temp_data <- future_map_dfr(heat_files, load_era5_data, .options = furrr_options(seed = TRUE)) %>%
  filter(!is.na(date), !is.na(temp_max), !is.na(site),
         between(date, as.Date("2000-01-01"), as.Date("2023-12-31")))

# Log data summary
cat(sprintf("Loaded temperature data for %d unique sites\n", length(unique(all_temp_data$site))))
cat(sprintf("Date range: %s to %s\n", format(min(all_temp_data$date, na.rm = TRUE)), 
            format(max(all_temp_data$date, na.rm = TRUE))))
cat(sprintf("Temperature range: %.1f to %.1f °C\n", min(all_temp_data$temp_max, na.rm = TRUE), 
            max(all_temp_data$temp_max, na.rm = TRUE)))
cat(sprintf("Sample site IDs: %s\n", paste(head(unique(all_temp_data$site), 5), collapse = ", ")))
if (nrow(all_temp_data) == 0) {
  stop("No valid temperature data loaded")
}

# Load and filter prison facility data
cat(sprintf("\nLoading prison facility data from %s...\n", prison_csv))
prison_data <- read_csv(prison_csv, col_types = cols(.default = "?"))
ca_prisons <- prison_data %>% filter(STATE == "CA")
if (nrow(ca_prisons) == 0) {
  stop("No California prison facilities found in prison data")
}
cat(sprintf("Loaded %d California prison facilities\n", nrow(ca_prisons)))

# Match temperature data with prison sites
temp_sites <- unique(all_temp_data$site)
prison_sites <- ca_prisons$ID_CROH
matched_sites <- intersect(temp_sites, prison_sites)
cat(sprintf("Sites in temperature data: %d\n", length(temp_sites)))
cat(sprintf("Prison sites (CA): %d\n", length(prison_sites)))
cat(sprintf("Matched sites: %d\n", length(matched_sites)))

# Prepare facility locations and join with temperature data
facility_locations <- ca_prisons %>%
  filter(ID_CROH %in% matched_sites, !is.na(Lat), !is.na(Lon)) %>%
  select(facility_name = NAME, site_id = ID_CROH, Lat, Lon)
temp_data <- all_temp_data %>%
  left_join(select(facility_locations, site_id, facility_name),
            by = c("site" = "site_id")) %>%
  filter(!is.na(facility_name))

# Save processed temperature data
cat(sprintf("\nSaving processed temperature data to %s...\n", output_csv))
write_csv(temp_data, output_csv)

# Create missing_years_table (2006–2017)
cat("\nCreating missing years table (2006–2017)...\n")
counties <- unique(temp_data$site)
years <- 2006:2017
complete_grid <- CJ(MM_FIPS_COUNTY_NAME = counties, year = years)
existing_data <- unique(temp_data[, .(site, year = year(date))])
missing_years_dt <- complete_grid[!existing_data, on = .(MM_FIPS_COUNTY_NAME = site, year)]
missing_years_table <- missing_years_dt[, .(missing_years = paste(year, collapse = ",")),
                                        by = MM_FIPS_COUNTY_NAME]
missing_years_table <- merge(CJ(MM_FIPS_COUNTY_NAME = counties), 
                             missing_years_table, 
                             by = "MM_FIPS_COUNTY_NAME", 
                             all.x = TRUE)
missing_years_table[is.na(missing_years), missing_years := ""]
cat(sprintf("Missing years computed for %d sites\n", nrow(missing_years_table)))

# Function to calculate heat exposure metrics
calculate_heat_metrics <- function(temp_data, facility_info) {
  results <- tibble()
  
  for (facility in facility_info$facility_name) {
    facility_temps <- temp_data %>%
      filter(facility_name == facility, month(date) >= 5, month(date) <= 10) %>%
      mutate(year = year(date), hot_day = temp_max >= 28) %>%
      filter(between(year, 2000, 2023))
    
    if (nrow(facility_temps) == 0) next
    
    annual_heat <- facility_temps %>%
      group_by(year) %>%
      summarise(
        total_days = n(),
        hot_days = sum(hot_day, na.rm = TRUE),
        heat_pct = (hot_days / total_days) * 100,
        .groups = "drop"
      ) %>%
      filter(total_days >= 300)
    
    if (nrow(annual_heat) < 15) next
    
    current_heat_level <- annual_heat %>%
      filter(year >= 2019) %>%
      summarise(avg_heat = mean(heat_pct, na.rm = TRUE)) %>%
      pull(avg_heat)
    
    trend_stats <- if (nrow(annual_heat) >= 20) {
      model <- lm(heat_pct ~ year, data = annual_heat)
      list(total_trend_change = coef(model)[2] * 23,
           trend_pvalue = summary(model)$coefficients[2, 4])
    } else {
      list(total_trend_change = NA, trend_pvalue = NA)
    }
    
    results <- bind_rows(results, tibble(
      facility_name = facility,
      site = facility_info$site_id[facility_info$facility_name == facility],
      Lat = facility_info$Lat[facility_info$facility_name == facility],
      Lon = facility_info$Lon[facility_info$facility_name == facility],
      current_heat_level = current_heat_level,
      total_trend_change = trend_stats$total_trend_change,
      trend_pvalue = trend_stats$trend_pvalue,
      years_analyzed = nrow(annual_heat),
      hot_days_total = sum(annual_heat$hot_days, na.rm = TRUE)
    ))
  }
  
  return(results)
}

# Perform heat exposure analysis
cat("\nAnalyzing heat exposure...\n")
facility_robust_trend_data <- calculate_heat_metrics(temp_data, facility_locations) %>%
  filter(!is.na(current_heat_level), !is.na(total_trend_change)) %>%
  mutate(
    current_heat_level = pmax(0, pmin(100, current_heat_level)),
    total_trend_change = pmax(-50, pmin(50, total_trend_change))
  )
cat(sprintf("Analysis complete for %d facilities\n", nrow(facility_robust_trend_data)))
if (nrow(facility_robust_trend_data) == 0) {
  stop("No facilities with sufficient data for analysis")
}

# Join with missing_years_table and create kable table
heat_summary <- facility_robust_trend_data %>%
  left_join(missing_years_table, by = c("site" = "MM_FIPS_COUNTY_NAME")) %>%
  select(
    site, facility_name, current_heat_level, missing_years,
    total_trend_change, trend_pvalue, years_analyzed, hot_days_total
  ) %>%
  mutate(across(.cols = where(is.numeric), .fns = ~ round(.x, 2)))

cat("\nGenerating formatted table...\n")
kable_table <- heat_summary %>%
  kable(
    caption = "Heat Exceedances for California Prisons (2000–2023, May–October)",
    format = "html",
    col.names = c("Site ID", "Facility Name", "Current Heat Level (% Days ≥ 28°C, 2019–2023)",
                  "Missing Years (2006–2017)", "23-Year Trend Change (%)",
                  "Trend P-Value", "Years Analyzed", "Total Hot Days"),
    align = c("l", "l", "r", "l", "r", "r", "r", "r")
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

# Create heat exposure map
cat("\nCreating heat exposure map...\n")
ca_map <- map_data("state") %>% filter(region = "california")
create_heat_map <- function(data, map_data) {
  ggplot() +
    geom_polygon(data = map_data, aes(x = long, y = lat, group = group),
                 fill = "white", color = "gray50", linewidth = 0.5) +
    geom_point(data = data, aes(x = Lon, y = Lat,
                                color = current_heat_level,
                                size = abs(total_trend_change)),
               alpha = 0.8) +
    scale_color_viridis_c(
      option = "plasma",
      name = "Current (2019-2023)\n% Days ≥ 28°C",
      labels = scales::percent_format(scale = 1, accuracy = 0.1),
      trans = "sqrt", breaks = c(0, 5, 10, 20, 30, 40)
    ) +
    scale_size_continuous(
      name = "23-Year Trend\nChange Magnitude\n(2000-2023)",
      range = c(0.5, 5), breaks = c(0, 2, 5, 10, 15, 20),
      labels = scales::percent_format(scale = 1),
      limits = c(0, 25),
      guide = guide_legend(override.aes = list(alpha = 0.8))
    ) +
    coord_fixed(ratio = 1.3) +
    theme_minimal() +
    labs(
      title = "Heat Exposure at California Prison Facilities (2000-2023)",
      subtitle = "Current heat levels (2019-2023) and 23-year trend"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, color = "gray60"),
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.key.height = unit(1.0, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
    )
}

# Generate and display outputs
heat_map <- create_heat_map(facility_robust_trend_data, ca_map)
print(heat_map)
cat(sprintf("\nSaving trend data to %s...\n", output_trends_csv))
write_csv(heat_summary, output_trends_csv)
print(kable_table)
return(heat_summary)
