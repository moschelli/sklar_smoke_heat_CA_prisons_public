
# Bivariate Mapping of Smoke and Heat Exposure Trends for California Prisons
# Creates a combined map of smoke (2006–2023) and heat (2000–2023) trends

# Load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  tidyverse,     # Data manipulation and visualization
  stringr,       # String operations
  cowplot,       # Plot composition
  biscale,       # Bivariate color scales
  tigris,        # Geographic data
  scales,        # Label formatting
  sf,            # Spatial data handling
  grid,          # Plot layout
  furrr,         # Parallel processing
  here,          # Robust file paths
  kableExtra     # Table formatting
)
options(tigris_use_cache = TRUE)

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)

# Define file paths using here()
smoke_csv <- here("data/raw/smoke_prison_trends_CA.csv") XxXXXXXXXXXX
heat_csv <- here("data/raw/heat_prison_trends_CA.csv")  XxXXXXXXXXXX
out_path <- here("data/processed/smoke_heat_trends_map.png").   XxXXXXXXXXXX

# Validate input files
cat("Checking input files...\n")
if (!file.exists(smoke_csv) || !file.exists(heat_csv)) {
  stop(sprintf("Input file(s) missing: %s, %s", smoke_csv, heat_csv))
}

# Function to standardize latitude/longitude columns
std_latlon <- function(df) {
  nm <- names(df)
  df <- df %>%
    rename_with(~ "Lat", .cols = any_of(c("latitude", "lat"))) %>%
    rename_with(~ "Lon", .cols = any_of(c("longitude", "long", "lon")))
  return(df)
}

# Function to pick first matching column
pick_first_col <- function(df, candidates, default = NA_real_) {
  found <- intersect(candidates, names(df))
  if (length(found)) as.numeric(df[[found[1]]]) else rep(default, nrow(df))
}

# Load and process smoke data
cat(sprintf("Loading smoke data from %s...\n", smoke_csv))
smoke <- read_csv(smoke_csv, col_types = cols(.default = "?", current_smoke_level = "n", total_trend_change = "n")) %>%
  std_latlon() %>%
  mutate(
    trend_pvalue = pick_first_col(cur_data(), c("trend_pvalue", "p_value", "trend_pval", "pval")),
    current_smoke_level = as.numeric(current_smoke_level),
    total_trend_change = as.numeric(total_trend_change),
    trend_pval_cat = case_when(
      is.na(trend_pvalue) ~ "≥0.05",
      trend_pvalue < 0.01 ~ "<0.01",
      trend_pvalue < 0.05 ~ "<0.05",
      TRUE ~ "≥0.05"
    ),
    trend_pval_cat = factor(trend_pval_cat, levels = c("<0.01", "<0.05", "≥0.05"), ordered = TRUE)
  ) %>%
  filter(!is.na(Lat), !is.na(Lon), !is.na(current_smoke_level), !is.na(total_trend_change))
cat(sprintf("Loaded smoke data for %d facilities\n", nrow(smoke)))

# Load and process heat data
cat(sprintf("Loading heat data from %s...\n", heat_csv))
heat <- read_csv(heat_csv, col_types = cols(.default = "?", current_heat_level = "n", total_trend_change = "n")) %>%
  std_latlon() %>%
  mutate(
    trend_pvalue = pick_first_col(cur_data(), c("trend_pvalue", "p_value", "trend_pval", "pval")),
    current_heat_level = as.numeric(current_heat_level),
    total_trend_change = as.numeric(total_trend_change),
    trend_pval_cat = case_when(
      is.na(trend_pvalue) ~ "≥0.05",
      trend_pvalue < 0.01 ~ "<0.01",
      trend_pvalue < 0.05 ~ "<0.05",
      TRUE ~ "≥0.05"
    ),
    trend_pval_cat = factor(trend_pval_cat, levels = c("<0.01", "<0.05", "≥0.05"), ordered = TRUE),
    col_cat = paste0(cut(current_heat_level, 3), "_", cut(total_trend_change, 3))
  ) %>%
  filter(!is.na(Lat), !is.na(Lon), !is.na(current_heat_level), !is.na(total_trend_change))
cat(sprintf("Loaded heat data for %d facilities\n", nrow(heat)))

# Create formatted table
cat("\nGenerating formatted table...\n")
heat_summary <- heat %>%
  select(site, facility_name, current_heat_level, missing_years, total_trend_change, trend_pvalue, years_analyzed, hot_days_total) %>%
  mutate(across(.cols = where(is.numeric), .fns = ~ round(.x, 2)))
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
print(kable_table)

# Load California map
cat("\nLoading California map data...\n")
ca <- tigris::states(cb = TRUE, year = 2024) %>%
  dplyr::filter(STATEFP == "06") %>%
  sf::st_transform(4326)

# Create smoke panel
cat("\nCreating smoke exposure panel...\n")
smoke_panel <- ggdraw() +
  draw_plot(
    ggplot(smoke, aes(x = Lon, y = Lat, color = current_smoke_level)) +
      geom_sf(data = ca, fill = "grey98", inherit.aes = FALSE) +
      geom_point(size = 2.5) +
      scale_color_viridis_c(
        option = "inferno",
        labels = label_number(suffix = "%"),
        guide = "none"
      ) +
      theme_void(),
    0, 0, 0.85, 0.95
  ) +
  draw_plot(
    ggplot(smoke, aes(x = current_smoke_level, y = total_trend_change,
                      size = trend_pval_cat, color = current_smoke_level)) +
      geom_point() +
      scale_color_viridis_c(option = "inferno", labels = label_number(suffix = "%"), guide = "none") +
      scale_size_manual(
        name = "trend p-val",
        values = c("<0.01" = 2.6, "<0.05" = 1.9, "≥0.05" = 1.3),
        breaks = c("<0.01", "<0.05", "≥0.05")
      ) +
      scale_x_continuous(name = "pct smoke days\n(2019–2023)", labels = label_number(suffix = "%")) +
      scale_y_continuous(name = "total change (2006–2023)", labels = label_number(suffix = "%")) +
      theme_half_open(12) +
      theme(
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "right",
        legend.box.spacing = unit(5, "pt"),
        legend.position.inside = c(0.8, 0.25)
      ),
    0.47, 0.25, 0.5, 0.68
  )

# Gridlines for heat inset
x_breaks <- heat$col_cat %>%
  str_split_i("_", 1) %>%
  str_remove_all("\\(|\\]") %>%
  strsplit(",") %>% list_c() %>% as.numeric() %>%
  unique() %>% sort() %>% .[-c(1, length(.))]
y_breaks <- heat$col_cat %>%
  str_split_i("_", 2) %>%
  str_remove_all("\\(|\\]") %>%
  strsplit(",") %>% list_c() %>% as.numeric() %>%
  unique() %>% sort() %>% .[-c(1, length(.))]

# Create heat panel
cat("\nCreating heat exposure panel...\n")
heat_axis_label <- "pct days ≥ 28°C temp\n(2019–2023)"
heat_panel <- ggdraw() +
  draw_plot(
    ggplot(heat, aes(x = Lon, y = Lat, color = col_cat)) +
      geom_sf(data = ca, fill = "grey98", inherit.aes = FALSE) +
      geom_point(size = 2.5) +
      scale_color_manual(values = bi_pal("DkViolet2", 3, preview = FALSE) %>% unname(), guide = "none") +
      theme_void(),
    0, 0, 0.85, 0.95
  ) +
  draw_plot(
    ggplot(heat, aes(x = current_heat_level, y = total_trend_change,
                     size = trend_pval_cat, color = col_cat)) +
      geom_vline(xintercept = x_breaks, alpha = 0.4) +
      geom_hline(yintercept = y_breaks, alpha = 0.4) +
      geom_point() +
      scale_color_manual(values = bi_pal("DkViolet2", 3, preview = FALSE) %>% unname(), guide = "none") +
      scale_size_manual(
        name = "trend p-val",
        values = c("<0.01" = 2.6, "<0.05" = 1.9, "≥0.05" = 1.3),
        breaks = c("<0.01", "<0.05", "≥0.05")
      ) +
      scale_x_continuous(name = heat_axis_label, labels = label_number(suffix = "%")) +
      scale_y_continuous(name = "total change (2000–2023)", labels = label_number(suffix = "%")) +
      theme_half_open(12) +
      theme(
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "right",
        legend.box.spacing = unit(5, "pt"),
        legend.position.inside = c(0.8, 0.25)
      ),
    0.47, 0.25, 0.5, 0.68
  )

# Compose and save figure
cat("\nComposing and saving figure...\n")
fig <- plot_grid(
  smoke_panel, heat_panel, nrow = 1,
  labels = c("a) smoke exposure trend (2006–2023)",
             "b) heat exposure trend (2000–2023)"),
  hjust = -0.05
)
print(fig)
ggsave(out_path, fig, width = 8 * 1.1, height = 4 * 1.1, bg = "white", dpi = 300)
cat(sprintf("Saved figure to %s\n", out_path))
