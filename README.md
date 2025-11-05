# Heat & Smoke Exposure Analysis for California Prisons

This repository contains a three-part R analysis pipeline to assess **heat exposure trends (2000–2023)** and **smoke exposure trends (2006–2023)** at California prison facilities, culminating in a visualization of both hazards.

---

## Final Output
A **combined map** (`smoke_heat_trends_map.png`) showing:
- **Left panel (a)**: Smoke exposure (current level + 18-year trend, 2006–2023)
- **Right panel (b)**: Heat exposure (current level + 23-year trend, 2000–2023)

Each panel includes:
- A main map of California with prison locations
- An inset scatterplot showing current exposure vs. trend magnitude, with point size indicating trend significance (p-value)

---

## Pipeline Overview (Run Scripts in Order)

### **1. Script #1: `heat_prison_trends_CA.R`**
**Purpose**: Load and clean ERA5 temperature data, match to CA prisons, compute heat metrics, and generate a preliminary heat map.

**Key Steps**:
- Loads daily ERA5 `temp_max` data (2000–2023) from `data/raw/ERA5_iWBGT_WBGT_2000_2023/`
- Matches prison centroids (`prison_boundaries_centroids_fl.csv`) using `ID_CROH`
- Calculates:
  - % of days ≥28°C per year (May–Oct focus in later scripts)
  - Current heat level (2019–2023 average)
  - 23-year linear trend in heat days
- Saves processed data: `daily_era5_heat_data_ca_prisons.csv`
- Outputs initial **univariate heat map**

> **Output**: Cleaned daily temperature dataset + baseline heat map

---

### **2. Script #2: `smoke_prison_trends_CA.R`**
**Purpose**: Refine heat analysis, compute missing data (2006–2017), and produce a formatted results table.

**Key Steps**:
- Reloads and reprocesses ERA5 data with **parallel loading** (`furrr`)
- Focuses on **May–October** season
- Computes same heat metrics with stricter data thresholds
- Identifies **missing years (2006–2017)** per site
- Generates:
  - `heat_prison_trends_ca.csv` – final heat metrics
  - HTML table with current heat, trends, p-values, and missing years

> **Output**: Final heat trend dataset + publication-ready table

---

### **3. Script #3: `smoke_heat_trend_maps.R`**
**Purpose**: Combine smoke and heat trend data into a **bivariate figure**.

**Key Steps**:
- Loads:
  - `smoke_prison_trends_CA.csv` → smoke metrics (2006–2023)
  - `heat_prison_trends_CA.csv` → heat metrics (from Script #2)
- Standardizes column names (`Lat`, `Lon`, `trend_pvalue`, etc.)
- Creates **two-panel figure** using `cowplot` and `biscale`:
  - **Smoke panel**: univariate viridis scale
  - **Heat panel**: 3×3 bivariate color grid
  - Both include **inset trend scatterplots** with p-value-sized points
- Saves final figure: `smoke_heat_trends_map.png`

> **Output**: Final visualization
