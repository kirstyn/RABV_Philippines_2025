## Explore the full sequence dataset
library(dplyr)
library(ggplot2)
library(janitor)
library(lubridate)
library(scales)
library(rnaturalearth)
library(sf)
library(ggspatial)
library(rmapshaper)
library(viridis)
library(unikn)

full_data=read.csv("processed_data/processed_metadata/gathered_metadata/24Oct25_gathered_metadata_n797_raddl_and_manual_Corrected_stdGeo.csv")
phl <- st_read("raw_data/gis_data/phl_adm_psa_namria_20231106_shp/phl_admbnda_adm2_psa_namria_20231106.shp",
                    quiet = TRUE) %>%
  select(ADM1_EN, ADM2_EN, geometry)  # drop extra columns
# adm2 = province, adm1 = region
# -----------------------------
# 2. Host Distribution
# -----------------------------
# By scientific name
full_data %>% 
  tabyl(Host_scientific) %>% 
  arrange(desc(n))

# By common name
full_data %>% 
  tabyl(Host_common) %>% 
  arrange(desc(n))

# By host type
full_data %>% 
  tabyl(Host_type) %>% 
  adorn_pct_formatting(digits=1)

# Visualisation: Host type
ggplot(full_data, aes(x=Host_type, fill=Host_type)) +
  geom_bar() +
  theme_minimal() +
  labs(title="Number of samples by Host Type", x="Host Type", y="Count") +
  scale_fill_brewer(palette="Set2")

# -----------------------------
# 3. Temporal Coverage
# -----------------------------
# Convert Preferred_date to Date
full_data$Preferred_date <- dmy(full_data$Preferred_date)

# Samples per year
full_data %>%
  mutate(Year = year(Preferred_date)) %>%
  tabyl(Year) %>%
  ggplot(aes(x=factor(Year), y=n)) +
  geom_col(fill="steelblue") +
  theme_minimal() +
  labs(title="Number of samples per year", x="Year", y="Count")

# -----------------------------
# 4. Geographic Coverage
# -----------------------------

# Count samples per Region
region_counts <- full_data %>%
  tabyl(Region_std) %>%
  arrange(desc(n)) %>%
  filter(!is.na(Region_std))  # optional, remove NA
# Plot
ggplot(region_counts, aes(x = reorder(Region_std, n), y = n, fill = Region_std)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # horizontal bars for readability
  labs(
    x = "Region",
    y = "Number of samples",
    title = "Sample Distribution Across Regions"
  ) +
  theme_minimal(base_size = 14)

# Count samples per province and keep Region_std
province_counts <- full_data %>%
  filter(!is.na(Province_std)) %>%  # optional, remove NA
  group_by(Province_std, Region_std) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
# Plot
ggplot(province_counts, aes(x = reorder(Province_std, n), y = n, fill = Region_std)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # horizontal bars for readability
  labs(
    x = "Region",
    y = "Number of samples",
    title = "Sample Distribution Across Regions"
  ) +
  theme_minimal(base_size = 14)
# -----------------------------
# 4. Geographic Coverage- adm levels 
# -----------------------------
# Cases by province (ADM2)
cases_province <- full_data %>%
  group_by(Province_std) %>%
  summarise(n_cases = n(), .groups = "drop") 

# Cases by region (ADM1)
cases_region <- full_data %>%
  group_by(Region_std) %>%
  summarise(n_cases = n())

# Province-level map
phl_prov <- phl %>%
  left_join(cases_province, by = c("ADM2_EN" = "Province_std"))

phl_prov$geometry <- st_simplify(st_geometry(phl_prov), dTolerance = 2000)
# Region-level map
phl_region <- phl %>%
  group_by(ADM1_EN) %>%
  summarise() %>%
  left_join(cases_region, by = c("ADM1_EN" = "Region_std"))


# Province map
# ggplot(phl_prov_simp) +
#   geom_sf(aes(fill = n_cases), color = "white",size = 0.2) +  # white borders between provinces
#   scale_fill_viridis(option = "C", na.value = "grey90") +
#   theme_minimal() +
#   theme(
#     axis.text = element_blank(),      # remove axis labels
#     axis.ticks = element_blank(),     # remove axis ticks
#     panel.grid = element_blank(),     # remove gridlines
#     panel.background = element_rect(fill = "#f0e4d1", color = NA)  # light earth background
#   ) +
#   labs(title = "Rabies Cases by Province", fill = "Number of Cases")

# Region map
ggplot(phl_region) +
  geom_sf(aes(fill = n_cases), color = NA) +
  scale_fill_viridis(option = "C", na.value = "grey90") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),      # remove axis labels
    axis.ticks = element_blank(),     # remove axis ticks
    panel.grid = element_blank(),     # remove gridlines
    panel.background = element_rect(fill = "#f0e4d1", color = NA)  # light earth background
  ) +
  labs(title = "Rabies Cases by Region", fill = "Number of Cases")


# Visualisation: map scatter (requires Latitude/Longitude)
ggplot(full_data, aes(x=Longitude, y=Latitude, color=Region_std)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  labs(title="Sample locations by Region", x="Longitude", y="Latitude", color="Region") +
  scale_color_brewer(palette="Set3")

