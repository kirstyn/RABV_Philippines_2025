# Install packages if not already installed
# install.packages(c("sf", "ggplot2"))

library(sf)
library(ggplot2)
library(dplyr)
library(readr)

# Load the shapefile
# Replace the path with your shapefile path (you only need the .shp file, sf will read the rest)
shp <- st_read("~/GitHub/speedier-analysis/gis_data/PHL/PHL_adm1.shp")
shp <- st_read("raw_data/gis_data/phl_adm_psa_namria_20231106_shp/phl_admbnda_adm2_psa_namria_20231106.shp")
shp <- st_read("raw_data/gis_data/phl_adm_psa_namria_20231106_shp/phl_admbnda_adm1_psa_namria_20231106.shp")

# Basic plot
plot(shp)

# Or a nicer ggplot2 version
ggplot(shp) +
  geom_sf(fill = "lightblue", color = "black") +
  theme_minimal() +
  ggtitle("Shapefile Plot")

# Centroids for each polygon
adm_centroids <- shp %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  )

# View
head(adm_centroids)

# Save if needed
st_write(adm_centroids, "processed_data/adm_centroids.shp")

ggplot() +
  geom_sf(data = shp, aes(fill = ADM1_EN), colour = "grey50") +   # polygons
  geom_point(data = adm_centroids, aes(x = lon, y = lat), 
             colour = "red", size = 2) +                         # centroids
  theme_minimal() +
  labs(title = "Administrative Boundaries with Predefined Centroids",
       fill = "Region")


# 1. Load existing centroid file
prov_centroids <- read_csv("raw_data/gis_data/PHL_centroids.csv")

# 2. Extract ADM1 (region) centroids from shapefile output
region_centroids <- adm_centroids %>%
  st_drop_geometry() %>%       # drop polygons
  select(Longitude = lon,
         Latitude = lat,
         Loc_ID = ADM1_EN) %>% # region name as location ID
  mutate(Type = "Region")      # tag them as regions

# 3. Bind together
all_centroids <- bind_rows(prov_centroids, region_centroids)

# 4. Save combined file
write_csv(all_centroids, "raw_data/gis_data/PHL_all_centroids.csv")

# using adm2 as shapefile import, extract a list of provinces and corresponding regions

adm2_table <- shp %>%
  st_drop_geometry() %>%              # drop the geometry, keep only attributes
  select(ADM2_EN, ADM1_EN) %>%        # keep only ADM2 and ADM1
  distinct()                          # ensure uniqueness

head(adm2_table)
write_csv(adm2_table, "raw_data/gis_data/PHL_provinceTo_region_mapping.csv")

