library(rmapshaper)
library(sf)
library(dplyr)
library(ggplot2)

library(sf)

phl_prov <- st_read("raw_data/gis_data/phl_adm_psa_namria_20231106_shp/phl_admbnda_adm2_psa_namria_20231106.shp") %>%
  select(ADM2_EN, ADM1_EN, geometry)

# Reduce coordinate precision to help with memory
st_geometry(phl_prov) <- st_set_precision(st_geometry(phl_prov), 1e5)

prov_list <- split(phl_prov, phl_prov$ADM1_EN)

prov_simp_list <- lapply(prov_list, function(x) {
  ms_simplify(x, keep = 0.05, keep_shapes = TRUE)
})

phl_prov_simp <- do.call(rbind, prov_simp_list)

ggplot(phl_prov_simp) +
  geom_sf(aes(fill = ADM1_EN), color = "white", size = 0.2) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Philippines Provinces (Simplified)", fill = "Region")

saveRDS(phl_prov_simp, "processed_data/gis_data/phl_prov_simp.rds")

# Some polygons are invalid, which means things like self-intersections, duplicate vertices, or degenerate edges
# Find which polygons are invalid
invalid_idx <- which(!st_is_valid(phl_prov))
invalid_idx

# Then fix them 
phl_prov_valid <- st_make_valid(phl_prov)

# Drop remaining problems
phl_prov_valid <- phl_prov_valid[st_is_valid(phl_prov_valid), ]

phl_region <- phl_prov_valid %>%
  group_by(ADM1_EN) %>%
  summarise()

# simplify each region
region_list <- split(phl_region, phl_region$ADM1_EN)
region_simp_list <- lapply(region_list, function(x) {
  ms_simplify(x, keep = 0.05, keep_shapes = TRUE)
})
# combine them back together again
phl_region_simp <- do.call(rbind, region_simp_list)

# plot simplified regions (and check how long it takes!)
ggplot(phl_region_simp) +
  geom_sf(aes(fill = ADM1_EN), color = "white", size = 0.3) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Philippines Regions (Simplified)", fill = "Region")

saveRDS(phl_region_simp, "processed_data/gis_data/phl_region_simp.rds")
