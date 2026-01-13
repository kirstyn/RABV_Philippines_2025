# 1) Load libraries
# install.packages(c("terra","ggplot2","viridis"))
library(terra)
library(ggplot2)
library(viridis)

# 2) Read the raster
pop_raster <- rast("raw_data/gis_data/doh_data/worldpop/phl_pd_2020_1km.tif")

# 3) Inspect raster
pop_raster
summary(values(pop_raster))
plot(pop_raster)   # Quick base R plot

# Optional: check raster CRS
crs(pop_raster)

# 4) Convert to data frame for ggplot
# Only do this for smaller rasters; large rasters may require downsampling
pop_df <- as.data.frame(pop_raster, xy = TRUE)
names(pop_df) <- c("x", "y", "pop_density")

# Remove NA values
pop_df <- pop_df[!is.na(pop_df$pop_density), ]

# 5) Plot using ggplot2
ggplot(pop_df, aes(x = x, y = y, fill = pop_density)) +
  geom_raster() +
  scale_fill_viridis(
    option = "magma", 
    trans = "log10",   # population density is often skewed
    name = "Population Density\n(people per km²)"
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Philippines Population Density (2020, 1km)",
    caption = "Source: PHL_PD_2020_1km.tif"
  )
