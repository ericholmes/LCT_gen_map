# Load required libraries
library(tidyverse)
library(sf)
library(nhdplusTools)
library(elevatr)
library(hillshader)
library(raster)
library(terra)
library(tidyterra)
library(ggspatial)
library(ggpattern)
library(ggnewscale)
library(scales)

saveplot = F

# Get US state map data
us_states <- map_data("state")
nevada <- us_states[us_states$region == "nevada",]
polygon <- st_polygon(list(as.matrix(us_states[us_states$region == "nevada",c("long", "lat")])))
nevada_sf <- st_sf(id = nevada$region[1], geometry = st_sfc(polygon), crs = 4326)

#Load Pluvial Lake Lahontan shapefile
max_pl <- read_sf("data/MF-2323_Projected_Shapefiles/max_pl.shp")
max_pl <- st_as_sf(max_pl)
max_pl_wgs <- st_transform(max_pl, crs = 4326)

#Load LCT sample points
lct_points_all <- read.csv("data/RAPTURE data_Table_Thetas_hyb_PVA_03062020_withlatlons1.csv")

lct_sf_all <- st_as_sf(lct_points_all[is.na(lct_points_all$X) == F,], coords = c("X", "Y"), crs = 4326)

#Load NHD HUC8s
huc8_all <- nhdplusTools::get_huc(AOI = nevada_sf, type = "huc08")
huc8_all <- st_transform(huc8_all, st_crs(nevada_sf))
huc8_all$huc4 <- substr(huc8_all$huc8, 1, 4)
huc8_all$huc6 <- substr(huc8_all$huc8, 1, 6)

huc8_dissolve <- huc8_all %>%
  group_by(huc6) %>%  # Group by the categorical variable
  summarise(geometry = st_union(geometry)) %>%  # Combine geometries
  st_as_sf()  # Ensure the result is an sf object

huc8_single <- st_union(huc8_all)

# Filter HUC-8 polygons by intersection with Nevada boundary
huc8_nv <- st_intersection(huc8_all, nevada_sf)

#Load rivers
nev_nhd <- nhdplusTools::get_nhdplus(AOI = huc8_single, realization = "flowline", streamorder = 4)
nev_nhd <- st_transform(nev_nhd, st_crs(nevada_sf))

#load lakes
nev_nhd_lakes <- nhdplusTools::get_waterbodies(AOI = huc8_single)
nev_nhd_lakes <- st_transform(nev_nhd_lakes, st_crs(nevada_sf))

# Download elevation data using the 'elevatr' package
elevation_data <- get_elev_raster(locations = huc8_all, z = 7, prj = st_crs(4326)$proj4string)
elevation_data <- mask(elevation_data, huc8_dissolve)

# Convert raster to a data frame for ggplot
# elevation_df <- as.data.frame(as(elevation_data, "SpatialPixelsDataFrame"))
# elev_sf <- st_as_sf(elevation_df, coords = c("x", "y"))
# colnames(elevation_df) <- c("elevation", "x", "y")
# 
# hs <- hillshader(elevation_data)

r <- rast(elevation_data)
## Create hillshade effect
slope <- terrain(r, "slope", unit = "radians")
aspect <- terrain(r, "aspect", unit = "radians")
hill <- shade(slope, aspect, 30, 270)

# normalize names
names(hill) <- "shades"

# Hillshading, but we need a palette
pal_greys <- hcl.colors(1000, "Grays")

# Use a vector of colors
index <- hill %>%
  mutate(index_col = rescale(shades, to = c(1, length(pal_greys)))) %>%
  mutate(index_col = round(index_col)) %>%
  pull(index_col)

##Save map
if(saveplot == T){png("output/LCT_map_drape_rivers_hatch2.png",
          height = 7, width = 5, units = "in", res = 1000, family = "serif")}

ggplot(data = nevada_sf) +
  geom_spatraster(data = hill, fill = vector_cols, maxcell = Inf, alpha = 1) +
  geom_spatraster(data = r, maxcell = Inf, show.legend = F) +
  scale_fill_hypso_tint_c( limits = c(-200, 4500),
                           palette = "wiki-2.0_hypso",
                           alpha = 0.5,
                           labels = label_comma(),
                           breaks = c(
                             seq(-200, 1000, 200),
                             seq(1100, 2500, 100),
                             2600)) +
  geom_sf_pattern(data = max_pl_wgs[is.na(max_pl_wgs$LAKENAME) == F, ], alpha = .3,
                  pattern = "crosshatch",
                  pattern_density = 0.01,       # Adjust density of the pattern
                  pattern_angle = 45,         # Angle of the crosshatch
                  pattern_spacing = 0.005,     # Spacing between lines
                  pattern_fill = "white",     # Color of the pattern
                  pattern_color = "white", 
                  color = "white",fill = NA, pattern_size = .1) +
  geom_sf(data = nev_nhd_lakes[nev_nhd_lakes$ftype %in% "LakePond",], color = NA, fill = "darkslategrey") +
  geom_sf(data = nev_nhd, color = "darkslategrey", linewidth = .2) +
  geom_sf(data = huc8_all, fill = NA, color = "grey40") +
  geom_sf(data = huc8_dissolve, fill = NA, color = "black") +
  geom_sf(fill = NA, color = "black", linetype = 2) +
  # geom_sf(data = lct_sf) +
  # geom_sf(data = lct_sf_all, aes(color = Basin), size=1) +
  geom_sf(data = lct_sf_all, color = "black", size=1) +
  theme_bw() + #theme(legend.position.inside = c(0.8, 0.2)) +
  coord_sf(xlim = c(-121.1, -112.5), ylim = c(43.2, 34.2), expand = FALSE) +
  labs(title = "LCT genetics map", x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.2, line_width = 1,
                   pad_x = unit(.35, "in")) + 
  annotation_north_arrow(location = "bl", which_north = "false", 
                         style = north_arrow_fancy_orienteering(),
                         height = unit(0.3,"in"), width = unit(0.3,"in"),
                         pad_x = unit(.02, "in"),pad_y = unit(.02, "in"))

if(saveplot == T){dev.off()}


# Play with hillshade -----------------------------------------------------

## Libraries

library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(scales)

# Get the data
library(geodata)



r <- rast(elevation_data)
## Create hillshade effect
slope <- terrain(r, "slope", unit = "radians")
aspect <- terrain(r, "aspect", unit = "radians")
hill <- shade(slope, aspect, 30, 270)

# normalize names
names(hill) <- "shades"

# Hillshading, but we need a palette
pal_greys <- hcl.colors(1000, "Grays")

# Use a vector of colors
index <- hill %>%
  mutate(index_col = rescale(shades, to = c(1, length(pal_greys)))) %>%
  mutate(index_col = round(index_col)) %>%
  pull(index_col)

# Get cols
vector_cols <- pal_greys[index]
hill_plot <- ggplot() +
  geom_spatraster(
    data = hill, fill = vector_cols, maxcell = Inf,
    alpha = 1
  )

hill_plot

# Try some options, but we need to be aware of the values of our raster

r_limits <- minmax(r) %>% as.vector()

# Rounded to lower and upper 500
r_limits <- c(floor(r_limits[1] / 500), ceiling(r_limits[2] / 500)) * 500

# And making min value to 0.
r_limits <- pmax(r_limits, 0)

# Regular gradient
grad <- hypso.colors(10, "dem_poster")

if(saveplot == T){png("output/LCT_map_drape_rivers_new.png",
                      height = 7, width = 5, units = "in", res = 1000, family = "serif")}
# final plot
base_plot <- ggplot() +
  geom_spatraster(
    data = hill, fill = vector_cols, maxcell = Inf,
    alpha = 1
  ) +
  # Avoid resampling with maxcell
  geom_spatraster(data = r, maxcell = Inf, show.legend = F) +
  scale_fill_hypso_tint_c(
    limits = c(-200, 4500),
    palette = "wiki-2.0_hypso",
    alpha = 0.5,
    labels = label_comma(),
    # For the legend I use custom breaks
    breaks = c(
      seq(-200, 1000, 200),
      seq(1100, 2500, 100),
      2600
    )) +
  geom_sf_pattern(data = max_pl_wgs[is.na(max_pl_wgs$LAKENAME) == F, ], alpha = .3,
                  pattern = "crosshatch",
                  pattern_density = 0.01,       # Adjust density of the pattern
                  pattern_angle = 45,         # Angle of the crosshatch
                  pattern_spacing = 0.005,     # Spacing between lines
                  pattern_fill = "grey90",     # Color of the pattern
                  pattern_color = "grey90", 
                  color = "grey90",fill = NA, pattern_size = .1) +
  geom_sf(data = nev_nhd_lakes[nev_nhd_lakes$ftype %in% "LakePond",], color = NA, fill = "darkslategrey") +
  geom_sf(data = nev_nhd, color = "darkslategrey", linewidth = .15) +
  geom_sf(data = huc8_all, fill = NA, color = "grey40") +
  geom_sf(data = huc8_dissolve, fill = NA, color = "black") +
  geom_sf(fill = NA, color = "black", linetype = 2) +
  # geom_sf(data = lct_sf) +
  # geom_sf(data = lct_sf_all, aes(color = Basin), size=1) +
  geom_sf(data = lct_sf_all, color = "black", size=1) +
  theme_bw() + #theme(legend.position.inside = c(0.8, 0.2)) +
  coord_sf(xlim = c(-121.1, -112.5), ylim = c(43.2, 34.2), expand = FALSE) +
  labs(title = "LCT genetics map", x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.2, line_width = 1,
                   pad_x = unit(.35, "in")) + 
  annotation_north_arrow(location = "bl", which_north = "false", 
                         style = north_arrow_fancy_orienteering(),
                         height = unit(0.3,"in"), width = unit(0.3,"in"),
                         pad_x = unit(.02, "in"),pad_y = unit(.02, "in"))

base_plot

dev.off()

base_plot <- hill_plot +
  # Avoid resampling with maxcell
  geom_spatraster(data = r, maxcell = Inf) +
  scale_fill_hypso_tint_c(
    limits = r_limits,
    palette = "dem_poster",
    # palette = "wiki-2.0_hypso",
    # palette = "usgs-gswa2",
    # palette = "gmt_globe_hypso",
    alpha = 0.5,
    labels = label_comma()#,
    # For the legend I use custom breaks
  #   breaks = c(
  #     -400, 250, 500, 750, 1000,
  #     1250, 1500, 1750, 2000, 
  #     2500, 3000, 4500
  #     # seq(0, 4500, 500)
  # )
)
base_plot
