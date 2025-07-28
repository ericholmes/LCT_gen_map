# Load required libraries
library(tidyverse)
library(sf)
library(nhdplusTools)
library(ggnewscale)
library(elevatr)
library(hillshader)
library(raster)
library(tidyterra)
library(ggspatial)

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

huc8_dissolve <- huc8_all %>%
  group_by(huc4) %>%  # Group by the categorical variable
  summarise(geometry = st_union(geometry)) %>%  # Combine geometries
  st_as_sf()  # Ensure the result is an sf object

# Filter HUC-8 polygons by intersection with Nevada boundary
huc8_nv <- st_intersection(huc8_all, nevada_sf)

# Download elevation data using the 'elevatr' package
elevation_data <- get_elev_raster(locations = nevada_sf, z = 7, prj = st_crs(4326)$proj4string)
elevation_data <- mask(elevation_data, huc8_dissolve)

# Convert raster to a data frame for ggplot
elevation_df <- as.data.frame(as(elevation_data, "SpatialPixelsDataFrame"))
elev_sf <- st_as_sf(elevation_df, coords = c("x", "y"))
colnames(elevation_df) <- c("elevation", "x", "y")

hs <- hillshader(elevation_data)

##Save map
png("output/LCT_map_drape.png",
          height = 6, width = 4.5, units = "in", res = 1000, family = "serif")

ggplot(data = nevada_sf) +
  geom_spatraster(data = rast(hs), show.legend = F, alpha = .5) +
  scale_fill_gradient(low = "black", high = "white", name = "Hillshade", na.value = NA) +
  new_scale_fill() +
  geom_tile(data = elevation_df, aes(x = x, y = y, fill = elevation), alpha = .5, show.legend = F) +
  scale_fill_gradientn(colors=terrain.colors(10), name="Elevation (m)", na.value = NA) +
  geom_sf(data = max_pl_wgs[is.na(max_pl_wgs$LAKENAME) == F, ], fill = "navy", color = NA, alpha = .4) +
  geom_sf(data = huc8_dissolve, fill = NA, color = "grey40") +
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

dev.off()

