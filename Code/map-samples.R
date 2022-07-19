library(ggmap)
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggsn)

# define box

sbbox <- make_bbox(lon = c(-98.7, -97.9), lat = c(19.13,19.5), f = .1)
#sbbox <- make_bbox(lon = c(-98.4, -98.01), lat = c(19.1,19.3), f = .1)
library(tidyverse)
# get map
bris<- get_map(location=sbbox, color = "bw" )
#myMap<-get_map(location=sbbox,  source="osm",  color="bw")
# create map
brisbanemap <- ggmap(bris)
brisbanemap
air<- read.csv("/home/yendi/Downloads/coord.csv") %>% mutate_at(c(1,2,3), as.factor)
# display map
maps<-brisbanemap +
  geom_point(data = air, mapping = aes(
    x = Longitude, y = Latitude),#,color=Sitio, shape=Transecto), 
             size=7, color="red" ) +
  geom_text(data = air, 
            mapping = aes(x = Longitude+0.001,
                         y = Latitude,
                        label =pol),
            size = 4, color = "gray20", 
           fontface = "bold", 
          check_overlap = T)+
  coord_equal() + # needed for ggsn
  guides(scale="none") + 
  ggsn::north(x.min = -98.75, x.max = -98.65, 
              y.min = 19.45, y.max = 19.5, scale = 1.5) + 
  ggsn::scalebar(x.min = -98.85, x.max = -98.60, 
                 y.min = 19.43, y.max = 19.41, 
                 dist = 5, transform = TRUE, 
                 model = "WGS84", height = 0.5, 
                 st.dist = 0.5, dist_unit = "km")
maps


ggsave('fig_map_ok_scale.png',
       width = 12, height = 6, dpi = 300, plot = maps)

#other backgrounds

#watercolor
sbbox <- make_bbox(lon = c(-98.7, -97.9), lat = c(19.13,19.6), f = .1)
#sbbox <- make_bbox(lon = c(-98.63, -98.1), lat = c(19.23,19.48), f = .1)


map <- get_stamenmap( bbox = sbbox,  maptype = "terrain-background")
#map <- get_stamenmap( bbox = sbbox,  maptype = "watercolor")


wat<-ggmap(map) + 
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
  )+
  coord_equal() + # needed for ggsn
  guides(scale="none") + 
  ggsn::north(x.min = -98.75, x.max = -98.65, 
              y.min = 19.60, y.max = 19.64, scale = 1.5) + 
  ggsn::scalebar(x.min = -98.85, x.max = -98.60, 
                 y.min = 19.63, y.max = 19.61, 
                 dist = 5, transform = TRUE, 
                 model = "WGS84", height = 0.5, 
                 st.dist = 0.5, dist_unit = "km")+
  geom_point(data = air, mapping = aes(
    x = Longitude, y = Latitude),#,color=Sitio, shape=Transecto), 
    size=7, color="red" ) +
  geom_text(data = air, 
            mapping = aes(x = Longitude+0.001,
                          y = Latitude,
                          label =pol),
            size = 4, color = "gray20", 
            fontface = "bold", 
            check_overlap = T) #+

wat
ggsave('fig_map_ok_tb.png',
       width = 12, height = 6, dpi = 300, plot = wat)
#torner transparent
sbbox <- make_bbox(lon = c(-98.7, -97.9), lat = c(19.13,19.6), f = .1)

map <- get_stamenmap( bbox = sbbox,  maptype = "toner-background", scale=1, color="bw")
mapatt<- attributes(map)
map_transparent<- matrix(adjustcolor(map, alpha.f = 0.2), nrow = nrow(map))
attributes(map_transparent)<- mapatt

transp<-ggmap(map) + 
  theme_void() + 
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
  )+
  coord_equal() + # needed for ggsn
  guides(scale="none") + 
  ggsn::north(x.min = -98.75, x.max = -98.65, 
              y.min = 19.60, y.max = 19.64, scale = 1.5) + 
  ggsn::scalebar(x.min = -98.85, x.max = -98.60, 
                 y.min = 19.63, y.max = 19.61, 
                 dist = 5, transform = TRUE, 
                 model = "WGS84", height = 0.5, 
                 st.dist = 0.5, dist_unit = "km")+
  geom_point(data = air, mapping = aes(
    x = Longitude, y = Latitude),#,color=Sitio, shape=Transecto), 
    size=7, color="red" ) +
  geom_text(data = air, 
            mapping = aes(x = Longitude+0.001,
                          y = Latitude,
                          label =pol),
            size = 4, color = "gray20", 
            fontface = "bold", 
            check_overlap = T) #+
#ggmap(map_transparent) + 
 # theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white'))
ggsave('fig_map_ok_transp2.png',
       width = 12, height = 6, dpi = 300, plot = transp)

#world (long perspective)
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## [1] "sf"  
## [1] "data.frame"
library(ggspatial)
wor<-ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.1) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-100, -94), ylim = c(16, 20)) +
  geom_point(data = air, mapping = aes(x = Longitude, y = Latitude), 
             size=6, color="red") +
  geom_text(data = air, 
            mapping = aes(x = Longitude+0.001,
                          y = Latitude,
                          label =pol),
            size = 4, color = "gray20", 
            fontface = "bold", 
            check_overlap = T)
ggsave('fig_map_ok_world.png',
       width = 12, height = 6, dpi = 300, plot = wor)
