
library(ggplot2)
library(ggspatial)
YK = map_data("world")
YK2 = map_data("world") %>% filter(group==985)


sites <- st_as_sf(data.frame(longitude = c(-80.15, -80.1), 
                             latitude = c(26.5, 6.8)), 
                  coords = c("longitude", "latitude"), crs = 4326, 
                  agr = "constant")


florida<-ggplot() + 
  theme_bw() + 
  geom_polygon(data = YK, aes(x=long, y = lat, group = group), color = 'black', 
               fill = "cornsilk") + 
  geom_polygon(data = YK2, aes(x=long, y = lat, group = group), 
               color = 'black', fill = 'grey') +
  coord_sf(xlim = c(-110, -85), ylim = c(10, 30))+
    xlab("Longitude")+ ylab("Latitude")+
  theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(fill = NA))+
  geom_point(data = air, mapping = aes(
    x = Longitude, y = Latitude),#,color=Sitio, shape=Transecto), 
    size=1, color="black" )+ annotation_scale(location = "br", plot_unit = "km")

florida

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
library(ggrepel)
#register_google(key = "AIzaSyBCSQUO78ha22hVG3nABplHQmxnY-psL68")

air<- read.csv("/home/yendi/Documents/corredor_scripts/coord.csv") %>% mutate_at(c(1,2,3,7), as.factor)
air$Sites<- factor(air$Names, levels = c("1. Tetlanohcan 1", "2. Tetlanohcan 2",
                                         "3. La Ascención Huitzcolotepec", "4. Tepeticpac",
                                         "5. El Carmen Las Carrozas","6. La Caridad",
                                         "7. Tandi-Chala", "8. San Francisco Mitepec",
                                         "9. Piedra Canteada de San Felipe, Hidalgo",
                                         "10. Paraje El Madroño", "11. Ejido San Gabriel" ,                   
                                         "12. San Rafael Ixtapalucan 2" ))

#watercolor
sbbox <- make_bbox(lon = c(-98.7, -97.9), lat = c(19.13,19.6), f = .1)
#sbbox <- make_bbox(lon = c(-98.63, -98.1), lat = c(19.23,19.48), f = .1)
map <- get_stamenmap( bbox = sbbox,  maptype = "terrain-background")
#map <- get_stamenmap( bbox = sbbox,  maptype = "watercolor")
wat<-ggmap(map) + 
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2),
    panel.grid.major = element_line(colour = "red"), 
    panel.grid.minor = element_line(colour = "red", size = 0.25)
  )+
  coord_equal() + # needed for ggsn
  guides(scale="none") + 
  ggsn::north(x.min = -98.75, x.max = -98.65, 
              y.min = 19.60, y.max = 19.64, scale = 1.5) + 
  ggsn::scalebar(x.min = -98, x.max = -97.9, 
                 st.size = 2,
                 y.min = 19.63, y.max = 19.61, 
                 dist = 5, transform = TRUE, 
                 model = "WGS84", height = 0.5, 
                 st.dist = 0.5, dist_unit = "km") +
  geom_point(data = air%>% group_by(Sites) %>% summarise_if(is.numeric, mean), mapping = aes(
    x = Longitude, y = Latitude, fill=Sites), 
    size=4, pch=21 ) +
  geom_text_repel(data = air %>% group_by(Site) %>% summarise_if(is.numeric, mean), 
            mapping = aes(x = Longitude+0.002,
                          y = Latitude,
                          label =Site),
            size = 4, color = "black", box.padding = .3, segment.color = NA,
            fontface = "bold")+
  scale_fill_viridis_d(option ="turbo", name="Sites")+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))+
  geom_vline(xintercept = c(-98.75, -98.50, -98.25, -98 ),
             color = "grey", size = 0.3, linetype = "longdash")+
  geom_hline(yintercept = c(19.6,19.5,19.4,19.3,19.2,19.1 ),
             color = "grey", size = 0.3,  linetype = "longdash")+
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))


wat2<- wat+theme(legend.position = "none",
                 axis.text = element_text(colour = "black"))+ylab(
                   "Latitude")+xlab("Longitude")
wat2

library(cowplot)
leg<- get_legend(wat)
arrowA <- data.frame(x1 = 8, x2 = 20, y1 = 10, y2 = 14.5)
#arrowB <- data.frame(x1 = 18.5, x2 = 23, y1 = 8.5, y2 = 6.5)

maps<-ggplot() +
  coord_equal(xlim = c(0, 30), ylim = c(0, 20), expand = FALSE) +
  annotation_custom(ggplotGrob(florida), xmin = 0, xmax = 15, ymin = 0, 
                    ymax = 20) +
  annotation_custom(ggplotGrob(wat2), xmin = 15, xmax = 28, ymin =8, 
                    ymax =24) +
  annotation_custom(leg, xmin = 20, xmax = 28, ymin = 5, 
                    ymax = 7) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = arrowA, 
               arrow = arrow(), lineend = "round") +
  theme_void()
maps
ggsave('map_sites_new.png',
       width = 10, height = 6.5, dpi = 300, plot = maps)
