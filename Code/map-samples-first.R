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
  ggsn::scalebar(x.min = -98.85, x.max = -98.60, 
                 y.min = 19.63, y.max = 19.61, 
                 dist = 5, transform = TRUE, 
                 model = "WGS84", height = 0.5, 
                 st.dist = 0.5, dist_unit = "km") +
  geom_point(data = air, mapping = aes(
    x = Longitude, y = Latitude, color=Sites), 
    size=8 ) +
  geom_text(data = air, 
            mapping = aes(x = Longitude+0.002,
                          y = Latitude,
                          label =Site),
            size = 5, color = "#808B96", 
            fontface = "bold", 
            check_overlap = T)+
  scale_color_viridis_d(option ="turbo", name="Sites")+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))+
  geom_vline(xintercept = c(-98.75, -98.50, -98.25, -98 ),
             color = "grey", size = 0.3)+
  geom_hline(yintercept = c(19.6,19.5,19.4,19.3,19.2,19.1 ),
             color = "grey", size = 0.3)

wat
wat2<- wat+theme(legend.position = "none")
#wat
ggsave('map_sites.png',
       width = 15, height = 6.5, dpi = 300, plot = wat)
