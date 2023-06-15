#plot ccas
library(tidyverse)
library(ggvegan)
library(readxl)
map<- read.csv("Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))

fdat <- ggplot2::fortify(vares_cca1_veg)
map1<- map[match(rownames(spp_hell1), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()

map<- read.csv("Data/coord.csv") %>% mutate_at(
    c(1,2,3,7), as.factor) %>% mutate(SampleID= 
                                        paste0("P",pol, "S", Sitio,"T", Transecto ))


#a1
fdat <- ggplot2::fortify(vares_cca1_veg)
map1<- map[match(rownames(spp_hell1), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()
  
y<-ggordiplots::gg_ordiplot(vares_cca1_veg, groups = map1$Site, ellipse =  F, 
                            spiders = T, hull = F, plot = F )
z <- y$plot
a1<-z+theme_linedraw()+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  theme_linedraw()+
  scale_color_viridis_d(option ="turbo" )+#color of points +
    ggrepel::geom_text_repel(data = fdat %>% filter(
    Score=="biplot") %>%  
      mutate(CCA1=CCA1*2, CCA2=CCA2*2) ,
    aes(x=CCA1, y=CCA2, label= Label),
    color="red", segment.color = NA)+
  #geom_segment(data = fdat %>% filter(Score=="biplot") %>% 
   #         mutate(CCA1=CCA1*5, CCA2=CCA2*5) ,
    #     aes(x=0, y=0, xend=CCA1, yend=CCA2), 
     #   arrow=arrow(length=unit(0.15,"cm")),
      # size= 0.5,color="red")+
  theme_linedraw() +
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
#  ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
 #                             inner_join(taxonomy, 
  #                                       by = c("Label"="Feature.ID")) %>%
   #                           mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
    #                          filter(Label %in% varscho)  %>% 
     #                         mutate(tax=case_when(
      #                          tax=="pseudograminearum" ~ "Fusarium",
       #                         tax=="higginsianum" ~ "Colletotrichum",
        #                        tax=="oryzae" ~ "Aspergillus",
         #                       tax=="Colletitruchum fioriniae" ~ "Colletitruchum",
          #                      tax=="Aspergillus bombycis" ~ "Aspergillus",
           #                     tax=="Aspergillus glaucus" ~ "Aspergillus",
            #                    tax=="Aspergillus oryzae" ~ "Aspergillus",
             #                   tax=="Laccaria bicolor" ~ "Laccaria",
              #                  tax=="[Candida] auris" ~ "Clavispora",
               #                 tax=="Penicillium arizonense" ~ "Penicillium",
                #                TRUE ~ as.character(tax))) %>% 
                 #             mutate(CCA1=CCA1*scale, CCA2=CCA2*scale) ,
                  #          aes(x=CCA1, y=CCA2, label= tax), fontface="bold.italic", 
                   #         color="DarkBlue",segment.color = NA)+
  geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites"))+ylab("CCA2")+xlab("CCA1")
  
#a2
fdat <- ggplot2::fortify(vares_cca2_veg)
map1<- map[match(rownames(spp_hell2), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()

y<-ggordiplots::gg_ordiplot(vares_cca2_veg, groups = map1$Site, ellipse =  F, 
                            spiders = T, hull = F, plot = F )
z <- y$plot
a2<-z+theme_linedraw()+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  theme_linedraw()+
  scale_color_viridis_d(option ="turbo" )+#color of points +
  ggrepel::geom_text_repel(data = fdat %>% filter(
    Score=="biplot") %>%  
      mutate(CCA1=CCA1*1, CCA2=CCA2*1) ,
    aes(x=CCA1, y=CCA2, label= Label),
    color="red", segment.color = NA)+
  #geom_segment(data = fdat %>% filter(Score=="biplot") %>% 
  #         mutate(CCA1=CCA1*5, CCA2=CCA2*5) ,
  #     aes(x=0, y=0, xend=CCA1, yend=CCA2), 
  #   arrow=arrow(length=unit(0.15,"cm")),
  # size= 0.5,color="red")+
  theme_linedraw() +
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  #  ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
  #                             inner_join(taxonomy, 
  #                                       by = c("Label"="Feature.ID")) %>%
  #                           mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
  #                          filter(Label %in% varscho)  %>% 
  #                         mutate(tax=case_when(
  #                          tax=="pseudograminearum" ~ "Fusarium",
  #                         tax=="higginsianum" ~ "Colletotrichum",
  #                        tax=="oryzae" ~ "Aspergillus",
  #                       tax=="Colletitruchum fioriniae" ~ "Colletitruchum",
  #                      tax=="Aspergillus bombycis" ~ "Aspergillus",
#                     tax=="Aspergillus glaucus" ~ "Aspergillus",
#                    tax=="Aspergillus oryzae" ~ "Aspergillus",
#                   tax=="Laccaria bicolor" ~ "Laccaria",
#                  tax=="[Candida] auris" ~ "Clavispora",
#                 tax=="Penicillium arizonense" ~ "Penicillium",
#                TRUE ~ as.character(tax))) %>% 
#             mutate(CCA1=CCA1*scale, CCA2=CCA2*scale) ,
#          aes(x=CCA1, y=CCA2, label= tax), fontface="bold.italic", 
#         color="DarkBlue",segment.color = NA)+
geom_label(
  data = y$df_mean.ord,
  aes(x = x, y = y, label=Group), 
  label.padding = unit(0.15, "lines"),label.size = 0.4,
)+guides(
  color=guide_legend(title="Sites"))+ylab("CCA2")+xlab("CCA1")




#a3
fdat <- ggplot2::fortify(vares_cca3_veg)
map1<- map[match(rownames(spp_hell3), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()

y<-ggordiplots::gg_ordiplot(vares_cca3_veg, groups = map1$Site, ellipse =  F, 
                            spiders = T, hull = F, plot = F )
z <- y$plot
a3<-z+theme_linedraw()+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  theme_linedraw()+
  scale_color_viridis_d(option ="turbo" )+#color of points +
  ggrepel::geom_text_repel(data = fdat %>% filter(
    Score=="biplot") %>%  
      mutate(CCA1=CCA1*1, CCA2=CCA2*1) ,
    aes(x=CCA1, y=CCA2, label= Label),
    color="red", segment.color = NA)+
  #geom_segment(data = fdat %>% filter(Score=="biplot") %>% 
  #         mutate(CCA1=CCA1*5, CCA2=CCA2*5) ,
  #     aes(x=0, y=0, xend=CCA1, yend=CCA2), 
  #   arrow=arrow(length=unit(0.15,"cm")),
  # size= 0.5,color="red")+
  theme_linedraw() +
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  #  ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
  #                             inner_join(taxonomy, 
  #                                       by = c("Label"="Feature.ID")) %>%
  #                           mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
  #                          filter(Label %in% varscho)  %>% 
  #                         mutate(tax=case_when(
  #                          tax=="pseudograminearum" ~ "Fusarium",
  #                         tax=="higginsianum" ~ "Colletotrichum",
  #                        tax=="oryzae" ~ "Aspergillus",
  #                       tax=="Colletitruchum fioriniae" ~ "Colletitruchum",
  #                      tax=="Aspergillus bombycis" ~ "Aspergillus",
#                     tax=="Aspergillus glaucus" ~ "Aspergillus",
#                    tax=="Aspergillus oryzae" ~ "Aspergillus",
#                   tax=="Laccaria bicolor" ~ "Laccaria",
#                  tax=="[Candida] auris" ~ "Clavispora",
#                 tax=="Penicillium arizonense" ~ "Penicillium",
#                TRUE ~ as.character(tax))) %>% 
#             mutate(CCA1=CCA1*scale, CCA2=CCA2*scale) ,
#          aes(x=CCA1, y=CCA2, label= tax), fontface="bold.italic", 
#         color="DarkBlue",segment.color = NA)+
geom_label(
  data = y$df_mean.ord,
  aes(x = x, y = y, label=Group), 
  label.padding = unit(0.15, "lines"),label.size = 0.4,
)+guides(
  color=guide_legend(title="Sites"))+ylab("CCA2")+xlab("CCA1")



#a4
fdat <- ggplot2::fortify(vares_cca4_veg)
map1<- map[match(rownames(spp_hell4), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()

y<-ggordiplots::gg_ordiplot(vares_cca4_veg, groups = map1$Site, ellipse =  F, 
                            spiders = T, hull = F, plot = F )
z <- y$plot
a4<-z+theme_linedraw()+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  theme_linedraw()+
  scale_color_viridis_d(option ="turbo" )+#color of points +
  ggrepel::geom_text_repel(data = fdat %>% filter(
    Score=="biplot") %>%  
      mutate(CCA1=CCA1*0.2, CCA2=CCA2*0.2) ,
    aes(x=CCA1, y=CCA2, label= Label),
    color="red", segment.color = NA)+
  #geom_segment(data = fdat %>% filter(Score=="biplot") %>% 
  #         mutate(CCA1=CCA1*5, CCA2=CCA2*5) ,
  #     aes(x=0, y=0, xend=CCA1, yend=CCA2), 
  #   arrow=arrow(length=unit(0.15,"cm")),
  # size= 0.5,color="red")+
  theme_linedraw() +
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  #  ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
  #                             inner_join(taxonomy, 
  #                                       by = c("Label"="Feature.ID")) %>%
  #                           mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
  #                          filter(Label %in% varscho)  %>% 
  #                         mutate(tax=case_when(
  #                          tax=="pseudograminearum" ~ "Fusarium",
  #                         tax=="higginsianum" ~ "Colletotrichum",
  #                        tax=="oryzae" ~ "Aspergillus",
  #                       tax=="Colletitruchum fioriniae" ~ "Colletitruchum",
  #                      tax=="Aspergillus bombycis" ~ "Aspergillus",
#                     tax=="Aspergillus glaucus" ~ "Aspergillus",
#                    tax=="Aspergillus oryzae" ~ "Aspergillus",
#                   tax=="Laccaria bicolor" ~ "Laccaria",
#                  tax=="[Candida] auris" ~ "Clavispora",
#                 tax=="Penicillium arizonense" ~ "Penicillium",
#                TRUE ~ as.character(tax))) %>% 
#             mutate(CCA1=CCA1*scale, CCA2=CCA2*scale) ,
#          aes(x=CCA1, y=CCA2, label= tax), fontface="bold.italic", 
#         color="DarkBlue",segment.color = NA)+
geom_label(
  data = y$df_mean.ord,
  aes(x = x, y = y, label=Group), 
  label.padding = unit(0.15, "lines"),label.size = 0.4,
)+guides(
  color=guide_legend(title="Sites"))+ylab("CCA2")+xlab("CCA1")




option2<-plot_grid(a1+theme(aspect.ratio =6/10, legend.position = "none"),
                   a2+theme(aspect.ratio =6/10, legend.position = "none"),
                   a3+theme(aspect.ratio =6/10, legend.position = "none"),
                   a4+theme(aspect.ratio =6/10, legend.position = "none"), nrow = 2, ncol = 2)
ggsave("option2_cca.png",width = 10, height =6, dpi = 300, plot = option2, device = "png")
