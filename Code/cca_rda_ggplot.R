#plot ccas
library(tidyverse)
library(ggvegan)
vars_choosing1<- fortify(vares_cca1_veg) %>% filter(Score=="species")%>%
inner_join(taxonomy_qiime2, by = c("Label"="Feature.ID")) %>% mutate(
tax= str_extract(Taxon, "[^_]+$"))#%>% filter(Label %in% varscho1)

vars_choosing2<- fortify(vares_cca2) %>% filter(Score=="species")%>%
  inner_join(taxonomy_single_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))#%>% filter(Label %in% varscho2)
vars_choosing3<- fortify(vares_cca3) %>% filter(Score=="species")%>%
  inner_join(taxonomy_paired_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))#%>% filter(Label %in% varscho3)
vars_choosing4<- fortify(vares_cca4) %>% filter(Score=="species")%>%
  inner_join(taxonomy_fungi, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))#%>% filter(Label %in% varscho4)


envschose<-c("P", "K", "Ca", "Mg", "P", "moisture", "WHC", "Silt")
varscho1<- c("62b57b7052f748d8888528f055a56166eab2eaf2","42c34c5cae8e20681d23626859950c1658f52183",
            "4dcb134727786b6b7cb9548cee1ac3171a946846")
varscho2<- c("Eukaryota__Ascomycota__Sordariomycetes__Glomerellales__Glomerellaceae__Colletotrichum__Colletotrichum fioriniae",
             "Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Aspergillus__Aspergillus bombycis",
             "Eukaryota__Basidiomycota__Agaricomycetes__Agaricales__Tricholomataceae__Laccaria__Laccaria bicolor")

varscho3<-c("Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Aspergillus__Aspergillus glaucus",
            "Eukaryota__Ascomycota__Saccharomycetes__Saccharomycetales__Metschnikowiaceae__Clavispora/Candida clade__[Candida] auris",
            "Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Penicillium__Penicillium arizonense")
varscho4<- c("101028","80884","5062")


cca_plot<- function(vares_cca, spp, varscho, envschose,taxonomy,scale, scal){
require(tidyverse)
require(readxl)
map<- read.csv("Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
require(ggvegan)
fdat <- ggplot2::fortify(vares_cca)
map1<- map[match(rownames(spp_hell1), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()

cca1_plot<-ggplot()+
  geom_point(data = fdat %>% filter(Score == "sites") %>% rename(
    SampleID=Label) %>% inner_join(map), aes(x =CCA1, y=CCA2), size=1, color="white")+
  ggforce::geom_mark_ellipse(data = fdat %>% filter(Score=="sites") %>% rename(
    SampleID=Label) %>% inner_join(map),aes(x =CCA1, y=CCA2,fill = Site))+
    scale_fill_viridis_d(option ="turbo")+#color of points 
 # scale_color_viridis_d(option ="turbo" )+
  ggrepel::geom_label_repel(data = fdat %>% filter(
    Score=="biplot")%>%mutate(Label=case_when(
                                Label=="ARCILLA"~"Silt", 
                                TRUE~as.character(Label))) %>% 
                              filter(Label %in% envschose) %>% 
                              mutate(CCA1=CCA1*scal, CCA2=CCA2*scal) ,
                            aes(x=CCA1, y=CCA2, label= Label),
                            segment.colour = "red",color="red", fontface="bold")+
 geom_segment(data = fdat %>% filter(Score=="biplot")%>%
               mutate(Label=case_when(
                Label=="ARCILLA"~"Silt", 
               TRUE~as.character(Label))) %>% 
            filter(Label %in% envschose) %>% 
           mutate(CCA1=CCA1*scal, CCA2=CCA2*scal) ,
        aes(x=0, y=0, xend=CCA1, yend=CCA2), 
       arrow=arrow(length=unit(0.15,"cm")),
      size= 0.5,color="red")+
  theme_linedraw() +
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 5), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  ggrepel::geom_label_repel(data=fdat %>% filter(Score == "sites")%>% rename(
    SampleID=Label) %>% inner_join(map) %>% group_by(Site) %>% 
      summarise_if(is.numeric, mean) %>% inner_join(map1_type),
    aes(x=CCA1, y=CCA2, label=Site, #, color=Type), 
        fontface = 'bold'))+
  ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
                             inner_join(taxonomy, 
                                        by = c("Label"="Feature.ID")) %>%
                             mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
                             filter(Label %in% varscho)  %>% 
                              mutate(tax=case_when(
                                tax=="pseudograminearum" ~ "Fusarium",
                                tax=="higginsianum" ~ "Colletotrichum",
                                tax=="oryzae" ~ "Aspergillus",
                                tax=="Colletitruchum fioriniae" ~ "Colletitruchum",
                                tax=="Aspergillus bombycis" ~ "Aspergillus",
                                tax=="Aspergillus glaucus" ~ "Aspergillus",
                                tax=="Aspergillus oryzae" ~ "Aspergillus",
                                tax=="Laccaria bicolor" ~ "Laccaria",
                                tax=="[Candida] auris" ~ "Clavispora",
                                tax=="Penicillium arizonense" ~ "Penicillium",
                                TRUE ~ as.character(tax))) %>% 
                             mutate(CCA1=CCA1*scale, CCA2=CCA2*scale) ,
                           aes(x=CCA1, y=CCA2, label= tax), fontface="bold.italic",
                           color="DarkBlue",segment.color = NA)

cca1_plot  }


cca_plot2<- function(vares_cca, spp, varscho, envschose,taxonomy, scale, scal){
  require(tidyverse)
  require(readxl)
  map<- read.csv("Data/coord.csv") %>% mutate_at(
    c(1,2,3,7), as.factor) %>% mutate(SampleID= 
                                        paste0("P",pol, "S", Sitio,"T", Transecto ))
  require(ggvegan)
  fdat <- ggplot2::fortify(vares_cca)
  map1<- map[match(rownames(spp_hell1), map$SampleID),]
  map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()
  
y<-ggordiplots::gg_ordiplot(vares_cca, groups = map1$Site, ellipse =  F, 
                            spiders = T, hull = F, plot = F )
z <- y$plot
a<-z+theme_linedraw()+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  theme_linedraw()+
 # scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
  scale_color_viridis_d(option ="turbo" )+#color of points +
  #scale_y_continuous(limits = c(-0.2,0.2))+
  #scale_x_continuous(limits = c(-0.2,0.2))+
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 5), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggrepel::geom_label_repel(data = fdat %>% filter(
    Score=="biplot")%>%mutate(Label=case_when(
      Label=="ARCILLA"~"Silt", 
      TRUE~as.character(Label))) %>% 
      filter(Label %in% envschose) %>% 
      mutate(CCA1=CCA1*scal, CCA2=CCA2*scal) ,
    aes(x=CCA1, y=CCA2, label= Label),
    color="red", fontface="bold", segment.color = NA)+
#  geom_segment(data = fdat %>% filter(Score=="biplot")%>%
 #                mutate(Label=case_when(
  #                 Label=="ARCILLA"~"Silt", 
   #                TRUE~as.character(Label))) %>% 
    #             filter(Label %in% envschose) %>% 
     #            mutate(CCA1=CCA1*scal, CCA2=CCA2*scal) ,
      #         aes(x=0, y=0, xend=CCA1, yend=CCA2), 
       #        arrow=arrow(length=unit(0.15,"cm")),
        #       size= 0.5,color="red")+
  theme_linedraw() +
  theme(axis.text = element_text(colour = "black", size = 5),
        axis.title = element_text(colour = "black", size = 5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 5), 
        legend.position = "right", 
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
  geom_hline(yintercept = 0, linetype = 2) +
  ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
                              inner_join(taxonomy, 
                                         by = c("Label"="Feature.ID")) %>%
                              mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
                              filter(Label %in% varscho)  %>% 
                              mutate(tax=case_when(
                                tax=="pseudograminearum" ~ "Fusarium",
                                tax=="higginsianum" ~ "Colletotrichum",
                                tax=="oryzae" ~ "Aspergillus",
                                tax=="Colletitruchum fioriniae" ~ "Colletitruchum",
                                tax=="Aspergillus bombycis" ~ "Aspergillus",
                                tax=="Aspergillus glaucus" ~ "Aspergillus",
                                tax=="Aspergillus oryzae" ~ "Aspergillus",
                                tax=="Laccaria bicolor" ~ "Laccaria",
                                tax=="[Candida] auris" ~ "Clavispora",
                                tax=="Penicillium arizonense" ~ "Penicillium",
                                TRUE ~ as.character(tax))) %>% 
                              mutate(CCA1=CCA1*scale, CCA2=CCA2*scale) ,
                            aes(x=CCA1, y=CCA2, label= tax), fontface="bold.italic", 
                            color="DarkBlue",segment.color = NA)+
  geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites"))
  #stat_ellipse(data = fdat %>% filter(Score=="sites") %>% rename(
   # SampleID=Label) %>% inner_join(map), aes(x =CCA1, y=CCA2, color=Type))
a}
cap_plot2<- function(vares_cap, spp, varscho, envschose,taxonomy, scale, scal){
  require(tidyverse)
  require(readxl)
  map<- read.csv("Data/coord.csv") %>% mutate_at(
    c(1,2,3,7), as.factor) %>% mutate(SampleID= 
                                        paste0("P",pol, "S", Sitio,"T", Transecto ))
  require(ggvegan)
  fdat <- ggplot2::fortify(vares_cap)
  map1<- map[match(rownames(spp_hell1), map$SampleID),]
  map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()
  
  y<-ggordiplots::gg_ordiplot(vares_cap, groups = map1$Site, ellipse =  F, 
                              spiders = T, hull = F, plot = F )
  z <- y$plot
  a<-z+theme_linedraw()+
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    # scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points +
    #scale_y_continuous(limits = c(-0.2,0.2))+
    #scale_x_continuous(limits = c(-0.2,0.2))+
    theme(axis.text = element_text(colour = "black", size = 5),
          axis.title = element_text(colour = "black", size = 5),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 5), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggrepel::geom_label_repel(data = fdat %>% filter(
      Score=="biplot")%>%mutate(Label=case_when(
        Label=="ARCILLA"~"Silt", 
        TRUE~as.character(Label))) %>% 
        filter(Label %in% envschose) %>% 
        mutate(CAP1=CAP1*scal, CAP2=CAP2*scal) ,
      aes(x=CAP1, y=CAP2, label= Label),
      color="red", fontface="bold", segment.color = NA)+
    #  geom_segment(data = fdat %>% filter(Score=="biplot")%>%
    #                mutate(Label=case_when(
    #                 Label=="ARCILLA"~"Silt", 
    #                TRUE~as.character(Label))) %>% 
    #             filter(Label %in% envschose) %>% 
    #            mutate(CAP1=CAP1*scal, CAP2=CAP2*scal) ,
    #         aes(x=0, y=0, xend=CAP1, yend=CAP2), 
    #        arrow=arrow(length=unit(0.15,"cm")),
    #       size= 0.5,color="red")+
    theme_linedraw() +
    theme(axis.text = element_text(colour = "black", size = 5),
          axis.title = element_text(colour = "black", size = 5),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 5), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    ggrepel::geom_label_repel(data = fdat %>% filter(Score=="species")%>%
                                inner_join(taxonomy, 
                                           by = c("Label"="Feature.ID")) %>%
                                mutate(tax= str_extract(Taxon, "[^_]+$")) %>% 
                                filter(Label %in% varscho)  %>% 
                                mutate(tax=case_when(
                                  tax=="pseudograminearum" ~ "Fusarium",
                                  tax=="higginsianum" ~ "Colletotrichum",
                                  tax=="oryzae" ~ "Aspergillus",
                                  tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
                                  tax=="Aspergillus bombycis" ~ "Aspergillus",
                                  tax=="Aspergillus glaucus" ~ "Aspergillus",
                                  tax=="Aspergillus oryzae" ~ "Aspergillus",
                                  tax=="Laccaria bicolor" ~ "Laccaria",
                                  tax=="[Candida] auris" ~ "Clavispora",
                                  tax=="Penicillium arizonense" ~ "Penicillium",
                                  TRUE ~ as.character(tax))) %>% 
                                mutate(CAP1=CAP1*scale, CAP2=CAP2*scale) ,
                              aes(x=CAP1, y=CAP2, label= tax), fontface="bold.italic", 
                              color="DarkBlue",segment.color = NA)+
    geom_label(
      data = y$df_mean.ord,
      aes(x = x, y = y, label=Group), 
      label.padding = unit(0.15, "lines"),label.size = 0.4,
    )+guides(
      color=guide_legend(title="Sites"))
  #stat_ellipse(data = fdat %>% filter(Score=="sites") %>% rename(
  # SampleID=Label) %>% inner_join(map), aes(x =CAP1, y=CAP2, color=Type))
  a}



#plotting with ggplot2
p1<-cca_plot2(vares_cca = vares_cca1_veg, spp = spp.16S1, varscho = varscho1,
              envschose = envschose, taxonomy = taxonomy_qiime2, scale = 1, scal=4.5)
p1<-cap_plot2(vares_cap = cap.env1, spp = spp.16S1, varscho = varscho1,
              envschose = envschose, taxonomy = taxonomy_qiime2, scale = 1, scal=0.8)
p1
p2<-cca_plot2(vares_cca = vares_cca2, spp = spp.16S2, varscho = varscho2,
              envschose = envschose, taxonomy = taxonomy_single_micop, scale = 0.8,
              scal = 1)
p2<-cap_plot2(vares_cap = cap.env2, spp = spp.16S2, varscho = varscho2,
              envschose = envschose, taxonomy = taxonomy_single_micop, scale = 1.55,
              scal = 0.8)
p2
p3<-cca_plot2(vares_cca = vares_cca3, spp = spp.16S3, varscho = varscho3,
              envschose = envschose, taxonomy = taxonomy_paired_micop, scale = 1.8, scal=0.7)
p3<-cap_plot2(vares_cap = cap.env3, spp = spp.16S3, varscho = varscho3,
              envschose = envschose, taxonomy = taxonomy_paired_micop, scale = 1.5,
              scal = 0.8)
p3
p4<-cca_plot2(vares_cca = vares_cca4, spp = spp.16S4, varscho = varscho4,
              envschose = envschose, taxonomy = taxonomy_fungi, scale = 1.3, scal = 0.2)
p4<-cap_plot2(vares_cap = cap.env4, spp = spp.16S4, varscho = varscho4,
              envschose = envschose, taxonomy = taxonomy_fungi, scale = 0.5,
              scal = 0.5)
p4
library(cowplot)
option1<-plot_grid(p1+theme(aspect.ratio =6/10, legend.position = "none"),
                   p2+theme(aspect.ratio =6/10, legend.position = "none"),
                   p3+theme(aspect.ratio =6/10, legend.position = "none"),
                   p4+theme(aspect.ratio =6/10, legend.position = "none"), nrow = 2, ncol = 2)

pp1<-cca_plot(vares_cca = vares_cca1, spp = spp.16S1, varscho = varscho1,
              envschose = envschose, taxonomy = taxonomy_qiime2, scale = 1.5, scal=12)
pp2<-cca_plot(vares_cca = vares_cca2, spp = spp.16S2, varscho = varscho2,
              envschose = envschose, taxonomy = taxonomy_single_micop, scale = 4, scal = 8)
pp3<-cca_plot(vares_cca = vares_cca3, spp = spp.16S3, varscho = varscho3,
              envschose = envschose, taxonomy = taxonomy_paired_micop, scale =20, scal = 8)
pp4<-cca_plot(vares_cca = vares_cca4, spp = spp.16S4, varscho = varscho4,
              envschose = envschose, taxonomy = taxonomy_fungi, scale = 100, scal = 12)

option2<-plot_grid(pp1+theme(aspect.ratio =6/10, legend.position = "none"),
                   pp2+theme(aspect.ratio =6/10, legend.position = "none"),
                   pp3+theme(aspect.ratio =6/10, legend.position = "none"),
                   pp4+theme(aspect.ratio =6/10, legend.position = "none"), nrow = 2, ncol = 2)
ggsave("option2_cca.png",width = 10, height =7, dpi = 300, plot = option2, device = "png")
