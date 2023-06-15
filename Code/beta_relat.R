library(tidyverse)
library(readxl)

veget<- read_excel("Data/vegetacion.xlsx")


#species data
spe_data<- veget %>% group_by(SampleID, Especie) %>% count() %>% pivot_wider(
  names_from = "Especie", values_from = "n")%>% mutate_if(
    is.numeric, ~replace(., is.na(.), 0)) %>% column_to_rownames(var = "SampleID")
gen_data<-veget %>% separate(
  Especie, c("Genus", "Specie"), remove =F ) %>%   mutate(
    Genus=case_when(
      Genus=="Quecus"~ "Quercus",
      TRUE ~ as.character(Genus))) %>% group_by(SampleID, Genus) %>% count() %>% pivot_wider(
        names_from = "Genus", values_from = "n")%>% mutate_if(
          is.numeric, ~replace(., is.na(.), 0)) %>% column_to_rownames(var = "SampleID")

tipo_data<- veget %>% separate(
  Especie, c("Genus", "Specie"), remove =F ) %>%   mutate(
    Genus=case_when(
      Genus=="Quecus"~ "Quercus",
      TRUE ~ as.character(Genus))) %>% mutate(
        Tipo=case_when(
          Genus=="Pinus"~ "Coníferas",
          Genus=="Abies"~ "Coníferas",
          Genus=="Juniperus"~ "Coníferas",
          Genus=="Quercus"~ "Latifoliadas",
          Genus=="Alnus"~ "Latifoliadas",
          Genus=="Prunus"~ "Latifoliadas",
          Genus=="Salix"~ "Latifoliadas",
          Genus=="Arbutus"~ "Latifoliadas") ) %>% group_by(SampleID, Tipo) %>% count() %>% pivot_wider(
            names_from = "Tipo", values_from = "n")%>% mutate_if(
              is.numeric, ~replace(., is.na(.), 0)) %>% column_to_rownames(var = "SampleID")



#declare data

veg<-veget %>% separate(
  Especie, c("Genus", "Specie"), remove =F ) %>%   mutate(
    Genus=case_when(
      Genus=="Quecus"~ "Quercus",
      TRUE ~ as.character(Genus))) %>% mutate(
        Tipo=case_when(
          Genus=="Pinus"~ "Coníferas",
          Genus=="Abies"~ "Coníferas",
          Genus=="Juniperus"~ "Coníferas",
          Genus=="Quercus"~ "Latifoliadas",
          Genus=="Alnus"~ "Latifoliadas",
          Genus=="Prunus"~ "Latifoliadas",
          Genus=="Salix"~ "Latifoliadas",
          Genus=="Arbutus"~ "Latifoliadas") ) %>% mutate(
            Dominante=case_when(
              Genus=="Quercus"~ "Quercus",
              Genus=="Arbutus"~ "Arbutus",
              Genus=="Pinus"~ "Pinus",
              Genus=="Abies"~ "Abies",
              Genus=="Juniperus"~ "otras coníferas",
              Genus=="Alnus"~ "otras Latifoliadas",
              Genus=="Alnus"~ "otras Latifoliadas",
              Genus=="Prunus"~ "otras Latifoliadas",
              Genus=="Salix"~ "otras Latifoliadas")  )

#proporcion coniferas y latifoliadas
prop_n<-veg %>% group_by(SampleID,Tipo) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Tipo", values_from = "prop")  %>% mutate_if(is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0)) %>% rename(prop_conif="Coníferas",
                                                       prop_latif="Latifoliadas")

##proporcion generos
prop_n<-veg %>% group_by(SampleID,Genus) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop_gen<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Genus", values_from = "prop")  %>% mutate_if(is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
  rename_with(~str_c("prop_", .), .cols = -SampleID)

# data species

data_prop <- prop %>% inner_join(prop_gen) %>% column_to_rownames(
  var = "SampleID") 

data_spet<- spe_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_gent<- gen_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_spet<- tipo_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_propt<- data_prop %>% log_norm()

#load data of metagenomes
spp.16S1=data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
                     check.names = F) %>% t() %>% as.data.frame(
                     ) %>% rownames_to_column(
                       var = "SampleID") %>% separate(
                         SampleID, c(
                           "id_metagenome", "R", "unmap", "Paired"), 
                         sep = "_")%>% inner_join(
                           metadata) %>% dplyr::select(
                             -id_metagenome:-Paired, 
                             -id_sequence:-id_fisicoq, -Sites, 
                             -Names) %>% column_to_rownames(
                               var="SampleID") %>% t(
                               ) %>% as.data.frame() %>% mutate_all(as.numeric)
spp.16S1<- spp.16S1[,match(rownames(data_spet), colnames(spp.16S1))]

spp.16S2=read.delim("Data/table_micop_single.txt") 
spp.16S2<- spp.16S2[,match(rownames(data_spet), colnames(spp.16S2))]

spp.16S3=read.delim("Data/table_micop_paired.txt") 
spp.16S3<- spp.16S3[,match(rownames(data_spet), colnames(spp.16S3))]

spp.16S4=read.delim("Data/table_fungi_again.txt", 
                     skip = 1, row.names = 1, check.names = F) %>% dplyr::select_if(
                       is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                         var = "id_sequence") %>% separate(
                           ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                           sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                             metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names,  -id_fisicoq) %>% column_to_rownames(
                               var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)

spp.16S4<- spp.16S4[,match(rownames(data_spet), colnames(spp.16S4))]

spp.16S1=data.frame(t(spp.16S1), check.names = F) 
spp.16S2=data.frame(t(spp.16S2), check.names = F)
spp.16S3=data.frame(t(spp.16S3), check.names = F)
spp.16S4=data.frame(t(spp.16S4), check.names = F)

spp.16S1t=spp.16S1 %>%t() %>% rel_ab() %>% log_norm() %>% t()
spp.16S2t=spp.16S2 %>%t() %>%  rel_ab() %>% log_norm() %>% t()
spp.16S3t=spp.16S3 %>%t() %>%  rel_ab() %>% log_norm() %>% t()
spp.16S4t=spp.16S4 %>%t() %>%  rel_ab() %>% log_norm() %>% t()

#betadiv
library(vegan)
data_spet_dist<-vegdist(data_spet, method = "bray") %>% as.matrix()
data_16S1_dist<- vegdist(spp.16S1t, method = "bray")%>% as.matrix()
data_16S2_dist<- vegdist(spp.16S2t, method = "bray")%>% as.matrix()
data_16S3_dist<- vegdist(spp.16S3t, method = "bray")%>% as.matrix()
data_16S4_dist<- vegdist(spp.16S4t, method = "bray")%>% as.matrix()

#data_spet_dist<- vegdist(spe_data, method = "horn")
#data_16S1_dist<- vegdist(spp.16S1, method = "horn")
#data_16S2_dist<- vegdist(spp.16S2, method = "horn")
#data_16S3_dist<- vegdist(spp.16S3, method = "horn")
#data_16S4_dist<- vegdist(spp.16S4, method = "horn")


vegan::mantel(data_spet_dist,data_16S1_dist)
vegan::mantel(data_spet_dist,data_16S2_dist)
vegan::mantel(data_spet_dist,data_16S3_dist)
vegan::mantel(data_spet_dist,data_16S4_dist)



#beta diversity
library(reshape2)

data_spet_dist[upper.tri(data_spet_dist)] <- NA 
data_16S1_dist[upper.tri(data_16S1_dist)] <- NA 
data_16S2_dist[upper.tri(data_16S2_dist)] <- NA 
data_16S3_dist[upper.tri(data_16S3_dist)] <- NA 
data_16S4_dist[upper.tri(data_16S4_dist)] <- NA 


spet_dist<- data_spet_dist   %>% melt(., varnames = c(
  "SampleID.x", "SampleID.y")) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
  filter(!is.na(value)) %>% dplyr::rename(SpeDistance=value)

dist_16S1<- data_16S1_dist  %>% melt(., varnames = c(
  "SampleID.x", "SampleID.y")) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
  filter(!is.na(value)) %>% dplyr::rename(Distance=value)

dist_16S2<- data_16S2_dist  %>% melt(., varnames = c(
  "SampleID.x", "SampleID.y")) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
  filter(!is.na(value)) %>% dplyr::rename(Distance=value)

dist_16S3<- data_16S3_dist  %>% melt(., varnames = c(
  "SampleID.x", "SampleID.y")) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
  filter(!is.na(value)) %>% dplyr::rename(Distance=value)

dist_16S4<- data_16S4_dist  %>% melt(., varnames = c(
  "SampleID.x", "SampleID.y")) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
  filter(!is.na(value)) %>% dplyr::rename(Distance=value)

D1<- spet_dist %>% full_join(dist_16S1)
D2<- spet_dist %>% full_join(dist_16S2)
D3<- spet_dist %>% full_join(dist_16S3)
D4<- spet_dist %>% full_join(dist_16S4)


cor_D1<- cor.test(D1$Distance,D1$SpeDistance, method= "pearson", alternative = "two.sided") %>% tidy()
cor_D2<- cor.test(D2$Distance,D2$SpeDistance, method= "pearson", alternative = "two.sided") %>% tidy()
cor_D3<- cor.test(D3$Distance,D3$SpeDistance, method= "pearson", alternative = "two.sided") %>% tidy()
cor_D4<- cor.test(D4$Distance,D4$SpeDistance, method= "pearson", alternative = "two.sided") %>% tidy()


stats_D1 <- data.frame(label = paste0("correlation: r =", signif(cor_D1$estimate,2), ",",
                                       " p-value =", signif(cor_D1$p.value, 3)))
stats_D2 <- data.frame(label = paste0("correlation: r =", signif(cor_D2$estimate,2), ",",
                                      " p-value =", signif(cor_D2$p.value, 3)))
stats_D3 <- data.frame(label = paste0("correlation: r =", signif(cor_D3$estimate,2), ",",
                                      " p-value =", signif(cor_D3$p.value, 3)))
stats_D4 <- data.frame(label = paste0("correlation: r =", signif(cor_D4$estimate,2), ",",
                                      " p-value =", signif(cor_D4$p.value, 3)))
#PLOT
p1<-D1 %>% 
  ggplot(aes(SpeDistance, Distance)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = F) +
  xlab("Bray curtis fungal dissimilarity") +
  ylab("Bray curtis vegtation dissimilarity") +
  ylim(.2, 1) +
 # xlim(0,60)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.text = element_text(size = 12),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  ggtitle(stats_D1$label)+theme(plot.title = element_text(size = 10))

p2<-D2 %>% 
  ggplot(aes(SpeDistance, Distance)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = F) +
  xlab("Bray curtis fungal dissimilarity") +
  ylab("Bray curtis vegtation dissimilarity") +
  ylim(.2, 1) +
  # xlim(0,60)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.text = element_text(size = 12),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  ggtitle(stats_D2$label)+theme(plot.title = element_text(size = 10))


p3<-D3 %>% 
  ggplot(aes(SpeDistance, Distance)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = F) +
  xlab("Bray curtis fungal dissimilarity") +
  ylab("Bray curtis vegtation dissimilarity") +
  ylim(.2, 1) +
  # xlim(0,60)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.text = element_text(size = 12),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  ggtitle(stats_D3$label)+theme(plot.title = element_text(size = 10))

p4<-D4 %>% 
  ggplot(aes(SpeDistance, Distance)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = F) +
  xlab("Bray curtis fungal dissimilarity") +
  ylab("Bray curtis vegtation dissimilarity") +
  ylim(0, 1) +
  # xlim(0,60)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.text = element_text(size = 12),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  ggtitle(stats_D4$label)+theme(plot.title = element_text(size = 10))


ps<-plot_grid(p1,p2,p3,p4)

ggsave("Fig.beta_all.png",width =8, height = 7, dpi = 300, plot = ps, device = "png")

  
