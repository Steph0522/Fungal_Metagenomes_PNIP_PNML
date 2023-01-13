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

data_spet<- spe_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_gent<- gen_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_spet<- tipo_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_propt<- data_prop %>% log_norm()

#load data of metagenomes
spp.16S1=read.delim("Data/table_micop_single.txt") 
spp.16S1<- spp.16S1[,match(rownames(data_spet), colnames(spp.16S1))]
spp.16S2=read.delim("Data/table_micop_paired.txt") 
spp.16S2<- spp.16S2[,match(rownames(data_spet), colnames(spp.16S2))]
spp.16S3=read.delim("Data/table_kraken.txt")
spp.16S3<- spp.16S3[,match(rownames(data_spet), colnames(spp.16S3))]
spp.16S4=read.delim("Data/table_qiime2.txt")
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
data_spet_dist<-vegdist(data_spet, method = "bray")
data_16S1_dist<- vegdist(spp.16S1t, method = "bray")
data_16S2_dist<- vegdist(spp.16S2t, method = "bray")
data_16S3_dist<- vegdist(spp.16S3t, method = "bray")
data_16S4_dist<- vegdist(spp.16S4t, method = "bray")

#data_spet_dist<- vegdist(spe_data, method = "horn")
#data_16S1_dist<- vegdist(spp.16S1, method = "horn")
#data_16S2_dist<- vegdist(spp.16S2, method = "horn")
#data_16S3_dist<- vegdist(spp.16S3, method = "horn")
#data_16S4_dist<- vegdist(spp.16S4, method = "horn")


vegan::mantel(data_spet_dist,data_16S1_dist)
vegan::mantel(data_spet_dist,data_16S2_dist)
vegan::mantel(data_spet_dist,data_16S3_dist)
vegan::mantel(data_spet_dist,data_16S4_dist)



#richness calcs
library(hillR)
q0_plant<- hill_taxa(spe_data, MARGIN = 1, q = 0) %>% as.data.frame(
  ) %>% rename(plant_richness=".") %>% rownames_to_column(var = "SampleID")
q1_plant<- hill_taxa(spe_data, MARGIN = 1, q = 1)%>% as.data.frame(
) %>% rename(plant_richness=".") %>% rownames_to_column(var = "SampleID")
q2_plant<- hill_taxa(spe_data, MARGIN = 1, q = 2)%>% as.data.frame(
) %>% rename(plant_richness=".") %>% rownames_to_column(var = "SampleID")

q0_sp1<- hill_taxa(spp.16S1, MARGIN = 1, q = 0)
q0_sp2<- hill_taxa(spp.16S2, MARGIN = 1, q = 0)
q0_sp3<- hill_taxa(spp.16S3, MARGIN = 1, q = 0)
q0_sp4<- hill_taxa(spp.16S4, MARGIN = 1, q = 0)

q1_sp1<- hill_taxa(spp.16S1, MARGIN = 1, q = 1)
q1_sp2<- hill_taxa(spp.16S2, MARGIN = 1, q = 1)
q1_sp3<- hill_taxa(spp.16S3, MARGIN = 1, q = 1)
q1_sp4<- hill_taxa(spp.16S4, MARGIN = 1, q = 1)

q2_sp1<- hill_taxa(spp.16S1, MARGIN = 1, q = 2)
q2_sp2<- hill_taxa(spp.16S2, MARGIN = 1, q = 2)
q2_sp3<- hill_taxa(spp.16S3, MARGIN = 1, q = 2)
q2_sp4<- hill_taxa(spp.16S4, MARGIN = 1, q = 2)

#q0

q0_data<- cbind(
  q0_sp1, q0_sp2, q0_sp3, q0_sp4) %>% as.data.frame() %>%rownames_to_column(
    var = "SampleID") %>% pivot_longer(
    cols = -SampleID, names_to = "Origin", values_to = "richness") %>% mutate(
      Method=case_when(
        Origin=="q0_sp1" ~ "Micop Single",
        Origin=="q0_sp2" ~ "Micop Paired",
        Origin=="q0_sp3" ~ "Kraken2",
        Origin=="q0_sp4" ~ "QIIME2")) %>% inner_join(q0_plant)

q0_data %>% ggplot(aes(x = plant_richness, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(~Method, scales = "free_y")

#q1

q1_data<- cbind(
  q1_sp1, q1_sp2, q1_sp3, q1_sp4) %>% as.data.frame() %>%rownames_to_column(
    var = "SampleID") %>% pivot_longer(
      cols = -SampleID, names_to = "Origin", values_to = "richness") %>% mutate(
        Method=case_when(
          Origin=="q1_sp1" ~ "Micop Single",
          Origin=="q1_sp2" ~ "Micop Paired",
          Origin=="q1_sp3" ~ "Kraken2",
          Origin=="q1_sp4" ~ "QIIME2")) %>% inner_join(q1_plant)

q1_data %>% ggplot(aes(x = plant_richness, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(~Method, scales = "free_y")

q1_data %>% ggplot(aes(x = plant_richness, y = richness,
                       color = Method, shape = Method)) +
  geom_point() + 
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE)

q1_1<- q1_sp1 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q1_1$plant_richness, q1_1$.)

q1_2<- q1_sp2 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q1_2$plant_richness, q1_2$.)

q1_3<- q1_sp3 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q1_3$plant_richness, q1_3$.)

q1_4<- q1_sp4 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q1_4$plant_richness, q1_4$.)

#q2

q2_data<- cbind(
  q2_sp1, q2_sp2, q2_sp3, q2_sp4) %>% as.data.frame() %>%rownames_to_column(
    var = "SampleID") %>% pivot_longer(
      cols = -SampleID, names_to = "Origin", values_to = "richness") %>% mutate(
        Method=case_when(
          Origin=="q2_sp1" ~ "Micop Single",
          Origin=="q2_sp2" ~ "Micop Paired",
          Origin=="q2_sp3" ~ "Kraken2",
          Origin=="q2_sp4" ~ "QIIME2")) %>% inner_join(q2_plant)

q2_data %>% ggplot(aes(x = plant_richness, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(~Method, scales = "free_y")


q2_1<- q2_sp1 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q2_1$plant_richness, q2_1$.)

q2_2<- q1_sp2 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q2_2$plant_richness, q2_2$.)

q2_3<- q2_sp3 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q2_3$plant_richness, q2_3$.)


q2_4<- q2_sp4 %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>%  inner_join(q1_plant)
cor.test(q2_4$plant_richness, q2_4$.)

#probando proporciones, con q=2 da mejor
a<- q2_sp1 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop)

cor.test(a$prop_conif, a$.)

#probando proporciones, generos arbutus y abies
a<- q2_sp1 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop_gen)

cor.test(a$prop_Arbutus, a$.)
a %>% ggplot(aes(x = ., y = prop_Arbutus))+
  geom_point()+ geom_smooth(method = "lm")

#resuminedo alpha

q1_data<- cbind(
  q2_sp1, q2_sp2, q2_sp3, q2_sp4) %>% as.data.frame() %>%rownames_to_column(var = "SampleID") %>% pivot_longer(
      cols = -SampleID, names_to = "Origin", values_to = "richness") %>% mutate(
        Method=case_when(
          Origin=="q2_sp1" ~ "Micop Single",
          Origin=="q2_sp2" ~ "Micop Paired",
          Origin=="q2_sp3" ~ "Kraken2",
          Origin=="q2_sp4" ~ "QIIME2")) %>% inner_join(prop)

q1_data %>% ggplot(aes(x = prop_conif, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(~Method, scales = "free_y")

q1_data %>% ggplot(aes(x = prop_conif, y = richness,
  color = Method, shape = Method)) +
  geom_point() + 
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE)


#beta diversity




