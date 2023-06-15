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
spp.16S2=read.delim("Data/table_micop_single.txt") 
spp.16S2<- spp.16S2[,match(rownames(data_spet), colnames(spp.16S2))]
spp.16S3=read.delim("Data/table_micop_paired.txt") 
spp.16S3<- spp.16S3[,match(rownames(data_spet), colnames(spp.16S3))]
spp.16S4=read.delim("Data/table_kraken.txt")
spp.16S4<- spp.16S4[,match(rownames(data_spet), colnames(spp.16S4))]
spp.16S1=read.delim("Data/table_qiime2.txt")
spp.16S1<- spp.16S1[,match(rownames(data_spet), colnames(spp.16S1))]

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

q0_sp1<- hill_taxa(spp.16S1, MARGIN = 1, q = 0) %>% as.data.frame()
q0_sp2<- hill_taxa(spp.16S2, MARGIN = 1, q = 0) %>% as.data.frame()
q0_sp3<- hill_taxa(spp.16S3, MARGIN = 1, q = 0) %>% as.data.frame()
q0_sp4<- hill_taxa(spp.16S4, MARGIN = 1, q = 0) %>% as.data.frame()

q1_sp1<- hill_taxa(spp.16S1, MARGIN = 1, q = 1) %>% as.data.frame()
q1_sp2<- hill_taxa(spp.16S2, MARGIN = 1, q = 1) %>% as.data.frame()
q1_sp3<- hill_taxa(spp.16S3, MARGIN = 1, q = 1) %>% as.data.frame()
q1_sp4<- hill_taxa(spp.16S4, MARGIN = 1, q = 1) %>% as.data.frame()

q2_sp1<- hill_taxa(spp.16S1, MARGIN = 1, q = 2)%>% as.data.frame()
q2_sp2<- hill_taxa(spp.16S2, MARGIN = 1, q = 2)%>% as.data.frame()
q2_sp3<- hill_taxa(spp.16S3, MARGIN = 1, q = 2)%>% as.data.frame()
q2_sp4<- hill_taxa(spp.16S4, MARGIN = 1, q = 2)%>% as.data.frame()

#q0
corq0_qiime2<- cor.test(q0_sp1$., q0_plant$plant_richness, method = "spearman")
qiime2_q0 <- data.frame(label = paste0("Pearson: r = ", signif(corq0_qiime2$estimate,3), 
                                          ", p-value = ", signif(corq0_qiime2$p.value, 3)))
corq0_single<- cor.test(q0_sp2$., q0_plant$plant_richness, method = "spearman")
single_q0 <- data.frame(label = paste0("Pearson: r = ", signif(corq0_single$estimate,3), 
                                       ", p-value = ", signif(corq0_single$p.value, 3)))

corq0_paired<- cor.test(q0_sp3$., q0_plant$plant_richness, method = "spearman")
paired_q0 <- data.frame(label = paste0("Pearson: r = ", signif(corq0_paired$estimate,3), 
                                       ", p-value = ", signif(corq0_paired$p.value, 3)))


corq0_kraken2<- cor.test(q0_sp4$., q0_plant$plant_richness, method = "spearman")
kraken2_q0 <- data.frame(label = paste0("Pearson: r = ", signif(corq0_kraken2$estimate,3), 
                                       ", p-value = ", signif(corq0_kraken2$p.value, 3)))
ann_text<-data.frame(richness=c( 500,140,200,500),
                     plant_richness=c(2.5,2.5,2.5,2.5),
                     Origin=c("QIIME2", "Single Micop", "Paired Micop", "Kraken2"),
                     label=c(qiime2_q0$label, 
                             single_q0$label,
                             paired_q0$label,
                             kraken2_q0$label)) 

names_col<-rownames(q0_sp1)


q0_data<- cbind(
  q0_sp1$., q0_sp2$., q0_sp3$., q0_sp4$.) %>% as.data.frame() %>% 
  rename("QIIME2"=V1, "Single Micop"=V2, "Paired Micop"=V3, "Kraken2"=V4)
rownames(q0_data)<-names_col

q0_data<- q0_data%>% 
  rownames_to_column(
    var = "SampleID") %>% pivot_longer(
    cols = -SampleID, names_to = "Origin", values_to = "richness") %>% 
  inner_join(q0_plant)
q0_data$Origin<- factor(q0_data$Origin, 
                           levels=c("QIIME2","Single Micop","Paired Micop","Kraken2"))
ann_text$Origin<- factor(ann_text$Origin, 
                        levels=c("QIIME2","Single Micop","Paired Micop","Kraken2"))

q0_plot<-q0_data %>% ggplot(aes(x = plant_richness, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(
    ~Origin, scales = "free_y")+
  geom_text(data = ann_text,label=ann_text$label, size=3)+
  xlab("Plant richness")+ylab("Fungal richness")

ggsave("q0_plot.png",width =7.5, height =6.5, dpi = 300,
       plot = q0_plot, device = "png")
#q1
corq1_qiime2<- cor.test(q1_sp1$., q1_plant$plant_richness, method = "spearman")
qiime2_q1 <- data.frame(label = paste0("Pearson: r = ", signif(corq1_qiime2$estimate,3), 
                                       ", p-value = ", signif(corq1_qiime2$p.value, 3)))
corq1_single<- cor.test(q1_sp2$., q1_plant$plant_richness, method = "spearman")
single_q1 <- data.frame(label = paste0("Pearson: r = ", signif(corq1_single$estimate,3), 
                                       ", p-value = ", signif(corq1_single$p.value, 3)))

corq1_paired<- cor.test(q1_sp3$., q1_plant$plant_richness, method = "spearman")
paired_q1 <- data.frame(label = paste0("Pearson: r = ", signif(corq1_paired$estimate,3), 
                                       ", p-value = ", signif(corq1_paired$p.value, 3)))


corq1_kraken2<- cor.test(q1_sp4$., q1_plant$plant_richness, method = "spearman")
kraken2_q1 <- data.frame(label = paste0("Pearson: r = ", signif(corq1_kraken2$estimate,3), 
                                        ", p-value = ", signif(corq1_kraken2$p.value, 3)))
ann_text<-data.frame(richness=c( 300,70,11,11),
                     plant_richness=c(2.5,2.5,2.5,2.5),
                     Origin=c("QIIME2", "Single Micop", "Paired Micop", "Kraken2"),
                     label=c(qiime2_q1$label, 
                             single_q1$label,
                             paired_q1$label,
                             kraken2_q1$label)) 

names_col<-rownames(q1_sp1)


q1_data<- cbind(
  q1_sp1$., q1_sp2$., q1_sp3$., q1_sp4$.) %>% as.data.frame() %>% 
  rename("QIIME2"=V1, "Single Micop"=V2, "Paired Micop"=V3, "Kraken2"=V4)
rownames(q1_data)<-names_col

q1_data<- q1_data%>% 
  rownames_to_column(
    var = "SampleID") %>% pivot_longer(
      cols = -SampleID, names_to = "Origin", values_to = "richness") %>% 
  inner_join(q1_plant)
q1_data$Origin<- factor(q1_data$Origin, 
                        levels=c("QIIME2","Single Micop","Paired Micop","Kraken2"))
ann_text$Origin<- factor(ann_text$Origin, 
                         levels=c("QIIME2","Single Micop","Paired Micop","Kraken2"))

q1_plot<-q1_data %>% ggplot(aes(x = plant_richness, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(
    ~Origin, scales = "free_y")+
  geom_text(data = ann_text,label=ann_text$label, size=3)+
  xlab("Plant q=1")+ylab("Fungal q=1")

ggsave("q1_plot.png",width =7.5, height =6.5, dpi = 300,
       plot = q1_plot, device = "png")

#q2
corq2_qiime2<- cor.test(q2_sp1$., q2_plant$plant_richness, method = "spearman")
qiime2_q2 <- data.frame(label = paste0("Pearson: r = ", signif(corq2_qiime2$estimate,3), 
                                       ", p-value = ", signif(corq2_qiime2$p.value, 3)))
corq2_single<- cor.test(q2_sp2$., q2_plant$plant_richness, method = "spearman")
single_q2 <- data.frame(label = paste0("Pearson: r = ", signif(corq2_single$estimate,3), 
                                       ", p-value = ", signif(corq2_single$p.value, 3)))

corq2_paired<- cor.test(q2_sp3$., q2_plant$plant_richness, method = "spearman")
paired_q2 <- data.frame(label = paste0("Pearson: r = ", signif(corq2_paired$estimate,3), 
                                       ", p-value = ", signif(corq2_paired$p.value, 3)))


corq2_kraken2<- cor.test(q2_sp4$., q2_plant$plant_richness, method = "spearman")
kraken2_q2 <- data.frame(label = paste0("Pearson: r = ", signif(corq2_kraken2$estimate,3), 
                                        ", p-value = ", signif(corq2_kraken2$p.value, 3)))
ann_text<-data.frame(richness=c( 300,70,11,11),
                     plant_richness=c(2.5,2.5,2.5,2.5),
                     Origin=c("QIIME2", "Single Micop", "Paired Micop", "Kraken2"),
                     label=c(qiime2_q2$label, 
                             single_q2$label,
                             paired_q2$label,
                             kraken2_q2$label)) 

names_col<-rownames(q2_sp1)


q2_data<- cbind(
  q2_sp1$., q2_sp2$., q2_sp3$., q2_sp4$.) %>% as.data.frame() %>% 
  rename("QIIME2"=V1, "Single Micop"=V2, "Paired Micop"=V3, "Kraken2"=V4)
rownames(q2_data)<-names_col

q2_data<- q2_data%>% 
  rownames_to_column(
    var = "SampleID") %>% pivot_longer(
      cols = -SampleID, names_to = "Origin", values_to = "richness") %>% 
  inner_join(q2_plant)
q2_data$Origin<- factor(q2_data$Origin, 
                        levels=c("QIIME2","Single Micop","Paired Micop","Kraken2"))
ann_text$Origin<- factor(ann_text$Origin, 
                         levels=c("QIIME2","Single Micop","Paired Micop","Kraken2"))

q2_plot<-q2_data %>% ggplot(aes(x = plant_richness, y = richness))+
  geom_point()+ geom_smooth(method = "lm")+facet_wrap(
    ~Origin, scales = "free_y")+
  geom_text(data = ann_text,label=ann_text$label, size=3)+
  xlab("Plant q=2")+ylab("Fungal q=2")

ggsave("q2_plot.png",width =7.5, height =6.5, dpi = 300,
       plot = q2_plot, device = "png")



#probando proporciones, con q=2 da mejor
a<- q2_sp1 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop)
b<- q2_sp2 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop)
c<- q2_sp3 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop)
d<- q2_sp4 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop)

cor.test(a$prop_latif, a$.)
cor.test(b$prop_latif, b$.)
cor.test(c$prop_latif, c$.)
cor.test(d$prop_latif, d$.)


#probando proporciones, generos arbutus y abies
#qiime2
a<- q2_sp1 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop_gen)

cor.test(a$., a$prop_Abies)
cor.test(a$., a$prop_Juniperus)

a %>% ggplot(aes(x = ., y = prop_Abies))+
  geom_point()+ geom_smooth(method = "lm")

a %>% ggplot(aes(x = ., y = prop_Juniperus))+
  geom_point()+ geom_smooth(method = "lm")

#single

b<- q2_sp2 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop_gen)
cor.test(b$., b$prop_Abies)
cor.test(b$., b$prop_Arbutus)
cor.test(b$., b$prop_Quercus)

b %>% ggplot(aes(x = ., y = prop_Abies))+
  geom_point()+ geom_smooth(method = "lm")

b %>% ggplot(aes(x = ., y = prop_Arbutus))+
  geom_point()+ geom_smooth(method = "lm")


b %>% ggplot(aes(x = ., y = prop_Quercus))+
  geom_point()+ geom_smooth(method = "lm")

#paired
c<- q2_sp3 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop_gen)

cor.test(c$., c$prop_Abies)
cor.test(c$., c$prop_Pinus)


c %>% ggplot(aes(x = ., y = prop_Abies))+
  geom_point()+ geom_smooth(method = "lm")


c%>% ggplot(aes(x = ., y = prop_Pinus))+
  geom_point()+ geom_smooth(method = "lm")

#kraken

d<- q2_sp4 %>% as.data.frame() %>%rownames_to_column(
  var = "SampleID") %>%  inner_join(prop_gen)
cor.test(d$., d$prop_Salix)

d%>% ggplot(aes(x = ., y = prop_Salix))+
  geom_point()+ geom_smooth(method = "lm")



# matrices de correlación-------------------------------------------------------------------------


library(corrr)

data_all<- cbind(q0_sp1, q0_sp2, q0_sp3, q0_sp4)
colnames(data_all)<- c("QIIME2", "Single_Micop", "Paired_Micop", "Kraken2")

data_all<- data_all %>% rownames_to_column(var = "SampleID") %>% 
  inner_join(prop_gen)

res.cor<- correlate(data_all[-1], diagonal = 1)

mat<-res.cor %>% focus(QIIME2, Single_Micop, Paired_Micop) %>% 
  column_to_rownames(var = "term") %>% slice(-1)

library(ComplexHeatmap)
Heatmap(mat, cluster_rows = F, cluster_columns = F,
        column_title = "q0",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", mat[i, j]),
                    x, y, gp = gpar(fontsize = 10))
        })
