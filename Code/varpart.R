library(tidyverse)
library(qiime2R)
library(readxl)
metadata_secas<- read_excel("Data/Metadatos.xlsx", sheet = "secas-marzo")
fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
fq_secas2<- read.csv("Data/fisicoq-la.csv")
meta_fq_secas<- metadata_secas %>% full_join(fq_secas) 
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")
env<-meta_fq_secas_all 

vars<- c("MO", "Ca", "LIMO", "ARCILLA", "WHC", "Fe", "K", "Mg",
         "pH", "Cu", "P", "N", "moisture")

envs<- env %>% column_to_rownames(var = "SampleID") %>%  dplyr::select(vars)

table_single<-read.delim("Data/table_micop_single.txt") 
table_paired<-read.delim("Data/table_micop_paired.txt") 
table_kraken<-read.delim("Data/table_fungi_again.txt", 
                         skip = 1,
                         row.names = 1, check.names = F) %>% dplyr::select_if(
                           is.numeric)%>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "id_sequence") %>% separate(
                             ., "id_sequence", c(
                               "kraken", "fungi", 
                               "id_metagenome", "report", "bracken"), 
                             sep = "_") %>% dplyr::select(
                               -kraken, -fungi, -report, -bracken) %>% full_join(
                                 metadata) %>% dplyr::select(-id_sequence:-Transecto,-Names, -id_metagenome, -Sites, -id_new, -id_fisicoq) %>% column_to_rownames(var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)

table_qiime2<-data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
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

table_qiime2_ord<- table_qiime2[,match(env$SampleID, colnames(table_qiime2))]
table_kraken_ord<- table_kraken[,match(env$SampleID, colnames(table_kraken))]
table_single_ord<- table_single[,match(env$SampleID, colnames(table_single))]
table_paired_ord<- table_paired[,match(env$SampleID, colnames(table_paired))]


table_qiime2_t<-data.frame(t(table_qiime2_ord), check.names = F)
table_single_t<-data.frame(t(table_single_ord), check.names = F)
table_paired_t<-data.frame(t(table_paired_ord), check.names = F)
table_kraken_t<-data.frame(t(table_kraken_ord), check.names = F)



veget<- read_excel("Data/vegetacion.xlsx")%>% dplyr::select(
-DAP, -CopaNS, -CopaEW,-volumenamdera, -distanciacopas, -diamcopa, -basal_area )

veg<-veget %>% separate(
  Especie, c("Genus", "Specie"), remove =F ) %>%   mutate(
    Genus=case_when(
      Genus=="Quecus"~ "Quercus",
      TRUE ~ as.character(Genus))) %>% mutate(
        Type=case_when(
          Genus=="Pinus"~ "Conifer",
          Genus=="Abies"~ "Conifer",
          Genus=="Juniperus"~ "Conifer",
          Genus=="Quercus"~ "Broadleaf",
          Genus=="Alnus"~ "Broadleaf",
          Genus=="Prunus"~ "Broadleaf",
          Genus=="Salix"~ "Broadleaf",
          Genus=="Arbutus"~ "Broadleaf") ) %>% mutate(
            Dominante=case_when(
              Genus=="Quercus"~ "Broadleaf",
              Genus=="Arbutus"~ "other Broadleaf",
              Genus=="Pinus"~ "Pinus",
              Genus=="Abies"~ "Abies",
              Genus=="Juniperus"~ "other Conifer",
              Genus=="Alnus"~ "other Broadleaf",
              Genus=="Alnus"~ "other Broadleaf",
              Genus=="Prunus"~ "other Broadleaf",
              Genus=="Salix"~ "other Broadleaf")  )

##proportions
prop_n<-veg %>% group_by(SampleID,Type) %>% dplyr::count()
prop_total<- veg %>% group_by(SampleID) %>% dplyr::count() %>% dplyr::rename(total=n)
prop<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Type", values_from = "prop")  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0)) %>% dplyr::rename(prop_Conif="Conifer",prop_Broadleaf="Broadleaf")

prop_n<-veg %>% group_by(SampleID,Genus) %>% dplyr::count()
prop_total<- veg %>% group_by(SampleID) %>% dplyr::count() %>% dplyr::rename(total=n)
prop_gen<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Genus", values_from = "prop")  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0))%>% 
  rename_with(~str_c("prop_", .), .cols = -SampleID)

# mean data
means<-veg %>% group_by(SampleID) %>% summarise_if(is.numeric, ~mean(., na.rm = TRUE))%>% 
  rename_with(~str_c("total_mean_", .), .cols = -SampleID)
veg_group<- veg %>% group_by(SampleID,Genus)%>% summarise_if(  is.numeric, ~mean(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")

means_genus<-function(y){veg_group %>% dplyr::select(
  SampleID, Genus,y) %>%  pivot_wider(
    . , names_from = "Genus", values_from = y)  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("mean_",y, "_"), .), .cols = -SampleID)}

vect<-c( "Height",   "coverage"  )
mean.list.genus <- vector("list")


for(i in vect) {
  mean.list.genus [[i]] <- means_genus(y = i)
}

mean.df.genus<-   mean.list.genus [[1]] %>% inner_join(
  mean.list.genus [[2]])

data_prop <- prop %>% inner_join(prop_gen) %>% column_to_rownames(
  var = "SampleID") 

data_veg_mean<- means %>% inner_join(mean.df.genus) %>% dplyr::select_at(
  vars(-contains("sites"))) %>%   dplyr::select_at(vars(-contains("volumen")))%>% column_to_rownames(
    var = "SampleID")

data_propt<- data_prop %>% log_norm() #%>% as.matrix() %>%scale(scale = T, center = T) 

vegeta= merge(data_propt, data_veg_mean, by=0)%>%dplyr::rename(
  SampleID="Row.names") %>% column_to_rownames(
    var = "SampleID")  %>% as.matrix() %>%#scale(
     # scale = T, center = T) #%>% 
as.data.frame() %>% rownames_to_column(
        var = "SampleID") %>% dplyr::select(
          -prop_Conif, -prop_Quercus,
          -mean_Height_Abies, -mean_Height_Arbutus,
          -mean_Height_Juniperus,-mean_Height_Pinus,
          -mean_Height_Quercus,-mean_Height_Prunus,
          -mean_Height_Alnus, -mean_Height_Salix,
          -mean_coverage_Prunus,  -mean_coverage_Salix,
          -mean_coverage_Alnus, -mean_coverage_Abies  )%>% 
  select_all(~str_replace(., "prop_", ""))  %>% 
  select_all(~str_replace(., "mean_", ""))  %>% 
  select_all(~str_replace(., "coverage_", "cov_")) %>% mutate(SampleID=case_when(
    SampleID=="P6S1T1" ~"P6S2T1",
    SampleID=="P6S1T2" ~"P6S2T2",
    SampleID=="P6S1T3" ~"P6S2T3",
    SampleID=="P6S2T1" ~"P6S1T1",
    SampleID=="P6S2T2" ~"P6S1T2",
    SampleID=="P6S2T3" ~"P6S1T3",
    TRUE~as.character(SampleID))) %>% arrange(SampleID) %>% column_to_rownames(var = "SampleID")

coords<- read_csv("Data/coord.csv")

coords_mat<- coords%>% mutate(P=paste0("P",pol),
                              S=paste0("S", Sitio),
                              T=paste0("T", Transecto)) %>%
  unite("SampleID",     P:T, sep = "") %>%
  dplyr::select(SampleID, Longitude, Latitude) %>% 
  column_to_rownames( var = "SampleID") %>% as.matrix() 

coords_pcnm <- as.data.frame(scores(pcnm(dist(distance_complete))))


#library(geosphere)
#distance<- distm(coords_mat)/1000
#colnames(distance)<- rownames(coords_mat)
#rownames(distance)<- rownames(coords_mat)
#distance_complete<- distance
#distance[upper.tri(distance)] <- NA 
#library(reshape2)
#distance_dm<-melt(as.matrix(distance), varnames = c(
 # "SampleID.x", "SampleID.y")) %>% drop_na() %>% filter(!value==0)

env.z <- decostand(envs, method = "standardize")
veg.z <- decostand(vegeta, method = "standardize")
spe.hel <- decostand(table_qiime2_t, method = "rclr")

spe.part.all <- varpart(table_qiime2_t, env.z, veg.z, transfo="hell")
spe.part.all <- varpart(table_single_t, env.z, veg.z, transfo="hell")
spe.part.all <- varpart(table_paired_t, env.z, veg.z, transfo="hell")
spe.part.all <- varpart(table_kraken_t, env.z, veg.z, transfo="hell")

summary(spe.part.all)
plot(spe.part.all, digits = 1, Xnames = c("Environmental", "Vegetation"),
     bg = c( 'tomato', 'green'))
