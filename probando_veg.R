library(tidyverse)
library(readxl)

#load data and format

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

#datos totales promedios y medianas
means<-veg %>% group_by(SampleID) %>% summarise_if(is.numeric, ~mean(., na.rm = TRUE))%>% 
  rename_with(~str_c("total_mean_", .), .cols = -SampleID)
medians<- veg %>% group_by(SampleID)%>% summarise_if(is.numeric, ~median(., na.rm = TRUE))%>% 
  rename_with(~str_c("total_median_", .), .cols = -SampleID)

#datos totales promedios y medianas por generos
veg_group<- veg %>% group_by(SampleID,Genus)%>% summarise_if(  is.numeric,
                                                               ~mean(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")
veg_group2<- veg %>% group_by(SampleID,Genus)%>% summarise_if(  is.numeric,
                                                                ~median(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")

means_genus<-function(y){veg_group %>% dplyr::select(SampleID, Genus,y) %>%  pivot_wider(
  . , names_from = "Genus", values_from = y)  %>% mutate_if(
    is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("mean_",y, "_"), .), .cols = -SampleID)}

median_genus<-function(y){veg_group2 %>% dplyr::select(
  SampleID, Genus,y) %>%  pivot_wider(
    . , names_from = "Genus", values_from = y)  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("median_",y, "_"), .), .cols = -SampleID)}

vect<-c( "DAP", "Altura",  "CopaNS",  "distanciacopas", "diamcopa", "coberturacopa" ,
         "areabasal"    ,  "volumenamdera" )
mean.list.genus <- vector("list")
median.list.genus <- vector("list")

for(i in vect) {
  mean.list.genus [[i]] <- means_genus(y = i)
}
for(i in vect) {
  median.list.genus [[i]] <- median_genus(y = i)
}


mean.df.genus<-   mean.list.genus [[1]] %>% inner_join(
  mean.list.genus [[2]]) %>% inner_join(
    mean.list.genus [[3]])%>% inner_join(
      mean.list.genus [[4]])%>% inner_join(
        mean.list.genus [[5]])%>% inner_join(
          mean.list.genus [[6]])%>% inner_join(
            mean.list.genus [[7]])%>% inner_join(
              mean.list.genus [[8]])

median.df.genus<- median.list.genus[[1]] %>% inner_join(
  median.list.genus[[2]]) %>% inner_join(
    median.list.genus[[3]])%>% inner_join(
      median.list.genus[[4]])%>% inner_join(
        median.list.genus[[5]])%>% inner_join(
          median.list.genus[[6]])%>% inner_join(
            median.list.genus[[7]])%>% inner_join(
              median.list.genus[[8]])

#uniendo datos

data_prop <- prop %>% inner_join(prop_gen) %>% column_to_rownames(
  var = "SampleID") 
data_veg_mean<- means %>% inner_join(mean.df.genus) %>% dplyr::select_at(
  vars(-contains("sites"))) %>%   dplyr::select_at(vars(-contains("volumen")))%>% column_to_rownames(
    var = "SampleID")
data_veg_mean_transform<- data_veg_mean  %>% as.matrix() %>%scale(scale = T, center = T) 
#data_veg_mean_transforms<- data_veg_mean%>% as.matrix() %>%sqrt()

data_veg_medians<- means %>% inner_join(median.df.genus) %>% dplyr::select_at(
  vars(-contains("sites"))) %>%   dplyr::select_at(vars(-contains("volumen")))%>% column_to_rownames(
    var = "SampleID")
data_veg_median_transform<- data_veg_medians%>% as.matrix() %>%scale(scale = T, center = T) 
#data_veg_median_transforms<- data_veg_medians %>% as.matrix() %>%sqrt()

data_spet<- spe_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_gent<- gen_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_tipot<- tipo_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_propt<- data_prop %>% log_norm()

#data_me<-merge(data_prop, data_veg_mean, by = 0) %>% column_to_rownames(var = "Row.names")
#data_met<- merge(data_prop, data_veg_mean_transform, by = 0)%>% column_to_rownames(var = "Row.names")
#data_mets<- merge(data_prop, data_veg_mean_transforms, by = 0)%>% column_to_rownames(var = "Row.names")
#data_med<-merge(data_prop, data_veg_medians, by = 0)%>% column_to_rownames(var = "Row.names")
#data_medt<- merge(data_prop, data_veg_median_transform, by = 0)%>% column_to_rownames(var = "Row.names")
#data_medts<- merge(data_prop, data_veg_median_transforms, by = 0)%>% column_to_rownames(var = "Row.names")

#load data of metagenomes
spp.16S1=read.delim("Data/table_micop_single.txt") 
spp.16S2=read.delim("Data/table_micop_paired.txt") 
spp.16S3=read.delim("Data/table_kraken.txt")
spp.16S4=read.delim("Data/table_qiime2.txt")
spp.16S1=data.frame(t(spp.16S1), check.names = F)
spp.16S2=data.frame(t(spp.16S2), check.names = F)
spp.16S3=data.frame(t(spp.16S3), check.names = F)
spp.16S4=data.frame(t(spp.16S4), check.names = F)

#environmental vars
env.16S= merge(data_propt, data_veg_mean_transform, by=0)%>%rename(SampleID="Row.names")

#kraken
mm3=env.16S%>%  inner_join(spp.16S3 %>%rownames_to_column(var = "SampleID") )
dim(mm3)
head(mm3)
env.16S3=mm3[1:74] %>% column_to_rownames(var = "SampleID")
dim(env.16S3)
colnames(env.16S3)
spp.16S3=mm3[c(1,75:146)]%>% column_to_rownames(var = "SampleID")
colnames(spp.16S3)

#single
mm1=env.16S%>%  inner_join(spp.16S1 %>%rownames_to_column(var = "SampleID") )
dim(mm1)
head(mm1)
env.16S1=mm1[1:74] %>% column_to_rownames(var = "SampleID")
colnames(env.16S1)
spp.16S1=mm1[c(1,75:284)]%>% column_to_rownames(var = "SampleID")

#paired
mm2=env.16S%>%  inner_join(spp.16S2 %>%rownames_to_column(var = "SampleID") )
dim(mm2)
head(mm2)
env.16S2=mm2[1:74] %>% column_to_rownames(var = "SampleID")
colnames(env.16S2)
spp.16S2=mm2[c(1,75:308)]%>% column_to_rownames(var = "SampleID")

#qiime2
mm4=env.16S%>%  inner_join(spp.16S4 %>%rownames_to_column(var = "SampleID") )
dim(mm4)
head(mm4)
env.16S4=mm2[1:74] %>% column_to_rownames(var = "SampleID")
dim(env.16S4)
spp.16S2=mm4[c(1,75:3212)]%>% column_to_rownames(var = "SampleID")

#transoforming enviromental and species data
library(vegan)
spp.16S_hell1=decostand(spp.16S1, "hell") # Hellinger transformation
spp.16S_hell2=decostand(spp.16S2, "hell")
spp.16S_hell3=decostand(spp.16S3, "hell")
spp.16S_hell4=decostand(spp.16S4, "hell")

#env.16S_st1=data.frame(scale(env.16S1, scale=T, center=F)) # standardize env. data
#env.16S_st2=data.frame(scale(env.16S2, scale=T, center=F))
#env.16S_st3=data.frame(scale(env.16S3, scale=T, center=F))
#env.16S_st4=data.frame(scale(env.16S4, scale=T, center=F))

env.16S_st1=env.16S1
env.16S_st2=env.16S2
env.16S_st3=env.16S3
env.16S_st4=env.16S4

# Forward selection procedure
#accounting for abundance
cap.env=capscale(spp.16S_hell1~.,env.16S_st1,distance = "bray")
mod0.env=capscale(spp.16S_hell1~1,env.16S_st1,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
cap1<-step.env$anova %>% mutate(method="Single Micop") # These are the selected structuring variables that are included in the subsequent analyses 

cap.env=capscale(spp.16S_hell2~.,env.16S_st2,distance = "bray")
mod0.env=capscale(spp.16S_hell2~1,env.16S_st2,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
cap2<-step.env$anova %>% mutate(method="Paired Micop") 

cap.env=capscale(spp.16S_hell3~.,env.16S_st3,distance = "bray")
mod0.env=capscale(spp.16S_hell3~1,env.16S_st3,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
cap3<-step.env$anova %>% mutate(method="Kraken2") 


cap.env=capscale(spp.16S_hell4~.,env.16S_st4,distance = "bray")
mod0.env=capscale(spp.16S_hell4~1,env.16S_st4,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
cap4<-step.env$anova %>% mutate(method="QIIME2") 

#accounting for richness

cap.env=capscale(spp.16S_hell1~.,env.16S_st1,distance = "jaccard")
mod0.env=capscale(spp.16S_hell1~1,env.16S_st1,distance = "jaccard")
step.env=ordistep(mod0.env,scope = formula(cap.env))
capj1<-step.env$anova %>% mutate(method="Single Micop") # These are the selected structuring variables that are included in the subsequent analyses 

cap.env=capscale(spp.16S_hell2~.,env.16S_st2,distance = "jaccard")
mod0.env=capscale(spp.16S_hell2~1,env.16S_st2,distance = "jaccard")
step.env=ordistep(mod0.env,scope = formula(cap.env))
capj2<-step.env$anova %>% mutate(method="Paired Micop") 

cap.env=capscale(spp.16S_hell3~.,env.16S_st3,distance = "jaccard")
mod0.env=capscale(spp.16S_hell3~1,env.16S_st3,distance = "jaccard")
step.env=ordistep(mod0.env,scope = formula(cap.env))
capj3<-step.env$anova %>% mutate(method="Kraken2") 


cap.env=capscale(spp.16S_hell4~.,env.16S_st4,distance = "jaccard")
mod0.env=capscale(spp.16S_hell4~1,env.16S_st4,distance = "jaccard")
step.env=ordistep(mod0.env,scope = formula(cap.env))
cap4j<-step.env$anova %>% mutate(method="QIIME2") 

vars_caps<- rbind(cap1, cap2, cap3, cap4)
vars_capsj<- rbind(capj1, capj2, capj3, cap4j)

rownames(vars_caps)#vars to select
rownames(vars_capsj)#vars to select



vars_caps_filter<-vars_caps %>% filter(`Pr(>F)`<0.05)
vars_capsj_filter<-vars_capsj %>% filter(`Pr(>F)`<0.05)

#write.csv(vars_capsj_filter, "/home/yendi/Desktop/vars_capsj.csv", row.names = T)


#selecting using cca
vares_cca1 <- cca(spp.16S_hell1 ~ . ,
                  data=env.16S_st1)
env1<-envfit(vares_cca1 ~ . ,
             data=env.16S_st1)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Single Micop"=".")

vares_cca2 <- cca(spp.16S_hell2 ~ .,
                  data=env.16S_st2)
env2<-envfit(vares_cca2 ~ .,
             data=env.16S_st2)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Paired Micop"=".")

vares_cca3 <- cca(spp.16S_hell3 ~ . ,
                  data=env.16S_st3)
env3<-envfit(vares_cca3 ~ . ,
             data=env.16S_st3)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Kraken2"=".")

vares_cca4 <- cca(spp.16S_hell4 ~ . ,
                  data=env.16S_st4)
env4<-envfit(vares_cca4 ~ .,
             data=env.16S_st4)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("QIIME2"=".")

envs<- cbind(env4, env1, env2, env3)
write.csv(envs, "/home/yendi/Desktop/envs_t.csv")

#env.16S=env.16S
#vars_chose<-  c("prop_Alnus", "prop_Quercus",               
#"total_mean_DAP" , "mean_Altura_Abies", "mean_Altura_Pinus",
#"prop_Arbutus"  , "mean_distanciacopas_Quercus", "mean_distanciacopas_Pinus",
#"mean_CopaNS_Salix" ,"mean_areabasal_Abies" ,"mean_areabasal_Salix", "mean_Altura_Prunus")

#vars_chose<- c( "mean_CopaNS_Abies","mean_Altura_Pinus","mean_areabasal_Salix",
 #              "mean_DAP_Pinus","prop_Arbutus",  
  #             "mean_Altura_Abies", "mean_CopaNS_Quercus",    
   #            "mean_areabasal_Abies", "mean_CopaNS_Salix"   )

vars_chose<- c("mean_CopaNS_Abies", "mean_Altura_Pinus", "mean_areabasal_Salix",     
               "mean_distanciacopas_Pinus", 
             "prop_Arbutus","prop_conif",               
              "mean_Altura_Juniperus" , "prop_Arbutus",             
            "mean_coberturacopa_Abies" , "mean_CopaNS_Abies" ,
          "mean_Altura_Quercus",    
         "mean_areabasal_Abies"  ,    "prop_Salix",               
          "mean_Altura_Prunus" )

#c<- read.csv("/home/yendi/Desktop/envs_t.csv", row.names = 1)

df<- env.16S %>% dplyr::select(SampleID,vars_chose)
#df<- env.16S %>% dplyr::select(SampleID,rownames(c))

dfss<- merge(env.16S, dfsA, by = "SampleID") %>% dplyr::select(
  SampleID,vars_chose,P,K,Ca,Mg,moisture,WHC, LIMO)


#df1<- merge(env.16S1, dfs, by = 0)
#df1<- merge(env.16S1, dfsA, by = 0) %>% dplyr::select(Row.names,
 # prop_Alnus, prop_Quercus, total_mean_DAP, P:LIMO) %>% column_to_rownames(var = "Row.names") %>% 
  #dist(., method = "euclidean")
#df1<- merge(env.16S1, dfs, by = 0) %>% dplyr::select(
 # Row.names,mean_CopaNS_Abies,mean_Altura_Pinus,
  #mean_areabasal_Salix,mean_DAP_Pinus, P:LIMO) %>% column_to_rownames(var = "Row.names") %>% 
  #dist(., method = "euclidean")

df1<-env.16S1 %>% dplyr::select(mean_CopaNS_Abies,mean_Altura_Pinus,
                                mean_areabasal_Salix,mean_distanciacopas_Pinus)  %>% 
  dist(., method = "euclidean")
#df1<-env.16S1 %>% dplyr::select(rownames(c1))  %>% 
 # dist(., method = "euclidean")
df2<-env.16S2 %>% dplyr::select(prop_Arbutus,prop_conif, mean_Altura_Juniperus,
                               mean_Altura_Juniperus,mean_coberturacopa_Abies,
mean_CopaNS_Abies ,       mean_Altura_Quercus     )  %>% 
dist(., method = "euclidean")

#df2<- merge(env.16S2, dfs, by=0) %>% dplyr::select(Row.names,prop_Arbutus,
 #   mean_Altura_Pinus, mean_distanciacopas_Quercus,
  #                               mean_Altura_Abies, P:LIMO)%>% column_to_rownames(var = "Row.names") %>% 
  #dist(., method = "euclidean")
#df2<- merge(env.16S2, dfs, by=0) %>% dplyr::select(rownames(c2)) %>% 
 # dist(., method = "euclidean")

df3<- env.16S3 %>% dplyr::select(mean_CopaNS_Abies) %>% 
 dist(., method = "euclidean")
df4<- env.16S3 %>% dplyr::select(mean_areabasal_Abies,
prop_Salix, mean_Altura_Prunus  ) %>% 
 dist(., method = "euclidean")
#dfs3<- env.16S3 %>% dplyr::select("prop_Juniperus"  ,  "mean_DAP_Pinus" ,
 #                                 "mean_Altura_Juniperus"  , "mean_CopaNS_Pinus" ,
  #                                "mean_CopaNS_Salix" , "mean_distanciacopas_Pinus",
   #                               "mean_diamcopa_Pinus", "mean_areabasal_Pinus"  )  %>% 
  #dist(., method = "euclidean")

dfs=df %>% 
  column_to_rownames(var = "SampleID")
dfs_dist<- dist(dfs, method = "euclidean")

dfsA=dfss %>% 
 column_to_rownames(var = "SampleID")
dfs_distA<- dist(dfsA, method = "euclidean")



nut.dist<- dfs_dist %>% as.matrix()

#nut.dist[upper.tri(nut.dist)] <- NA 

map<- read.csv("Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)

otu <- list(table_single_micop, table_paired_micop, table_qiime2, table_fungi)
otu_match<- lapply(otu, otu.match) # matching to map
otu_single<- lapply(otu_match, otu.single) #remove singletons
otu_norm<- lapply(otu_single, otu.norm)#Normalize
bc.dist<- lapply(otu_norm, beta_div_dist)#Calculate Bray-Curtis dissimilarities
bc.dist2<- lapply(otu_single, beta_div_dist_hill, q=1)


allm1<-mantel(bc.dist2[[1]], dfs_distA)
allm2<-mantel(bc.dist2[[2]], dfs_distA)
allm3<-mantel(bc.dist2[[3]], dfs_distA)
allm4<-mantel(bc.dist2[[4]], dfs_distA)

mantel(bc.dist2[[1]], dfs_dist)
mantel(bc.dist2[[2]], dfs_dist)
mantel(bc.dist2[[3]], dfs_dist)
mantel(bc.dist2[[4]], dfs_dist)

callm1<-cor.test(bc.dist2[[1]], as.matrix(dfs_distA))
callm2<-cor.test(bc.dist2[[2]], as.matrix(dfs_distA))
callm3<-cor.test(bc.dist2[[3]], as.matrix(dfs_distA))
callm4<-cor.test(bc.dist2[[4]], as.matrix(dfs_distA))


m1<-mantel(bc.dist2[[1]], df1)
M1<-mantel(bc.dist2[[1]], df1)

m2<-mantel(bc.dist2[[2]], df2)
M2<-mantel(bc.dist2[[2]], df2)


m3<-mantel(bc.dist2[[3]], df3)
m4<-mantel(bc.dist2[[4]], df4)

m1b<-mantel(bc.dist[[1]], df1)
m2b<-mantel(bc.dist[[2]], df2)
m3b<-mantel(bc.dist[[3]], df3)
m4b<-mantel(bc.dist[[4]], df4)



nut.dist.tidy <- nut.dist %>% 
  melt(as.matrix(df1), varnames = c(
    "SampleID.x", "SampleID.y"), value.name = "EucDist") %>% drop_na(
          ) %>% filter(!EucDist==0) %>% 
  filter(!is.na(EucDist)) %>% 
  filter(SampleID.x != SampleID.y) 

  bc.dist[[1]][upper.tri(bc.dist[[1]])] <- NA 
  bc.dist.tidy<-bc.dist[[1]] %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
  b.dist.filt <- bc.dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value) %>% 
    full_join(nut.dist.tidy) %>% 
    mutate(Similarity = 1 - Distance)  

  b.dist.filt %>% 
    ggplot(aes(EucDist, Similarity)) +
    geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
    geom_smooth(method = "lm", color = "black", se = F) +
    xlab("Enviromental distance") +
    ylab("Bray-curtis similarity") +
    ylim(.2, max.sim) +
    xlim(0,8)+
    theme_linedraw()+theme(legend.position = "none", 
                           axis.text = element_text(size = 12),
                           strip.text = element_text(size = 12, face = "bold.italic"),
                           panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                           colour = "#E5E8E8"), 
                           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                           colour = "#E5E8E8"))  
  