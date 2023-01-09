library(tidyverse)
library(readxl)
veget<- read_excel("Data/vegetacion.xlsx")

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

#proporcion generos dominantes
prop_n<-veg %>% group_by(SampleID,Dominante) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop_domin<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Dominante", values_from = "prop")  %>% mutate_if(is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
  rename_with(~str_c("prop_domin", .), .cols = -SampleID)

#proporcion generos
prop_n<-veg %>% group_by(SampleID,Genus) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop_gen<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Genus", values_from = "prop")  %>% mutate_if(is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
  rename_with(~str_c("prop_", .), .cols = -SampleID)

#proporcion familias
prop_n<-veg %>% group_by(SampleID,Familia) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop_fam<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Familia", values_from = "prop")  %>% mutate_if(is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
  rename_with(~str_c("prop_", .), .cols = -SampleID)

#datos totales promedios y medianas
means<-veg %>% group_by(SampleID) %>% summarise_if(is.numeric, ~mean(., na.rm = TRUE))%>% 
  rename_with(~str_c("total_mean_", .), .cols = -SampleID)
medians<- veg %>% group_by(SampleID)%>% summarise_if(is.numeric, ~median(., na.rm = TRUE))%>% 
  rename_with(~str_c("total_median_", .), .cols = -SampleID)

#datos totales promedios y medianas por coniferas y latifoliadas
veg_group<- veg %>% group_by(SampleID,Tipo)%>% summarise_if(  is.numeric,
                                                              ~mean(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")
veg_group2<- veg %>% group_by(SampleID,Tipo)%>% summarise_if(  is.numeric,
                                                               ~median(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")

means_tipo<-function(y){veg_group %>% dplyr::select(SampleID, Tipo,y) %>%  pivot_wider(
    . , names_from = "Tipo", values_from = y)  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
  rename_with(~str_c(paste0("mean_",y, "_"), .), .cols = -SampleID)}

median_tipo<-function(y){veg_group2 %>% dplyr::select(
  SampleID, Tipo,y) %>%  pivot_wider(
  . , names_from = "Tipo", values_from = y)  %>% mutate_if(
    is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("median_",y, "_"), .), .cols = -SampleID)}

vect<-c( "DAP", "Altura",  "CopaNS",  "distanciacopas", "diamcopa", "coberturacopa" ,
         "areabasal"    ,  "volumenamdera" )
mean.list <- vector("list")
median.list <- vector("list")

for(i in vect) {
  mean.list [[i]] <- means_tipo(y = i)
}
for(i in vect) {
  median.list [[i]] <- median_tipo(y = i)
}


mean.df<- mean.list[[1]] %>% inner_join(
  mean.list[[2]]) %>% inner_join(
    mean.list[[3]])%>% inner_join(
      mean.list[[4]])%>% inner_join(
        mean.list[[5]])%>% inner_join(
          mean.list[[6]])%>% inner_join(
            mean.list[[7]])%>% inner_join(
              mean.list[[8]])

median.df<- median.list[[1]] %>% inner_join(
  median.list[[2]]) %>% inner_join(
    median.list[[3]])%>% inner_join(
      median.list[[4]])%>% inner_join(
        median.list[[5]])%>% inner_join(
          median.list[[6]])%>% inner_join(
            median.list[[7]])%>% inner_join(
              median.list[[8]])

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
#datos totales promedios y medianas por generos dominantes
veg_group<- veg %>% group_by(SampleID,Dominante)%>% summarise_if(  is.numeric,
                                                                   ~mean(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")
veg_group2<- veg %>% group_by(SampleID,Dominante)%>% summarise_if(  is.numeric,
                                                                    ~median(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")

means_Dominante<-function(y){veg_group %>% dplyr::select(SampleID, Dominante,y) %>%  pivot_wider(
  . , names_from = "Dominante", values_from = y)  %>% mutate_if(
    is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("mean_Dominante_",y, "_"), .), .cols = -SampleID)}

median_Dominante<-function(y){veg_group2 %>% dplyr::select(
  SampleID, Dominante,y) %>%  pivot_wider(
    . , names_from = "Dominante", values_from = y)  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("median_Dominante_",y, "_"), .), .cols = -SampleID)}

vect<-c( "DAP", "Altura",  "CopaNS",  "distanciacopas", "diamcopa", "coberturacopa" ,
         "areabasal"    ,  "volumenamdera" )
mean.list.Dominante <- vector("list")
median.list.Dominante <- vector("list")

for(i in vect) {
  mean.list.Dominante [[i]] <- means_Dominante(y = i)
}
for(i in vect) {
  median.list.Dominante [[i]] <- median_Dominante(y = i)
}



mean.df.Dominante<-   mean.list.Dominante [[1]] %>% inner_join(
  mean.list.Dominante [[2]]) %>% inner_join(
    mean.list.Dominante [[3]])%>% inner_join(
      mean.list.Dominante [[4]])%>% inner_join(
        mean.list.Dominante [[5]])%>% inner_join(
          mean.list.Dominante [[6]])%>% inner_join(
            mean.list.Dominante [[7]])%>% inner_join(
              mean.list.Dominante [[8]])

median.df.Dominante<- median.list.Dominante[[1]] %>% inner_join(
  median.list.Dominante[[2]]) %>% inner_join(
    median.list.Dominante[[3]])%>% inner_join(
      median.list.Dominante[[4]])%>% inner_join(
        median.list.Dominante[[5]])%>% inner_join(
          median.list.Dominante[[6]])%>% inner_join(
            median.list.Dominante[[7]])%>% inner_join(
              median.list.Dominante[[8]])

#datos totales promedios y medianas por  Familias
veg_group<- veg %>% group_by(SampleID,Familia)%>% summarise_if(  is.numeric,
                                                                 ~mean(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")
veg_group2<- veg %>% group_by(SampleID,Familia)%>% summarise_if(  is.numeric,
                                                                  ~median(., na.rm = TRUE)) %>% dplyr::select(-Sites) #%>% column_to_rownames(var = "SampleID")

means_Familia<-function(y){veg_group %>% dplyr::select(SampleID, Familia,y) %>%  pivot_wider(
  . , names_from = "Familia", values_from = y)  %>% mutate_if(
    is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("mean_",y, "_"), .), .cols = -SampleID)}

median_Familia<-function(y){veg_group2 %>% dplyr::select(
  SampleID, Familia,y) %>%  pivot_wider(
    . , names_from = "Familia", values_from = y)  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("median_",y, "_"), .), .cols = -SampleID)}

vect<-c( "DAP", "Altura",  "CopaNS",  "distanciacopas", "diamcopa", "coberturacopa" ,
         "areabasal"    ,  "volumenamdera" )
mean.list.Familia <- vector("list")
median.list.Familia <- vector("list")

for(i in vect) {
  mean.list.Familia [[i]] <- means_Familia(y = i)
}
for(i in vect) {
  median.list.Familia [[i]] <- median_Familia(y = i)
}



mean.df.Familia<-   mean.list.Familia [[1]] %>% inner_join(
  mean.list.Familia [[2]]) %>% inner_join(
    mean.list.Familia [[3]])%>% inner_join(
      mean.list.Familia [[4]])%>% inner_join(
        mean.list.Familia [[5]])%>% inner_join(
          mean.list.Familia [[6]])%>% inner_join(
            mean.list.Familia [[7]])%>% inner_join(
              mean.list.Familia [[8]])

median.df.Familia<- median.list.Familia[[1]] %>% inner_join(
  median.list.Familia[[2]]) %>% inner_join(
    median.list.Familia[[3]])%>% inner_join(
      median.list.Familia[[4]])%>% inner_join(
        median.list.Familia[[5]])%>% inner_join(
          median.list.Familia[[6]])%>% inner_join(
            median.list.Familia[[7]])%>% inner_join(
              median.list.Familia[[8]])

#uniendo datos

data_joined<- prop %>% inner_join(prop_domin) %>% inner_join(
  prop_gen) %>% inner_join(prop_fam) %>% inner_join(
    means) %>%inner_join(medians) %>%  inner_join(mean.df) %>% inner_join(
      median.df) %>% inner_join(mean.df.genus) %>% inner_join(
        median.df.genus) %>% inner_join(mean.df.Dominante) %>% inner_join(
          median.df.Dominante) %>% inner_join(mean.df.Familia) %>% inner_join(median.df.Familia)
# checking collinearity
library(picante)
library(Rmisc)

corss<- cor(data_joined[-1], method = "pearson")
library(corrplot)
#par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
corrplot(corss, type="upper")
#png(filename = "cor.png", width = 830, height = 630, res = 100)
corrplot(corss,  method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
#por genero
data_gen<- data_joined %>% dplyr::select(SampleID,prop_conif, prop_latif, 
                                         starts_with("total_mean"), -total_mean_Sites,
                                         starts_with("total_median"), -total_median_Sites,
                                         34:65,66:193 ) %>% dplyr::select(
                                           -contains("volumen")) %>% dplyr::select(
                                             -contains("Sites"))
#colnames(data_gen)<- c(1:171)
corss<- cor(data_gen[-1], method = "pearson")
corrplot(corss,  method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
# RDA separately for 16S and 18S data
library(vegan)
library(qiime2R)
#setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")
spp.16S1=read.delim("Data/table_micop_single.txt") 
spp.16S2=read.delim("Data/table_micop_paired.txt") 

spp.16S3=read.delim("Data/table_fungi_again.txt", 
                    skip = 1, row.names = 1, check.names = F) %>% dplyr::select_if(
                      is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                        var = "id_sequence") %>% separate(
                          ., "id_sequence", 
                          c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                          sep = "_") %>% dplyr::select(
                            -kraken, -fungi, -report, -bracken) %>% full_join(
                            metadata) %>% dplyr::select(
                              -id_sequence:-Transecto, -id_metagenome, 
                              -Sites, -id_new,   -id_fisicoq, -Pol) %>% column_to_rownames(
                              var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)
spp.16S4=data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
                    check.names = F) %>% t() %>% as.data.frame(
                    ) %>% rownames_to_column(
                      var = "SampleID") %>% separate(
                        SampleID, c(
                          "id_metagenome", "R", "unmap", "Paired"), 
                        sep = "_")%>% inner_join(
                          metadata) %>% dplyr::select(
                            -id_metagenome:-Paired, 
                            -id_sequence:-id_fisicoq, -Sites, 
                            -Pol) %>% column_to_rownames(
                              var="SampleID") %>% t(
                              ) %>% as.data.frame() %>% mutate_all(as.numeric)



env.16S=data_gen

spp.16S1=data.frame(t(spp.16S1), check.names = F)
spp.16S2=data.frame(t(spp.16S2), check.names = F)
spp.16S3=data.frame(t(spp.16S3), check.names = F)
spp.16S4=data.frame(t(spp.16S4), check.names = F)


#select variables that did not show collinearity

#kraken
mm3=env.16S%>%  inner_join(spp.16S3 %>%rownames_to_column(var = "SampleID") )
#mm3=na.omit(mm3)
dim(mm3)
head(mm3)
env.16S3=mm3[,c(85:95,98,101,103)]
env.16S3=data_gen %>% column_to_rownames(var = "SampleID")
dim(env.16S3)
colnames(env.16S3)
spp.16S3=mm3[,1:73]%>% column_to_rownames(var = "SampleID")


#single
mm1=env.16S%>%  inner_join(spp.16S1 %>%rownames_to_column(var = "SampleID") )
#mm1=na.omit(mm1)
dim(mm1)
head(mm1)
#env.16S1=mm1[,c(221:231,234,235, 237:239)]
#env.16S1=mm1[,c(221:231,234, 237,239)]
env.16S1=data_gen %>% column_to_rownames(var = "SampleID")
colnames(env.16S1)
spp.16S1=mm1[,1:211]%>% column_to_rownames(var = "SampleID")

#paired
mm2=env.16S%>%  inner_join(spp.16S2 %>%rownames_to_column(var = "SampleID") )
dim(mm2)
head(mm2)
env.16S2=mm2[,c(245:255,258, 261,263)]
env.16S2=data_gen %>% column_to_rownames(var = "SampleID")
spp.16S=mm2[,1:235]%>% column_to_rownames(var = "SampleID")

#qiime2
mm4=env.16S%>%  inner_join(spp.16S4 %>%rownames_to_column(var = "SampleID") )
#mm4=na.omit(mm4)
dim(mm4)
head(mm4)
#env.16S4=mm4[,c(3151:3161, 3164,  3167,3169)]
env.16S4=data_gen %>% column_to_rownames(var = "SampleID")
dim(env.16S4)
spp.16S4=mm4[,1:3139]%>% column_to_rownames(var = "SampleID")


#dcas
#dca=decorana(spp.16S1)
#plot(dca) # RDA is okay to use

#transoforming enviromental and species data

spp.16S_hell1=decostand(spp.16S1, "hell") # Hellinger transformation
spp.16S_hell2=decostand(spp.16S2, "hell")
spp.16S_hell3=decostand(spp.16S3, "hell")
spp.16S_hell4=decostand(spp.16S4, "hell")

#write_csv(env.16S_st1, "/home/yendi/Desktop/viendo.csv")

env.16S_st1=data.frame(scale(env.16S1, scale=T, center=F)) # standardize env. data
env.16S_st2=data.frame(scale(env.16S2, scale=T, center=F))
env.16S_st3=data.frame(scale(env.16S3, scale=T, center=F))
env.16S_st4=data.frame(scale(env.16S4, scale=T, center=F))

#PC1=rda(spp.16S_hell1~., env.16S_st1)
#plot(PC1, display = "sites" ,type = c("text"))
#plot(PC1, type="n")
#text(PC1, dis="cn", lwd=1, cex=1)
#text(PC1, dis="sites", type=c("text"), cex=0.8)
#points(PC1, dis="sites", lwd=1, cex=0.5)
#PC1=rda(spp.16S1~., env.16S1)
#plot(PC1, type="n")
#text(PC1, dis="cn", lwd=2.5, cex=1.2)
#text(PC1, dis="sites", type=c("text"), cex=0.8)
#points(PC1, dis="sites", lwd=2.5)

#quartz.save("RDA0_16S.pdf", type="pdf")

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


vars_caps_filter<-vars_caps %>% filter(`Pr(>F)`<0.05)
vars_capsj_filter<-vars_capsj %>% filter(`Pr(>F)`<0.05)

write.csv(vars_capsj_filter, "/home/yendi/Desktop/vars_capsj.csv", row.names = T)


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

write.csv(envs, "/home/yendi/Desktop/envs.csv")

env.16S=data_gen 

df<- env.16S %>% dplyr::select_at(
  vars(-contains("mean")))%>% dplyr::select_at(
    vars(-contains("Alnus")))
dfs=data.frame(df[1],scale(df[,2:72], center = T, scale = T)) %>% 
  column_to_rownames(var = "SampleID")
dfs_dist<- dist(dfs, method = "euclidean")

nut.dist<- dfs_dist %>% as.matrix()

nut.dist[upper.tri(nut.dist)] <- NA 

