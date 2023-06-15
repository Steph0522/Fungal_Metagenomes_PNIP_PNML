library(tidyverse)
library(readxl)
library(qiime2R)

#load data and format
map<- read.csv("Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)


veget<- read_excel("Data/vegetacion.xlsx") %>% dplyr::select(
  -DAP, -CopaNS, -CopaEW,-volumenamdera, -distanciacopas, -diamcopa, -basal_area )

#species data
spe_data<- veget%>%   mutate(
 Especie=case_when(
    Especie=="Quecus glabrescens"~ "Quercus glabrescens",
    TRUE ~ as.character(Especie))) %>% group_by(SampleID, Especie) %>% count() %>% pivot_wider(
  names_from = "Especie", values_from = "n")%>% mutate_if(
    is.numeric, ~replace(., is.na(.), 0)) %>% column_to_rownames(var = "SampleID")
gen_data<-veget %>% separate(
  Especie, c("Genus", "Specie"), remove =F ) %>%   mutate(
    Genus=case_when(
      Genus=="Quecus"~ "Quercus",
      TRUE ~ as.character(Genus))) %>% group_by(SampleID, Genus) %>% count() %>% pivot_wider(
        names_from = "Genus", values_from = "n")%>% mutate_if(
          is.numeric, ~replace(., is.na(.), 0)) %>% column_to_rownames(var = "SampleID")

Type_data<- veget %>% separate(
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
          Genus=="Arbutus"~ "Broadleaf") ) %>% group_by(SampleID, Type) %>% count() %>% pivot_wider(
            names_from = "Type", values_from = "n")%>% mutate_if(
              is.numeric, ~replace(., is.na(.), 0)) %>% column_to_rownames(var = "SampleID")



#declare data

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


#proporcion coniferas y Broadleaf
prop_n<-veg %>% group_by(SampleID,Type) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop<- prop_n %>% inner_join(prop_total)%>% mutate(
prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Type", values_from = "prop")  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0)) %>% rename(prop_conif="Conifer",
                                                       prop_broadleaf="Broadleaf")


##proporcion generos
prop_n<-veg %>% group_by(SampleID,Genus) %>% count()
prop_total<- veg %>% group_by(SampleID) %>% count() %>% rename(total=n)
prop_gen<- prop_n %>% inner_join(prop_total)%>% mutate(
  prop=n/total*100)%>% dplyr::select(-n,-total) %>%  pivot_wider(
    . , names_from = "Genus", values_from = "prop")  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
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

means_genus<-function(y){veg_group %>% dplyr::select(
  SampleID, Genus,y) %>%  pivot_wider(
  . , names_from = "Genus", values_from = y)  %>% mutate_if(
    is.numeric, ~round(.,digits = 2))  %>% mutate_if(
      is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("mean_",y, "_"), .), .cols = -SampleID)}

median_genus<-function(y){veg_group2 %>% dplyr::select(
  SampleID, Genus,y) %>%  pivot_wider(
    . , names_from = "Genus", values_from = y)  %>% mutate_if(
      is.numeric, ~round(.,digits = 2))  %>% mutate_if(
        is.numeric, ~replace(., is.na(.), 0))%>% 
    rename_with(~str_c(paste0("mean_",y, "_"), .), .cols = -SampleID)}

vect<-c( "DAP", "Height",  "CopaNS",  "distanciacopas", "diamcopa", "coverage" ,
         "basal_area"    ,  "volumenamdera" )

vect<-c( "Height",   "coverage"    )
mean.list.genus <- vector("list")
median.list.genus <- vector("list")

for(i in vect) {
  mean.list.genus [[i]] <- means_genus(y = i)
}
for(i in vect) {
  median.list.genus [[i]] <- median_genus(y = i)
}


mean.df.genus<-   mean.list.genus [[1]] %>% inner_join(
  mean.list.genus [[2]]) #%>% inner_join(
    #mean.list.genus [[3]])#%>% inner_join(
    #  mean.list.genus [[4]])%>% inner_join(
   #     mean.list.genus [[5]])#%>% inner_join(
      #    mean.list.genus [[6]])%>% inner_join(
       #     mean.list.genus [[7]])%>% inner_join(
        #      mean.list.genus [[8]])

median.df.genus<- median.list.genus[[1]] %>% inner_join(
  median.list.genus[[2]])# %>% inner_join(
  #  median.list.genus[[3]])#%>% inner_join(
      #median.list.genus[[4]])%>% inner_join(
       # median.list.genus[[5]])#%>% inner_join(
  #        median.list.genus[[6]])%>% inner_join(
 #           median.list.genus[[7]])%>% inner_join(
#              median.list.genus[[8]])

#uniendo datos

data_prop <- prop %>% inner_join(prop_gen) %>% column_to_rownames(
  var = "SampleID") 
data_veg_mean<- means %>% inner_join(mean.df.genus) %>% dplyr::select_at(
  vars(-contains("sites"))) %>%   dplyr::select_at(vars(-contains("volumen")))%>% column_to_rownames(
    var = "SampleID")
data_veg_mean_transform<- data_veg_mean  %>% as.matrix() %>%scale(scale = T, center = T) 
#data_veg_mean_transforms<- data_veg_mean%>% as.matrix() %>%sqrt()

data_veg_medians<- medians %>% inner_join(median.df.genus) %>% dplyr::select_at(
  vars(-contains("sites"))) %>%   dplyr::select_at(vars(-contains("volumen")))%>% column_to_rownames(
    var = "SampleID")
data_veg_median_transform<- data_veg_medians%>% as.matrix() %>%scale(scale = T, center = T) 
#data_veg_median_transforms<- data_veg_medians %>% as.matrix() %>%sqrt()

data_spet<- spe_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_gent<- gen_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_Typet<- Type_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
data_propt<- data_prop %>% log_norm() #%>% as.matrix() %>%scale(scale = T, center = T) 

#data_me<-merge(data_prop, data_veg_mean, by = 0) %>% column_to_rownames(var = "Row.names")
#data_met<- merge(data_prop, data_veg_mean_transform, by = 0)%>% column_to_rownames(var = "Row.names")
#data_mets<- merge(data_prop, data_veg_mean_transforms, by = 0)%>% column_to_rownames(var = "Row.names")
#data_med<-merge(data_prop, data_veg_medians, by = 0)%>% column_to_rownames(var = "Row.names")
#data_medt<- merge(data_prop, data_veg_median_transform, by = 0)%>% column_to_rownames(var = "Row.names")
#data_medts<- merge(data_prop, data_veg_median_transforms, by = 0)%>% column_to_rownames(var = "Row.names")

#load data of metagenomes
library(qiime2R)
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
map<-metadata
spp.16S2t=read.delim("Data/table_micop_single.txt") 
spp.16S3t=read.delim("Data/table_micop_paired.txt") 
spp.16S1t=data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
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
spp.16S1t=spp.16S1t[, match(colnames(spp.16S2t), colnames(spp.16S1t))]
spp.16S4t=read.delim("Data/table_fungi_again.txt", 
                    skip = 1, row.names = 1, check.names = F) %>% dplyr::select_if(
                      is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                        var = "id_sequence") %>% separate(
                          ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                          sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                            metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names,  -id_fisicoq) %>% column_to_rownames(
                              var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)


spp.16S1=data.frame(t(spp.16S1t), check.names = F)
spp.16S2=data.frame(t(spp.16S2t), check.names = F)
spp.16S3=data.frame(t(spp.16S3t), check.names = F)
spp.16S4=data.frame(t(spp.16S4t), check.names = F)

#environmental vars
env.16S= merge(data_propt, data_veg_mean, by=0)%>%rename(
  SampleID="Row.names") %>% column_to_rownames(
    var = "SampleID")  %>% as.matrix() %>%scale(
      scale = T, center = T) %>% as.data.frame() %>% rownames_to_column(
        var = "SampleID") %>% dplyr::select(
          -prop_conif, -prop_Quercus,
          -mean_Height_Abies, -mean_Height_Arbutus,
          -mean_Height_Juniperus,-mean_Height_Pinus,
          -mean_Height_Quercus,-mean_Height_Prunus,
          -mean_Height_Alnus, -mean_Height_Salix,
          -mean_coverage_Prunus,  -mean_coverage_Salix,
          -mean_coverage_Alnus, -mean_coverage_Abies  )
#env.16S= merge(data_prop, data_veg_mean, by=0)%>%rename(SampleID="Row.names")
#correlation

#names<- colnames(env.16S[-1])
#cor<- env.16S[-1]
#colnames(cor)<- 1:14
#corss<- cor(cor, method = "pearson") 
#corss<- cor(env.16S[-1], method = "pearson") 

#library(corrplot)
#png(filename = "cor_veg.png", width = 750, height = 730, res = 100)
#corrplot(corss,  method = 'circle', type = 'lower', insig='blank',
 #        addCoef.col ='black', number.cex = 0.5, order = 'AOE', diag=FALSE)
#dev.off()
#corss[upper.tri(corss)] <- NA
#cor2<-corss %>%as.data.frame() %>%mutate(SampleID.x = row.names(.)) %>%
#  gather(key = "SampleID.y", value = "VarDist", -SampleID.x) %>%
 # filter(!is.na(VarDist))

#write.csv(as.data.frame(cor2), "Data/cor_veg_df5.csv")
#write_tsv(env.16S, "Data/vegetation_data_notransformed.tsv")

#qiime2
mm1=env.16S%>%  inner_join(spp.16S1 %>%rownames_to_column(var = "SampleID") )
dim(mm1)
head(mm1)
env.16S1=mm1[1:15] %>% column_to_rownames(var = "SampleID")
colnames(env.16S1)
spp.16S1=mm1[c(1,16:3153)]%>% column_to_rownames(var = "SampleID")
colnames(spp.16S1)


#single
mm2=env.16S%>%  inner_join(spp.16S2 %>%rownames_to_column(var = "SampleID") )
dim(mm2)
head(mm2)
env.16S2=mm2[1:15] %>% column_to_rownames(var = "SampleID")
colnames(env.16S2)
spp.16S2=mm2[c(1,16:225)]%>% column_to_rownames(var = "SampleID")

#paired
mm3=env.16S%>%  inner_join(spp.16S3 %>%rownames_to_column(var = "SampleID") )
dim(mm3)
head(mm3)
env.16S3=mm3[1:15] %>% column_to_rownames(var = "SampleID")
dim(env.16S3)
colnames(env.16S3)
spp.16S3=mm3[c(1,16:249)]%>% column_to_rownames(var = "SampleID")
colnames(spp.16S3)

#kraken
mm4=env.16S%>%  inner_join(spp.16S4 %>%rownames_to_column(var = "SampleID") )
dim(mm4)
head(mm4)
env.16S4=mm2[1:15] %>% column_to_rownames(var = "SampleID")
colnames(env.16S4)
spp.16S2=mm4[c(1,16:87)]%>% column_to_rownames(var = "SampleID")

#transoforming enviromental and species data
library(vegan)
set.seed(125)
spp.16S_hell1=decostand(spp.16S1, "hell") # Hellinger transformation
spp.16S_hell2=decostand(spp.16S2, "hell")
spp.16S_hell3=decostand(spp.16S3, "hell")
spp.16S_hell4=decostand(spp.16S4, "hell")

spp.16S_hell1=decostand(otu.single(spp.16S1), "clr", pseudocount=0.5) # Hellinger transformation
spp.16S_hell2=decostand(otu.single(spp.16S2),  "clr", pseudocount=0.5)
spp.16S_hell3=decostand(otu.single(spp.16S3),  "clr", pseudocount=0.5)
spp.16S_hell4=decostand(otu.single(spp.16S4),  "clr", pseudocount=0.5)
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
set.seed(128)
cap.env1=capscale(spp.16S_hell1~.,env.16S_st1,distance = "euclidean")
mod0.env=capscale(spp.16S_hell1~1,env.16S_st1,distance = "euclidean")
step.env=ordistep(mod0.env,scope = formula(cap.env1))
cap1<-step.env$anova %>% mutate(method="Single Micop") # These are the selected structuring variables that are included in the subsequent analyses 

cap.env2=capscale(spp.16S_hell2~.,env.16S_st2,distance = "euclidean")
mod0.env=capscale(spp.16S_hell2~1,env.16S_st2,distance = "euclidean")
step.env=ordistep(mod0.env,scope = formula(cap.env2))
cap2<-step.env$anova %>% mutate(method="Paired Micop") 

cap.env3=capscale(spp.16S_hell3~.,env.16S_st3,distance = "euclidean")
mod0.env=capscale(spp.16S_hell3~1,env.16S_st3,distance = "euclidean")
step.env=ordistep(mod0.env,scope = formula(cap.env3))
cap3<-step.env$anova %>% mutate(method="Kraken2") 


cap.env4=capscale(spp.16S_hell4~.,env.16S_st4,distance = "euclidean")
mod0.env=capscale(spp.16S_hell4~1,env.16S_st4,distance = "euclidean")
step.env=ordistep(mod0.env,scope = formula(cap.env4))
cap4<-step.env$anova %>% mutate(method="QIIME2")#ninguno 
set.seed(129)
envs<- rbind(cap1,cap2, cap3, cap4)

#capscale
#plot
pdf(file="cap_envs_veg.pdf", width =8.5, height =6)
library(viridis)
pal<- viridis(6, option = "D") 
#plots
par(mfrow=c(2,2),mar=c(2, 2, 0.5, 0.5))
color=rgb(0,0,0,alpha=0.5) 
ordiplot(cap.env4,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env4, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env1,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env4, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)


ordiplot(cap.env1,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env1, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env2,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env1, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(cap.env2,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env2, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env3,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env2, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(cap.env3,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env3, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env4,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env3, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)



dev.off()



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
write.csv(envs, "/home/yendi/Desktop/envs_veg_FINAL.csv")

#plot cca
pdf(file="cca_envs_veg.pdf", width =8.5, height =6)
library(viridis)
pal<- viridis(6, option = "D") 
#plots
par(mfrow=c(2,2),mar=c(2, 2, 0.5, 0.5))
color=rgb(0,0,0,alpha=0.5) 
ordiplot(vares_cca4,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca4, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env1,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca4, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)


ordiplot(vares_cca1,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca1, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env2,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca1, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(vares_cca2,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca2, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env3,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca2, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(vares_cca3,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca3, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env4,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca3, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)



dev.off()


#############################################################
#QUEDE POR AQUI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###################################################################

#env.16S=env.16S
#vars_chose<-  c("prop_Alnus", "prop_Broadleaf",               
#"total_mean_DAP" , "mean_Height_Abies", "mean_Height_Pinus",
#"prop_Arbutus"  , "mean_distanciacopas_Broadleaf", "mean_distanciacopas_Pinus",
#"mean_CopaNS_Salix" ,"mean_basal_area_Abies" ,"mean_basal_area_Salix", "mean_Height_Prunus")

#vars_chose<- c( "mean_CopaNS_Abies","mean_Height_Pinus","mean_basal_area_Salix",
 #              "mean_DAP_Pinus","prop_Arbutus",  
  #             "mean_Height_Abies", "mean_CopaNS_Broadleaf",    
   #            "mean_basal_area_Abies", "mean_CopaNS_Salix"   )

vars_chose<- c("total_mean_coverage",
               "prop_Abies",
               "total_mean_Height",
               "prop_Pinus",
               "mean_coverage_Pinus",
               "mean_coverage_Broadleaf",
               "prop_Arbutus",
               "prop_Salix",
               "prop_latif",
               "prop_Alnus",
               "prop_Juniperus",
               "mean_coverage_Arbutus")

vars_chose<- c("total_mean_Height", "prop_Abies" ,       "prop_Alnus" ,
  "prop_Arbutus"     ,"prop_Pinus"     )

#c<- read.csv("/home/yendi/Desktop/envs_t.csv", row.names = 1)

df<- env.16S %>% dplyr::select(SampleID,vars_chose)
#df<- env.16S %>% dplyr::select(SampleID,rownames(c))

metadata_secas<- read_excel("Data/Metadatos.xlsx", sheet = "secas-marzo")
fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
fq_secas2<- read.csv("Data/fisicoq-la.csv")
meta_fq_secas<- metadata_secas %>% full_join(fq_secas) #%>%  select(Sites:id_new,pH, MO, N, P) %>% mutate(Season="Dry")
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")
env.16SA=meta_fq_secas_all
dfA<- env.16SA%>% dplyr::select(SampleID, pH:Mn, moisture, WHC:CONDUC, ARCILLA:ARENA)
dfsA=data.frame(dfA[1],scale(dfA[,2:17], center = T, scale = T)) %>% dplyr::select(
  SampleID,K,Mg,moisture, LIMO)
dfss<- merge(env.16S, dfsA, by = "SampleID") %>% dplyr::select(
  SampleID,vars_chose,K,Mg,moisture, LIMO)


#df1<- merge(env.16S1, dfs, by = 0)
#df1<- merge(env.16S1, dfsA, by = 0) %>% dplyr::select(Row.names,
 # prop_Alnus, prop_Broadleaf, total_mean_DAP, P:LIMO) %>% column_to_rownames(var = "Row.names") %>% 
  #dist(., method = "euclidean")
#df1<- merge(env.16S1, dfs, by = 0) %>% dplyr::select(
 # Row.names,mean_CopaNS_Abies,mean_Height_Pinus,
  #mean_basal_area_Salix,mean_DAP_Pinus, P:LIMO) %>% column_to_rownames(var = "Row.names") %>% 
  #dist(., method = "euclidean")

#df1<-env.16S1 %>% dplyr::select(mean_CopaNS_Abies,mean_Height_Pinus,
 #                              mean_basal_area_Salix,mean_distanciacopas_Pinus)  %>% 
  #dist(., method = "euclidean")
#df1<-env.16S1 %>% dplyr::select(rownames(c1))  %>% 
 # dist(., method = "euclidean")
#df2<-env.16S2 %>% dplyr::select(prop_Arbutus,prop_conif, mean_Height_Juniperus,
 #                              mean_Height_Juniperus,mean_coverage_Abies,
#mean_CopaNS_Abies ,       mean_Height_Broadleaf     )  %>% 
#dist(., method = "euclidean")

#df2<- merge(env.16S2, dfs, by=0) %>% dplyr::select(Row.names,prop_Arbutus,
 #   mean_Height_Pinus, mean_distanciacopas_Broadleaf,
  #                               mean_Height_Abies, P:LIMO)%>% column_to_rownames(var = "Row.names") %>% 
  #dist(., method = "euclidea￼n")
#df2<- merge(env.16S2, dfs, by=0) %>% dplyr::select(rownames(c2)) %>% 
 # dist(., method = "euclidean")

#df3<- env.16S3 %>% dplyr::select(mean_CopaNS_Abies) %>% 
# dist(., method = "euclidean")
#df4<- env.16S3 %>% dplyr::select(mean_basal_area_Abies,
#prop_Salix, mean_Height_Prunus  ) %>% 
# dist(., method = "euclidean")
#dfs3<- env.16S3 %>% dplyr::select("prop_Juniperus"  ,  "mean_DAP_Pinus" ,
 #                                 "mean_Height_Juniperus"  , "mean_CopaNS_Pinus" ,
  #                                "mean_CopaNS_Salix" , "mean_distanciacopas_Pinus",
   #                               "mean_diamcopa_Pinus", "mean_basal_area_Pinus"  )  %>% 
  #dist(., method = "euclidean")

dfs=df %>% 
  column_to_rownames(var = "SampleID")
dfs_dist<- dist(dfs, method = "euclidean")

dfsA=dfss %>% 
 column_to_rownames(var = "SampleID")
dfs_distA<- dist(dfsA, method = "euclidean")



nut.dist<- dfs_distA %>% as.matrix()
nut.dist.tidy <- nut.dist %>% 
  melt(as.matrix(distance), varnames = c(
    "SampleID.x", "SampleID.y"), value.name = "EucDist") %>% drop_na() %>% filter(!EucDist==0) %>% 
  filter(!is.na(EucDist)) %>% 
  filter(SampleID.x != SampleID.y) 

#nut.dist[upper.tri(nut.dist)] <- NA 

map<- read.csv("Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)

otu <- list(spp.16S4t, spp.16S1t, spp.16S2t, spp.16S3t)
otu_match<- lapply(otu, otu.match) # matching to map
otu_single<- lapply(otu_match, otu.single) #remove singletons
otu_norm<- lapply(otu_single, otu.norm)#Normalize
bc.dist2<- lapply(otu_norm, beta_div_dist)#Calculate Bray-Curtis dissimilarities
bc.dist2<- lapply(otu_single, beta_div_dist_hill, q=1)

vegan::mantel(bc.dist2[[1]], dfs_distA)
vegan::mantel(bc.dist2[[2]], dfs_distA)
vegan::mantel(bc.dist2[[3]], dfs_distA)
vegan::mantel(bc.dist2[[4]], dfs_distA)

mantel(bc.dist2[[1]], dfs_dist)
mantel(bc.dist2[[2]], dfs_dist)
mantel(bc.dist2[[3]], dfs_dist)
mantel(bc.dist2[[4]], dfs_dist)

#plot
bc.dist.tidy.filt2<- lapply(bc.dist2, bc.dist.tidy.filter)
bc.dist.tidy.filt2<- lapply(bc.dist2, bc.dist.tidy.filter.hill)


mantel_qiime2_e<-vegan::mantel(bc.dist2[[1]], dfs_distA)
mantel_single_e<-vegan::mantel(bc.dist2[[2]], dfs_distA)
mantel_paired_e<-vegan::mantel(bc.dist2[[3]], dfs_distA)
mantel_kraken_e<-vegan::mantel(bc.dist2[[4]], dfs_distA)

library(generics)
max.sim<-1
cor_test_e<- lapply(bc.dist.tidy.filt2, cor.e)
lm_test_e<- lapply(bc.dist.tidy.filt2, lm.e)

stats_qiime2_e <- data.frame(label = paste("Pearson: r = ", signif(cor_test_e[[1]]$estimate,3), 
                                           ", p-value = ", signif(cor_test_e[[1]]$p.value, 3),
                                           "\nRegression: lope = ", signif(lm_test_e[[1]]$estimate, 3),
                                           "\nMantel: r = ", signif(mantel_qiime2_e$statistic,3), 
                                           ", p-value = ", signif(mantel_qiime2_e$signif, 3)))
stats_single_e <- data.frame(label = paste("Pearson: r = ", signif(cor_test_e[[2]]$estimate,3), 
                                           ", p-value = ", signif(cor_test_e[[2]]$p.value, 3),
                                           "\nRegression: slope = ", signif(lm_test_e[[2]]$estimate, 3),
                                           "\nMantel: r = ", signif(mantel_single_e$statistic,3), 
                                           ", p-value = ", signif(mantel_single_e$signif, 3)))
stats_paired_e <- data.frame(label = paste("Pearson: r = ", signif(cor_test_e[[3]]$estimate,3), 
                                           ", p-value = ", signif(cor_test_e[[3]]$p.value, 3),
                                           "\nRegression: slope = ", signif(lm_test_e[[3]]$estimate, 3),
                                           "\nMantel: r = ", signif(mantel_paired_e$statistic,3), 
                                           ", p-value = ", signif(mantel_paired_e$signif, 3)))
stats_kraken_e <- data.frame(label = paste("Pearson: r = ", signif(cor_test_e[[4]]$estimate,3), 
                                           ", p-value = ", signif(cor_test_e[[4]]$p.value, 3),
                                           "\nRegression: slope = ", signif(lm_test_e[[4]]$estimate, 3),
                                           "\nMantel: r = ", signif(mantel_kraken_e$statistic,3), 
                                           ", p-value = ", signif(mantel_kraken_e$signif, 3)))



plot_distance_env<- lapply(bc.dist.tidy.filt2, distance.plot1)


fourth<- plot_grid(
  plot_distance_env[[1]]+  ggtitle(stats_qiime2_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title= element_text(size = 10,hjust = 0.1,vjust = -100), axis.title.x = element_blank()),
  plot_distance_env[[2]]+ ggtitle(stats_single_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),
  plot_distance_env[[3]]+ ggtitle(stats_paired_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title.x = element_blank()),  
  plot_distance_env[[4]]+  ggtitle(stats_kraken_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),  
  ncol=1)

fifth<- plot_grid(
  plot_distance_env[[1]]+  ggtitle(stats_qiime2_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title= element_text(size = 10,hjust = 0.1,vjust = -100), axis.title.x = element_blank()),
  plot_distance_env[[2]]+ ggtitle(stats_single_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),
  plot_distance_env[[3]]+ ggtitle(stats_paired_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title.x = element_blank()),  
  plot_distance_env[[4]]+  ggtitle(stats_kraken_e$label)+theme(plot.title= element_text(size = 7)),#+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),  
  ncol=1)


veg<-plot_grid(fourth, fifth, ncol = 2, labels = c("HILL", "BRAY"), label_x = -.08)

ggsave("veg21.png",width = 5.7, height = 7.5, dpi = 300, plot = veg, device = "png")


#probando genero  y especies
  data_spet<- spe_data%>%t() %>% rel_ab() %>% log_norm() %>% t()
  data_gent<- gen_data%>%t() %>% rel_ab() %>% log_norm() %>% t()

  data_spet_dist<- data_spet%>% vegdist(., method = "bray")
  data_gent_dist<- data_gent %>% vegdist(., method = "bray")
  
  
  
  mantel(bc.dist2[[1]], data_spet_dist)
  cor.test(bc.dist2[[1]], as.matrix(data_spet_dist))
  
  
  mantel(bc.dist2[[2]], data_spet_dist)
  cor.test(bc.dist2[[2]], as.matrix(data_spet_dist))
  
  mantel(bc.dist2[[3]], data_spet_dist)
  cor.test(bc.dist2[[3]], as.matrix(data_spet_dist))

  mantel(bc.dist2[[4]], data_spet_dist)
  cor.test(bc.dist2[[4]], as.matrix(data_spet_dist))
  
  
  mantel(bc.dist2[[1]], data_gent_dist)
  cor.test(bc.dist2[[1]], as.matrix(data_gent_dist))
  
  
  mantel(bc.dist2[[2]], data_gent_dist)
  cor.test(bc.dist2[[2]], as.matrix(data_gent_dist))
  
  mantel(bc.dist2[[3]], data_gent_dist)
  cor.test(bc.dist2[[3]], as.matrix(data_gent_dist))

  mantel(bc.dist2[[4]], data_gent_dist)
  cor.test(bc.dist2[[4]], as.matrix(data_gent_dist))
  
  
  
