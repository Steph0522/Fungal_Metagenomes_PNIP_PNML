##################
# Data analyses #
#################


# Assessing environmental conditions and selecting structuring variables
library(reshape2)
library(tidyverse)
library(picante)
library(Rmisc)
metadata_secas<- read_excel("Data/Metadatos.xlsx", sheet = "secas-marzo")
fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
fq_secas2<- read.csv("Data/fisicoq-la.csv")


meta_fq_secas<- metadata_secas %>% full_join(fq_secas) #%>%  select(Sites:id_new,pH, MO, N, P) %>% mutate(Season="Dry")
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")

env.16S=meta_fq_secas_all
df<- env.16S%>% dplyr::select(Sites, pH:Mn, moisture, WHC:CONDUC, ARCILLA:ARENA) 
dfs=data.frame(scale(df[,2:17], center = T, scale = T)) #standardize environmental data

dfs$Sites<-NA
dfs$Sites<- df$Sites
melted_df=melt(dfs, id="Sites")
mean_table=summarySE(melted_df, measurevar = "value", groupvars=c("Sites","variable"), na.rm = T)

# checking collinearity
cotable<-cor.table(na.omit(dfs[,1:16]))$r
cors<- dfs %>% dplyr::select(pH:Mn, Humidity=moisture, WHC, Ec=CONDUC, Silt=ARCILLA,
                                           Lime=LIMO, Sand=ARENA) 
corss<- cor(cors, method = "pearson")
library(corrplot)
#par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
corrplot(corss, type="upper")
testRes = cor.mtest(cors, conf.level = 0.95)

#png(filename = "cor.png", width = 830, height = 630, res = 100)
corrplot(corss,  method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
#dev.off()

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
                         ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                         sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                           metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names,  -id_fisicoq) %>% column_to_rownames(
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
                           -Names) %>% column_to_rownames(
                             var="SampleID") %>% t(
                             ) %>% as.data.frame() %>% mutate_all(as.numeric)



env.16S=meta_fq_secas_all

spp.16S1=data.frame(t(spp.16S1))
spp.16S2=data.frame(t(spp.16S2))
spp.16S3=data.frame(t(spp.16S3))
spp.16S4=data.frame(t(spp.16S4))

#View(mm)


#select variables that did not show collinearity

#kraken
mm3=cbind(spp.16S3,env.16S)
#mm3=na.omit(mm3)
dim(mm3)
head(mm3)
env.16S3=mm3[,c(85:95,98,101,103)]
dim(env.16S3)
colnames(env.16S3)
spp.16S3=mm3[,1:72]


#single
mm1=cbind(spp.16S1,env.16S)
#mm1=na.omit(mm1)
dim(mm1)
head(mm1)
#env.16S1=mm1[,c(221:231,234,235, 237:239)]
env.16S1=mm1[,c(221:231,234, 237,239)]
colnames(env.16S1)
spp.16S1=mm1[,1:210]

#paired
mm2=cbind(spp.16S2,env.16S)
#mm2=na.omit(mm2)
dim(mm2)
head(mm2)
env.16S2=mm2[,c(245:255,258, 261,263)]
spp.16S=mm2[,1:234]

#qiime2
mm4=cbind(spp.16S4,env.16S)
#mm4=na.omit(mm4)
dim(mm4)
head(mm4)
env.16S4=mm4[,c(3151:3161, 3164,  3167,3169)]
dim(env.16S4)
spp.16S4=mm4[,1:3138]


#dcas
#dca=decorana(spp.16S1)
#plot(dca) # RDA is okay to use

#transoforming enviromental and species data

spp.16S_hell1=decostand(spp.16S1, "hell") # Hellinger transformation
spp.16S_hell2=decostand(spp.16S2, "hell")
spp.16S_hell3=decostand(spp.16S3, "hell")
spp.16S_hell4=decostand(spp.16S4, "hell")

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
cap1<-step.env1$anova %>% mutate(method="Single Micop") # These are the selected structuring variables that are included in the subsequent analyses 

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
#cap4<-step.env$anova %>% mutate(method="QIIME2") 

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
#cap4j<-step.env$anova %>% mutate(method="QIIME2") 

vars_caps<- rbind(cap1, cap2, cap3)
vars_capsj<- rbind(capj1, capj2, capj3)

rownames(vars_caps)#vars to select


#selecting using cca
vares_cca1 <- cca(spp.16S_hell1 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                 data=env.16S_st1)
env1<-envfit(vares_cca1 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
       data=env.16S_st1)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Single Micop"=".")

vares_cca2 <- cca(spp.16S_hell2 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                  data=env.16S_st2)
env2<-envfit(vares_cca2 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
             data=env.16S_st2)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Paired Micop"=".")

vares_cca3 <- cca(spp.16S_hell3 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                  data=env.16S_st3)
env3<-envfit(vares_cca3 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
             data=env.16S_st3)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Kraken2"=".")

vares_cca4 <- cca(spp.16S_hell4 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                  data=env.16S_st4)
env4<-envfit(vares_cca4 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
             data=env.16S_st4)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("QIIME2"=".")

envs<- cbind(env4, env1, env2, env3)
