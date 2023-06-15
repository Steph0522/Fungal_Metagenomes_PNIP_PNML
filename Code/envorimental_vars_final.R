##################
# Data analyses #
#################


# Assessing environmental conditions and selecting structuring variables
library(reshape2)
library(tidyverse)
library(picante)
library(Rmisc)
library(readxl)
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
spp.16S2=read.delim("Data/table_micop_single.txt") 
spp.16S3=read.delim("Data/table_micop_paired.txt") 

spp.16S4=read.delim("Data/table_fungi_again.txt", 
                   skip = 1, row.names = 1, check.names = F) %>% dplyr::select_if(
                     is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                       var = "id_sequence") %>% separate(
                         ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                         sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                           metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names,  -id_fisicoq) %>% column_to_rownames(
                             var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)
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


spp.16S1=spp.16S1[, match(colnames(spp.16S2), colnames(spp.16S1))]



env.16S=meta_fq_secas_all

spp.16S1=data.frame(t(spp.16S1), check.names = F)
spp.16S2=data.frame(t(spp.16S2), check.names = F)
spp.16S3=data.frame(t(spp.16S3), check.names = F)
spp.16S4=data.frame(t(spp.16S4), check.names = F)

#View(mm)


#select variables that did not show collinearity

#kraken
mm4=cbind(spp.16S4,env.16S)
#mm3=na.omit(mm3)
dim(mm4)
head(mm4)
env.16S4=mm4[,c(83:93,96,99,101)]
dim(env.16S4)
colnames(env.16S4)
spp.16S4=mm4[,1:72]
colnames(spp.16S4)

#single
mm2=cbind(spp.16S2,env.16S)
#mm1=na.omit(mm1)
dim(mm2)
head(mm2)
colnames(mm2)
#env.16S1=mm1[,c(221:231,234,235, 237:239)]
env.16S2=mm2[,c(221:231,234, 237,239)]
colnames(env.16S2)
spp.16S2=mm2[,1:210]
colnames(spp.16S2)

#paired
mm3=cbind(spp.16S3,env.16S)
#mm2=na.omit(mm2)
dim(mm3)
colnames(mm3)
head(mm3)
env.16S3=mm3[,c(245:255,258, 261,263)]
spp.16S3=mm3[,1:234]
colnames(env.16S3)
colnames(spp.16S3)

#qiime2
mm1=cbind(spp.16S1,env.16S)
#mm4=na.omit(mm4)
tail(colnames(mm1), n = 100)
dim(mm1)
head(mm1)
env.16S1=mm1[,c(3149:3159, 3162,  3165,3167)]
colnames(env.16S1)
dim(env.16S1)
spp.16S1=mm1[,1:3138]
colnames(spp.16S1)
#dcas
#dca=decorana(spp.16S1)
#plot(dca) # RDA is okay to use

#transoforming enviromental and species data

spp.16S_hell1=decostand(spp.16S1, "hell") # Hellinger transformation
spp.16S_hell2=decostand(spp.16S2, "hell")
spp.16S_hell3=decostand(spp.16S3, "hell")
spp.16S_hell4=decostand(spp.16S4, "hell")

spp.16S_hell1=decostand(otu.single(spp.16S1), "clr", pseudocount=0.5) # Hellinger transformation
spp.16S_hell2=decostand(otu.single(spp.16S2),  "clr", pseudocount=0.5)
spp.16S_hell3=decostand(otu.single(spp.16S3),  "clr", pseudocount=0.5)
spp.16S_hell4=decostand(otu.single(spp.16S4),  "clr", pseudocount=0.5)

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
cap.env1=capscale(spp.16S_hell4.,env.16S_st4,distance = "euclidean")
mod0.env=capscale(spp.16S_hell4~1,env.16S_st4,distance = "euclidean")
step.env1=ordistep(mod0.env,scope = formula(cap.env1))
cap1<-step.env1$anova %>% mutate(method="QIIME2") # These are the selected structuring variables that are included in the subsequent analyses 


cap.env2=capscale(spp.16S_hell2~.,env.16S_st2,distance = "bray")
mod0.env=capscale(spp.16S_hell2~1,env.16S_st2,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env2))
cap2<-step.env$anova %>% mutate(method="Single Micop") 

cap.env3=capscale(spp.16S_hell3~.,env.16S_st3,distance = "bray")
mod0.env=capscale(spp.16S_hell3~1,env.16S_st3,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env3))
cap3<-step.env$anova %>% mutate(method="Paired Micop") 


cap.env4=capscale(spp.16S_hell4~.,env.16S_st4,distance = "bray")
mod0.env=capscale(spp.16S_hell4~1,env.16S_st4,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env4))
cap4<-step.env$anova %>% mutate(method="Kraken2") 

#plot
#plot
pdf(file="cap_envs.pdf", width =8.5, height =6)
library(viridis)
pal<- viridis(6, option = "D") 
#plots
par(mfrow=c(2,2),mar=c(2, 2, 0.5, 0.5))
color=rgb(0,0,0,alpha=0.5) 
ordiplot(cap.env1,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env1, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env1,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env1, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)


ordiplot(cap.env2,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env2, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env2,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env2, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(cap.env3,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env3, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env3,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env3, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(cap.env4,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env4, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env4,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env4, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)



dev.off()



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
vares_cca1.1 <- cca(spp.16S_hell1 ~ P+K+Ca+Mg+moisture+WHC+ARCILLA ,
                  data=env.16S_st1)
env1<-envfit(vares_cca1 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
       data=env.16S_st1)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("QIIME2"=".")

vares_cca2 <- cca(spp.16S_hell2 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                  data=env.16S_st2)
vares_cca2.1 <- cca(spp.16S_hell2 ~ P+K+Ca+Mg+moisture+WHC+ARCILLA ,
                    data=env.16S_st2)
env2<-envfit(vares_cca2 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
             data=env.16S_st2)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Single Micop"=".")

vares_cca3 <- cca(spp.16S_hell3 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                  data=env.16S_st3)
vares_cca3.1 <- cca(spp.16S_hell3 ~ P+K+Ca+Mg+moisture+WHC+ARCILLA ,
                    data=env.16S_st3)
env3<-envfit(vares_cca3 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
             data=env.16S_st3)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Paired Micop"=".")

vares_cca4 <- cca(spp.16S_hell4 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
                  data=env.16S_st4)
vares_cca4.1 <- cca(spp.16S_hell4 ~ P+K+Ca+Mg+moisture+WHC+ARCILLA ,
                    data=env.16S_st4)
env4<-envfit(vares_cca4 ~ pH+MO+N+P+K+Ca+Mg+Fe+Cu+Mn+moisture+WHC+ARCILLA+ARENA ,
             data=env.16S_st4)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Kraken2"=".")

envs<- cbind(env4, env1, env2, env3)
rownames(envs)

pdf(file="cca_envs.pdf", width =8.5, height =6)
library(viridis)
pal<- viridis(6, option = "D") 
#plots
par(mfrow=c(2,2),mar=c(2, 2, 0.5, 0.5))
color=rgb(0,0,0,alpha=0.5) 
ordiplot(vares_cca1,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca1, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env1,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca1, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)


ordiplot(vares_cca2,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca2, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env2,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca2, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(vares_cca3,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca3, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env3,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca3, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(vares_cca4,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca4, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env4,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca4, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)



dev.off()

