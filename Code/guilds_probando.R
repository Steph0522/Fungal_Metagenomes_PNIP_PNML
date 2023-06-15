library(readxl)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggh4x)
metadata<-read_excel("../Data/Metadatos.xlsx") %>% mutate_if(
  is.numeric, as.factor)
meta<-metadata
guild_colors<- read_tsv("../Data/colors_guild")
trophic_colors<- read_tsv("../Data/colors_trophic")

names_col<- read.delim(
  "../Data/tabla_qiime2_blast_OTUS.guilds.txt", 
  check.names = F, row.names = 1)%>% t() %>% 
  as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% separate(
      SampleID, c(
        "id_metagenome", "R", "unmap", "Paired"), 
      sep = "_")%>% inner_join(meta) %>% dplyr::select(
        SampleID) %>% t() %>% as.vector()

guild_qiime2<- read.delim(
  "../Data/tabla_qiime2_blast_OTUS.guilds.txt", 
  check.names = F, row.names = 1) %>% filter(
    !Guild=="-",! `Confidence Ranking`=="Possible" ) %>% mutate(Method="QIIME2")

colnames(guild_qiime2)[1:36]<-names_col

guild_qiime2<- guild_qiime2 %>% rownames_to_column(var = "otu")

guild_qiime2_colap<- guild_qiime2 %>% group_by(Guild) %>% count()
guild_qiime2_colapsed<- guild_qiime2 %>% group_by(Guild) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Guild")
trophic_qiime2_colapsed<- guild_qiime2 %>% group_by(`Trophic Mode`) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Trophic Mode")


guild_single<- read.delim("../Data/table_micop_single.guilds.txt", 
                          check.names = F)%>% filter(
                            !Guild=="-",! `Confidence Ranking`=="Possible" ) %>% mutate(Method="MICOP SINGLE")

guild_single_colap<- guild_single %>% group_by(Guild) %>% count()
guild_single_colapsed<- guild_single %>% group_by(Guild) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Guild")
trophic_single_colapsed<- guild_single %>% group_by(`Trophic Mode`) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Trophic Mode")


guild_paired<- read.delim("../Data/table_micop_paired.guilds.txt", 
                          check.names = F)%>% filter(
                            !Guild=="-",! `Confidence Ranking`=="Possible" ) %>% mutate(Method="MICOP PAIRED")
guild_paired_colap<- guild_paired %>% group_by(Guild) %>% count()
guild_paired_colapsed<- guild_paired %>% group_by(Guild) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Guild")
trophic_paired_colapsed<- guild_paired %>% group_by(`Trophic Mode`) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Trophic Mode")



guild_fungi<- read.delim("../Data/table_fungi_again.guilds.txt", 
                         check.names = F, row.names = 1) %>% filter(
                           !Guild=="-",! `Confidence Ranking`=="Possible" ) %>% mutate(otu=paste0("otu", rownames(.))) %>% column_to_rownames(var = "otu")%>% select_if(is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                             var = "id_sequence") %>% separate(
                               ., "id_sequence", c(
                                 "kraken", "fungi", "id_metagenome", "report", "bracken"), 
                               sep = "_") %>% dplyr::select(
                                 -kraken, -fungi, -report, -bracken) %>% full_join(metadata) %>% dplyr::select(-id_sequence:-Names, -id_metagenome) %>% column_to_rownames(var = "SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "otu") %>% inner_join(guilds_fungi)%>% mutate(Method="KRAKEN2")

guild_fungi_colap<- guild_fungi %>% group_by(Guild) %>% count()
guild_fungi_colapsed<- guild_fungi %>% group_by(Guild)  %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Guild")
trophic_fungi_colapsed<- guild_fungi %>% group_by(`Trophic Mode`) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(var = "Trophic Mode")


relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}

guild_all<- rbind(guild_qiime2, guild_paired, guild_single, guild_fungi)


guilds_all<- guild_qiime2 %>% full_join(guild_paired) %>% full_join(guild_single) %>% full_join(guild_fungi)
#guilds_all<- guilds_all %>% dplyr::select_if(is.character) %>% dplyr::select(
 # otu, taxonomy, Taxon, Guild, `Trophic Mode`, `Confidence Ranking`) %>% filter(!`Confidence Ranking`=="Possible")
#write_tsv(guilds_all, "../Data/guilds_all.tsv")

guild_col<- guild_qiime2_colap %>% dplyr::rename("QIIME2"=n) %>% 
  full_join(guild_single_colap)%>% dplyr::rename("SINGLE MICOP"=n)%>% 
  full_join(guild_paired_colap)%>% dplyr::rename("PAIRED MICOP"=n)%>% 
  full_join(guild_fungi_colap)%>% dplyr::rename("KRAKEN2"=n)%>% 
  replace(is.na(.), 0)%>% mutate(
    total = rowSums(across(where(is.numeric)))) %>% arrange(-total) %>% 
  dplyr::select(-total) %>% filter(!Guild=="-" & !Guild=="NULL")

#write_csv(guild_col, "../Data/guild_colap.csv")


guild1<- trophic_qiime2_colapsed %>% relabunda() %>% log_norm()%>% t() %>% as.data.frame() %>%rownames_to_column(var = "SampleID") 
guild2<- trophic_single_colapsed %>% relabunda() %>% log_norm()%>% t() %>% as.data.frame()%>%rownames_to_column(var = "SampleID") 
guild3<- trophic_paired_colapsed %>% relabunda() %>% log_norm()%>% t() %>% as.data.frame()%>%rownames_to_column(var = "SampleID") 
guild4<- trophic_fungi_colapsed %>% relabunda() %>% log_norm()%>% t() %>% as.data.frame()%>%rownames_to_column(var = "SampleID") 

env_vars<- read.csv("../Data/cors_envs.csv", row.names = 1)
veg_vars<- read.delim("../Data/vegeta.tsv", row.names = 1)

env_guild1<- env_vars %>% rownames_to_column(var = "SampleID") %>% full_join(guild1)
cor1<- cor(env_guild1[-1], method = "pearson")
cor1m<- reshape2::melt(data = cor1, value.name = "cor")
cor1m<- cor1m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)


env_guild2<- env_vars %>% rownames_to_column(var = "SampleID") %>% 
  full_join(guild2)
cor2<- cor(env_guild2[-1], method = "pearson")
cor2m<- reshape2::melt(data = cor2, value.name = "cor")
cor2m<- cor2m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)

env_guild3<- env_vars %>% rownames_to_column(var = "SampleID") %>% 
  full_join(guild3)
cor3<- cor(env_guild3[-1], method = "pearson")
cor3m<- reshape2::melt(data = cor3, value.name = "cor")
cor3m<- cor3m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)


env_guild4<- env_vars %>% rownames_to_column(var = "SampleID") %>% 
  full_join(guild4)
cor4<- cor(env_guild4[-1], method = "pearson")
cor4m<- reshape2::melt(data = cor4, value.name = "cor")
cor4m<- cor4m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)

write_tsv(cor4m, "../Data/tcor4.tsv")


veg_guild1<- veg_vars %>% rownames_to_column(var = "SampleID") %>% full_join(guild1)
cor1<- cor(veg_guild1[-1], method = "pearson")
cor1m<- reshape2::melt(data = cor1, value.name = "cor")
cor1m<- cor1m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)


veg_guild2<- veg_vars %>% rownames_to_column(var = "SampleID") %>% 
  full_join(guild2)
cor2<- cor(veg_guild2[-1], method = "pearson")
cor2m<- reshape2::melt(data = cor2, value.name = "cor")
cor2m<- cor2m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)

veg_guild3<- veg_vars %>% rownames_to_column(var = "SampleID") %>% 
  full_join(guild3)
cor3<- cor(veg_guild3[-1], method = "pearson")
cor3m<- reshape2::melt(data = cor3, value.name = "cor")
cor3m<- cor3m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)


veg_guild4<- veg_vars %>% rownames_to_column(var = "SampleID") %>% 
  full_join(guild4)
cor4<- cor(veg_guild4[-1], method = "pearson")
cor4m<- reshape2::melt(data = cor4, value.name = "cor")
cor4m<- cor4m %>% filter(!cor==1&!cor==0) %>% mutate(cors=abs(cor)) %>% filter(cors>0.6)


write_tsv(cor4m, "../Data/tcor4_v.tsv")

#selecting using cca (no distance)

#qiime2
mm1<-guild_qiime2_colapsed%>%t() %>% as.data.frame() %>%  rownames_to_column(var = "SampleID") %>% inner_join(
  env) %>% column_to_rownames(var = "SampleID")
env1=mm1[,vars]
#env1=mm1[,c(3149:3159, 3162,  3165,3167)]
spp1=mm1[,1:39]

#single
mm2<-guild_single_colapsed%>%t() %>% as.data.frame() %>%  rownames_to_column(var = "SampleID") %>% inner_join(
  env) %>% column_to_rownames(var = "SampleID")
env2=mm2[,vars]
spp2=mm2[,1:27]

#paired
mm3<-guild_paired_colapsed%>%t() %>% as.data.frame() %>%  rownames_to_column(var = "SampleID") %>% inner_join(
  env) %>% column_to_rownames(var = "SampleID")
env3=mm3[,vars]
spp3=mm3[,1:27]

#kraken
mm4<-guild_fungi_colapsed%>%t() %>% as.data.frame() %>%  rownames_to_column(var = "SampleID") %>% inner_join(
  env) %>% column_to_rownames(var = "SampleID")
env4=mm4[,vars]
spp4=mm4[,1:11]


# Hellinger transformation

spp_hell1=decostand(spp1, "hell") 
spp_hell2=decostand(spp2, "hell")
spp_hell3=decostand(spp3, "hell")
spp_hell4=decostand(spp4, "hell")

# Compositional tranformation
spp_clr1=decostand(otu.single(spp1), "clr", pseudocount=0.5) 
spp_clr2=decostand(otu.single(spp2),  "clr", pseudocount=0.5)
spp_clr3=decostand(otu.single(spp3),  "clr", pseudocount=0.5)
spp_clr4=decostand(otu.single(spp4),  "clr", pseudocount=0.5)

# Transforming environmental data
envscho<-c("MO", "Ca", "Silt", "WHC", "Fe", "K", "Mg", "pH", "Cu", "Clay", "P", "Humidity", "N")
envst1=data.frame(scale(env1, scale=T, center=F)) %>% dplyr::rename(Silt=LIMO, Humidity=moisture, Clay=ARCILLA) %>% dplyr::select(all_of(envscho))
envst2=data.frame(scale(env2, scale=T, center=F))%>% dplyr::rename(Silt=LIMO, Humidity=moisture, Clay=ARCILLA) %>% dplyr::select(all_of(envscho))
envst3=data.frame(scale(env3, scale=T, center=F))%>% dplyr::rename(Silt=LIMO, Humidity=moisture, Clay=ARCILLA) %>% dplyr::select(all_of(envscho))
envst4=data.frame(scale(env4, scale=T, center=F))%>% dplyr::rename(Silt=LIMO, Humidity=moisture, Clay=ARCILLA) %>% dplyr::select(all_of(envscho))

# Forward selection procedure
#accounting for abundance
cap.env1=capscale(spp_hell1~.,envst1,distance = "bray")
mod0.env1=capscale(spp_hell1~1,envst1,distance = "bray")
step.env1=ordistep(mod0.env1,scope = formula(cap.env1))
cap1<-step.env1$anova %>% mutate(method="QIIME2") #no significants
meta1<- env1 %>% rownames_to_column(var = "SampleID") %>%  inner_join(metadata)

cap.env2=capscale(spp_hell2~.,envst2,distance = "bray")
mod0.env2=capscale(spp_hell2~1,envst2,distance = "bray")
step.env2=ordistep(mod0.env2,scope = formula(cap.env2))
cap2<-step.env2$anova %>% mutate(method="Single Micop") 
meta2<- env2 %>% rownames_to_column(var = "SampleID") %>%  inner_join(metadata)


cap.env3=capscale(spp_hell3~.,envst3,distance = "bray")
mod0.env3=capscale(spp_hell3~1,envst3,distance = "bray")
step.env3=ordistep(mod0.env3,scope = formula(cap.env3))
cap3<-step.env3$anova %>% mutate(method="Paired Micop") 
meta3<- env3 %>% rownames_to_column(var = "SampleID") %>%  inner_join(metadata)

cap.env4=capscale(spp_hell4~.,envst4,distance = "bray")
mod0.env4=capscale(spp_hell4~1,envst4,distance = "bray")
step.env4=ordistep(mod0.env4,scope = formula(cap.env4))
cap4<-step.env4$anova %>% mutate(method="Kraken2") 
meta4<- env4 %>% rownames_to_column(var = "SampleID") %>%  inner_join(metadata)

#compositional
cap.env1_c=capscale(spp_clr1~.,envst1,distance = "euclidean")
mod0.env1_c=capscale(spp_clr1~1,envst1,distance = "euclidean")
step.env1_c=ordistep(mod0.env1_c,scope = formula(cap.env1_c))
cap1_c<-step.env1_c$anova %>% mutate(method="QIIME2")

cap.env2_c=capscale(spp_clr2~.,envst2,distance = "euclidean")
mod0.env2_c=capscale(spp_clr2~1,envst2,distance = "euclidean")
step.env2_c=ordistep(mod0.env2_c,scope = formula(cap.env2_c))
cap2_c<-step.env2_c$anova %>% mutate(method="Single Micop") 

cap.env3_c=capscale(spp_clr3~.,envst3,distance = "euclidean")
mod0.env3_c=capscale(spp_clr3~1,envst3,distance = "euclidean")
step.env3_c=ordistep(mod0.env3_c,scope = formula(cap.env3_c))
cap3_c<-step.env3_c$anova %>% mutate(method="Paired Micop") 

cap.env4_c=capscale(spp_clr4~.,envst4,distance = "euclidean")
mod0.env4_c=capscale(spp_clr4~1,envst4,distance = "euclidean")
step.env4_c=ordistep(mod0.env4_c,scope = formula(cap.env4_c))
cap4_c<-step.env4_c$anova %>% mutate(method="Kraken2")  #no significant

#selecting using cca (no distance)
vares_cca1 <- cca(spp_hell1 ~., data=envst1)
envs1<-envfit(vares_cca1 ~ ., data=envst1)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("QIIME2"=".")

vares_cca2 <- cca(spp_hell2 ~., data=envst2)
envs2<-envfit(vares_cca2 ~ ., data=envst2)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Single Micop"=".")

vares_cca3 <- cca(spp_hell3 ~.,data=envst3)
envs3<-envfit(vares_cca3 ~ ., data=envst3)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Paired Micop"=".")

vares_cca4 <- cca(spp_hell4 ~.,data=envst4)
envs4<-envfit(vares_cca4 ~ ., data=envst4)$vectors[["pvals"]]%>% as.data.frame() %>% dplyr::select("Kraken2"=".")

envs<- cbind(envs1, envs2, envs3, envs4)
rownames(envs)
View(envs)
