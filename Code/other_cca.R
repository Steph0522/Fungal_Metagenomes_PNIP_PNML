#plot
library(ggvegan)
map<- read.csv("Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
map1<- map[match(rownames(spp.16S_hell1), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()
varscho1<- c("62b57b7052f748d8888528f055a56166eab2eaf2","42c34c5cae8e20681d23626859950c1658f52183",
             "4dcb134727786b6b7cb9548cee1ac3171a946846")
varscho2<- c("Eukaryota__Ascomycota__Sordariomycetes__Glomerellales__Glomerellaceae__Colletotrichum__Colletotrichum fioriniae",
             "Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Aspergillus__Aspergillus bombycis",
             "Eukaryota__Basidiomycota__Agaricomycetes__Agaricales__Tricholomataceae__Laccaria__Laccaria bicolor")

varscho3<-c("Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Aspergillus__Aspergillus glaucus",
            "Eukaryota__Ascomycota__Saccharomycetes__Saccharomycetales__Metschnikowiaceae__Clavispora/Candida clade__[Candida] auris",
            "Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Penicillium__Penicillium arizonense")
varscho4<- c("101028","80884","5062")


envschose<-c("P", "K", "Ca", "Mg", "P", "moisture", "WHC", "Silt")
varschose<- c("total_mean_coverage",
               "prop_Abies",
               "total_mean_Height",
               "prop_Pinus",
               "mean_coverage_Pinus",
               "mean_coverage_Quercus",
               "prop_Arbutus",
               "prop_Salix",
               "prop_Broadleaf",
               "prop_Alnus",
               "prop_Juniperus",
               "mean_coverage_Arbutus")

#ccas

vars_choosing1<- fortify(vares_cca1) %>% filter(Score=="species")%>%
  inner_join(taxonomy_qiime2, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho1) %>% mutate(
      CCA1=CCA1*1, CCA2=CCA2*1) %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing2<- fortify(vares_cca2) %>% filter(Score=="species")%>%
  inner_join(taxonomy_single_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho2) %>% mutate(
      CCA1=CCA1*2, CCA2=CCA2*2)  %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing3<- fortify(vares_cca3) %>% filter(Score=="species")%>%
  inner_join(taxonomy_paired_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho3)%>% mutate(
      CCA1=CCA1*10, CCA2=CCA2*10) %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing4<- fortify(vares_cca4) %>% filter(Score=="species")%>%
  inner_join(taxonomy_fungi, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho4) %>% mutate(
      CCA1=CCA1*52, CCA2=CCA2*52) %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

pdf(file="cca_envs_veg_signif.pdf", width =8.5, height =6)
library(viridis)
#plots
par(mfrow=c(2,2),mar=c(4,4, 0.5, 0.5))
#color=rgb(0,0,0,alpha=0.5) 
ordiplot(vares_cca1,type="n")
text(vares_cca1, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(vares_cca1, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing1$CCA1,y = vars_choosing1$CCA2, 
     cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)




ordiplot(vares_cca2,type="n")
text(vares_cca2, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(vares_cca2, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing2$CCA1,y = vars_choosing2$CCA2, 
     cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)



ordiplot(vares_cca3,type="n")
text(vares_cca3, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(vares_cca3, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing4$CCA1,y = vars_choosing4$CCA2, 
     cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



ordiplot(vares_cca4,type="n")
text(vares_cca4, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(vares_cca4, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing4$CCA1,y = vars_choosing4$CCA2, 
     cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()




#caps


#map1<- map[match(rownames(spp.16S_hell1), map$SampleID),]

vars_choosing1<- fortify(cap.env1) %>% filter(Score=="species")%>%
  inner_join(taxonomy_qiime2, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho1) %>% mutate(
      CAP1=CAP1*120, CAP2=CAP2*120) %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing2<- fortify(cap.env2) %>% filter(Score=="species")%>%
  inner_join(taxonomy_single_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho2) %>% mutate(
     CAP1=CAP1*15, CAP2=CAP2*15)  %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing3<- fortify(cap.env3) %>% filter(Score=="species")%>%
  inner_join(taxonomy_paired_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho3)%>% mutate(
      CAP1=CAP1*3, CCA2=CAP2*3) %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing4<- fortify(cap.env4) %>% filter(Score=="species")%>%
  inner_join(taxonomy_fungi, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho4) %>% mutate(
      CAP1=CAP1*2, CAP2=CAP2*2) %>% 
  mutate(tax=case_when(
    tax=="pseudograminearum" ~ "Fusarium",
    tax=="higginsianum" ~ "Colletotrichum",
    tax=="oryzae" ~ "Aspergillus",
    tax=="Colletotrichum fioriniae" ~ "Colletitruchum",
    tax=="Aspergillus bombycis" ~ "Aspergillus",
    tax=="Aspergillus glaucus" ~ "Aspergillus",
    tax=="Aspergillus oryzae" ~ "Aspergillus",
    tax=="Laccaria bicolor" ~ "Laccaria",
    tax=="[Candida] auris" ~ "Clavispora",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))


#plot
pdf(file="cap_envs_veg_sig_clr.pdf", width =8.5, height =6)
#jpeg(file="cap_envs_signif.jpeg")

library(viridis)
pal<- viridis(12, option = "H") 
#plots
par(mfrow=c(2,2),mar=c(4, 4, 0.5, 0.5))
color=rgb(0,0,0,alpha=1) 

ordiplot(cap.env1,type="n",  ylim = c(-2.5,2.5))
text(cap.env1, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(cap.env1, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, cex=0.6)
text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
     cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)


ordiplot(cap.env2,type="n")
text(cap.env2, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(cap.env2, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
     cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(cap.env3,type="n")
text(cap.env3, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(cap.env3, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
     cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(cap.env4,type="n")
text(cap.env4, dis="cn", scaling="sites", cex=0.8, col="red")
ordiellipse(cap.env4, groups =map1$Site, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.6)
text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
     cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)

dev.off()

