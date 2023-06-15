#plot
taxonomy_qiime2<- data.frame(read_qza("../Data/taxonomy_blast_dfc_0.98.qza")$data, check.names = F) %>% dplyr::select(Feature.ID,Taxon)
taxonomy_single_micop<- read.delim("../Data/table_micop_single.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_paired_micop<- read.delim("../Data/table_micop_paired.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_fungi<- read.delim("../Data/table_fungi_again.txt",
                            skip = 1, row.names = 1, check.names = F) %>% select_if(is.character) %>% rownames_to_column(
                              var = "#OTU ID") %>% dplyr::select(Feature.ID = "#OTU ID", Taxon= taxonomy)

map<- read.csv("../Data/coord.csv") %>% mutate_at(
  c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
map1<- map[match(rownames(spp.16S_hell1), map$SampleID),]
map1_type<- map1 %>% dplyr::select(Site, Type) %>% group_by(Site, Type) %>% count()
varscho1<- c("62b57b7052f748d8888528f055a56166eab2eaf2","42c34c5cae8e20681d23626859950c1658f52183",
             "4dcb134727786b6b7cb9548cee1ac3171a946846", "0d79ac9692c8006f7aa5165d05dede1eaa75d901")
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

vars_choosing1<- fortify(cap.env1) %>% filter(Score=="species")%>%
  inner_join(taxonomy_qiime2, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho1) %>% mutate(
      CAP1=CAP1*20, CAP2=CAP2*20) %>% 
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
    tax=="unidentified" ~ "Sebacina",
    tax=="Penicillium arizonense" ~ "Penicillium",
    TRUE ~ as.character(tax)))

vars_choosing2<- fortify(cap.env2) %>% filter(Score=="species")%>%
  inner_join(taxonomy_single_micop, by = c("Label"="Feature.ID")) %>% mutate(
    tax= str_extract(Taxon, "[^_]+$"))%>% filter(Label %in% varscho2) %>% mutate(
      CAP1=CAP1*5, CAP2=CAP2*10)  %>% 
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
      CAP1=CAP1*5, CAP2=CAP2*5) %>% 
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
      CAP1=CAP1*5, CAP2=CAP2*5) %>% 
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

png("cap_envs.png",width=6,height=4,units="in",res=1200)
library(viridis)
pal<- viridis(12, option = "H") 
pal2<- c( "#30123BFF","#30123BFF","#30123BFF","#4454C4FF","#4454C4FF","#4454C4FF",
                     "#4490FEFF","#4490FEFF","#4490FEFF" ,"#1FC8DEFF","#1FC8DEFF","#1FC8DEFF",
                     "#29EFA2FF","#29EFA2FF","#29EFA2FF", "#7DFF56FF","#7DFF56FF","#7DFF56FF",
                     "#C1F334FF","#C1F334FF","#C1F334FF", "#F1CA3AFF", "#F1CA3AFF", "#F1CA3AFF",
                     "#FE922AFF","#FE922AFF","#FE922AFF" , "#EA4F0DFF", "#EA4F0DFF", "#EA4F0DFF",
                     "#BE2102FF","#BE2102FF","#BE2102FF", "#7A0403FF","#7A0403FF","#7A0403FF")
                     #plots
par(mfrow=c(2,2),mar=c(2, 2, 1,1))

ordiplot(cap.env1,type="n", main="QIIME2", cex.main=1)
text(cap.env1, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env1,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env1, groups =meta1$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)
#text(cap.env1, dis="species", scaling="sites", cex=0.5, col="black")


ordiplot(cap.env2,type="n", main="SINGLE MICOP", cex.main=1)
text(cap.env2, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env2,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env2, groups =meta2$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(cap.env3,type="n", main="PAIRED MICOP", cex.main=1)
text(cap.env3, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env3,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env3, groups =meta3$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)

#text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(cap.env4,type="n", main="KRAKEN2", cex.main=1)
text(cap.env4, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env4,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env4, groups =meta4$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()


png("cap_envs_c.png",width=6,height=4,units="in",res=1200)
library(viridis)
pal<- viridis(12, option = "H") 
pal2<- c( "#30123BFF","#30123BFF","#30123BFF","#4454C4FF","#4454C4FF","#4454C4FF",
                     "#4490FEFF","#4490FEFF","#4490FEFF" ,"#1FC8DEFF","#1FC8DEFF","#1FC8DEFF",
                     "#29EFA2FF","#29EFA2FF","#29EFA2FF", "#7DFF56FF","#7DFF56FF","#7DFF56FF",
                     "#C1F334FF","#C1F334FF","#C1F334FF", "#F1CA3AFF", "#F1CA3AFF", "#F1CA3AFF",
                     "#FE922AFF","#FE922AFF","#FE922AFF" , "#EA4F0DFF", "#EA4F0DFF", "#EA4F0DFF",
                     "#BE2102FF","#BE2102FF","#BE2102FF", "#7A0403FF","#7A0403FF","#7A0403FF")
                     #plots
par(mfrow=c(2,2),mar=c(2, 2, 1, 1))

ordiplot(cap.env1_c,type="n", main="QIIME2", cex.main=1)
text(cap.env1_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env1_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env1_c, groups =meta1$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)
#text(cap.env1_c, dis="species", scaling="sites", cex=0.5, col="black")


ordiplot(cap.env2_c,type="n",main="SINGLE MICOP", cex.main=1)
text(cap.env2_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env2_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env2_c, groups =meta2$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(cap.env3_c,type="n",main="PAIRED MICOP", cex.main=1)
text(cap.env3_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env3_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env3_c, groups =meta3$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)

#text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(cap.env4_c,type="n",main="KRAKEN2", cex.main=1)
text(cap.env4_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.env4_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.env4_c, groups =meta4$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()


png("cca_envs.png",width=6,height=4,units="in",res=1200)
library(viridis)
pal<- viridis(12, option = "H") 
pal2<- c( "#30123BFF","#30123BFF","#30123BFF","#4454C4FF","#4454C4FF","#4454C4FF",
                     "#4490FEFF","#4490FEFF","#4490FEFF" ,"#1FC8DEFF","#1FC8DEFF","#1FC8DEFF",
                     "#29EFA2FF","#29EFA2FF","#29EFA2FF", "#7DFF56FF","#7DFF56FF","#7DFF56FF",
                     "#C1F334FF","#C1F334FF","#C1F334FF", "#F1CA3AFF", "#F1CA3AFF", "#F1CA3AFF",
                     "#FE922AFF","#FE922AFF","#FE922AFF" , "#EA4F0DFF", "#EA4F0DFF", "#EA4F0DFF",
                     "#BE2102FF","#BE2102FF","#BE2102FF", "#7A0403FF","#7A0403FF","#7A0403FF")
                     #plots
par(mfrow=c(2,2),mar=c(2, 2, 1, 1))

ordiplot(vares_cca1,type="n", main="QIIME2", cex.main=1)
text(vares_cca1, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca1,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca1, groups =meta1$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)
#text(vares_cca1, dis="species", scaling="sites", cex=0.5, col="black")


ordiplot(vares_cca2,type="n",main="SINGLE MICOP", cex.main=1)
text(vares_cca2, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca2,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca2, groups =meta2$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(vares_cca3,type="n",main="PAIRED MICOP", cex.main=1)
text(vares_cca3, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca3,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca3, groups =meta3$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)

#text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(vares_cca4,type="n",main="KRAKEN2", cex.main=1)
text(vares_cca4, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca4,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca4, groups =meta4$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()

png("cap_veg.png",width=6,height=4,units="in",res=1200)
library(viridis)
pal<- viridis(12, option = "H") 
pal2<- c( "#30123BFF","#30123BFF","#30123BFF","#4454C4FF","#4454C4FF","#4454C4FF",
                     "#4490FEFF","#4490FEFF","#4490FEFF" ,"#1FC8DEFF","#1FC8DEFF","#1FC8DEFF",
                     "#29EFA2FF","#29EFA2FF","#29EFA2FF", "#7DFF56FF","#7DFF56FF","#7DFF56FF",
                     "#C1F334FF","#C1F334FF","#C1F334FF", "#F1CA3AFF", "#F1CA3AFF", "#F1CA3AFF",
                     "#FE922AFF","#FE922AFF","#FE922AFF" , "#EA4F0DFF", "#EA4F0DFF", "#EA4F0DFF",
                     "#BE2102FF","#BE2102FF","#BE2102FF", "#7A0403FF","#7A0403FF","#7A0403FF")
                     #plots
par(mfrow=c(2,2),mar=c(2, 2, 1,1))

ordiplot(cap.veg1,type="n", main="QIIME2", cex.main=1)
text(cap.veg1, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg1,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg1, groups =meta1$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)
#text(cap.veg1, dis="species", scaling="sites", cex=0.5, col="black")


ordiplot(cap.veg2,type="n", main="SINGLE MICOP", cex.main=1)
text(cap.veg2, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg2,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg2, groups =meta2$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(cap.veg3,type="n", main="PAIRED MICOP", cex.main=1)
text(cap.veg3, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg3,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg3, groups =meta3$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)

#text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(cap.veg4,type="n", main="KRAKEN2", cex.main=1)
text(cap.veg4, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg4,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg4, groups =meta4$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()


png("cap_veg_c.png",width=6,height=4,units="in",res=1200)
library(viridis)
pal<- viridis(12, option = "H") 
pal2<- c( "#30123BFF","#30123BFF","#30123BFF","#4454C4FF","#4454C4FF","#4454C4FF",
                     "#4490FEFF","#4490FEFF","#4490FEFF" ,"#1FC8DEFF","#1FC8DEFF","#1FC8DEFF",
                     "#29EFA2FF","#29EFA2FF","#29EFA2FF", "#7DFF56FF","#7DFF56FF","#7DFF56FF",
                     "#C1F334FF","#C1F334FF","#C1F334FF", "#F1CA3AFF", "#F1CA3AFF", "#F1CA3AFF",
                     "#FE922AFF","#FE922AFF","#FE922AFF" , "#EA4F0DFF", "#EA4F0DFF", "#EA4F0DFF",
                     "#BE2102FF","#BE2102FF","#BE2102FF", "#7A0403FF","#7A0403FF","#7A0403FF")
                     #plots
par(mfrow=c(2,2),mar=c(2, 2, 1, 1))

ordiplot(cap.veg1_c,type="n", main="QIIME2", cex.main=1)
text(cap.veg1_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg1_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg1_c, groups =meta1$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)
#text(cap.veg1_c, dis="species", scaling="sites", cex=0.5, col="black")


ordiplot(cap.veg2_c,type="n",main="SINGLE MICOP", cex.main=1)
text(cap.veg2_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg2_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg2_c, groups =meta2$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(cap.veg3_c,type="n",main="PAIRED MICOP", cex.main=1)
text(cap.veg3_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg3_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg3_c, groups =meta3$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)

#text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(cap.veg4_c,type="n",main="KRAKEN2", cex.main=1)
text(cap.veg4_c, dis="cn", scaling="sites", cex=0.6, col="red")
#points(cap.veg4_c,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(cap.veg4_c, groups =meta4$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()


png("cca_veg.png",width=6,height=4,units="in",res=1200)
library(viridis)
pal<- viridis(12, option = "H") 
pal2<- c( "#30123BFF","#30123BFF","#30123BFF","#4454C4FF","#4454C4FF","#4454C4FF",
                     "#4490FEFF","#4490FEFF","#4490FEFF" ,"#1FC8DEFF","#1FC8DEFF","#1FC8DEFF",
                     "#29EFA2FF","#29EFA2FF","#29EFA2FF", "#7DFF56FF","#7DFF56FF","#7DFF56FF",
                     "#C1F334FF","#C1F334FF","#C1F334FF", "#F1CA3AFF", "#F1CA3AFF", "#F1CA3AFF",
                     "#FE922AFF","#FE922AFF","#FE922AFF" , "#EA4F0DFF", "#EA4F0DFF", "#EA4F0DFF",
                     "#BE2102FF","#BE2102FF","#BE2102FF", "#7A0403FF","#7A0403FF","#7A0403FF")
                     #plots
par(mfrow=c(2,2),mar=c(2, 2, 1, 1))

ordiplot(vares_cca1_veg,type="n", main="QIIME2", cex.main=1)
text(vares_cca1_veg, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca1_veg,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca1_veg, groups =meta1$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing1$CAP1,y = vars_choosing1$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing1$tax, font=3)
#text(vares_cca1_veg, dis="species", scaling="sites", cex=0.5, col="black")


ordiplot(vares_cca2_veg,type="n",main="SINGLE MICOP", cex.main=1)
text(vares_cca2_veg, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca2_veg,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca2_veg, groups =meta2$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing2$CAP1,y = vars_choosing2$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing2$tax, font=3)

ordiplot(vares_cca3_veg,type="n",main="PAIRED MICOP", cex.main=1)
text(vares_cca3_veg, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca3_veg,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca3_veg, groups =meta3$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)

#text(x = vars_choosing3$CAP1,y = vars_choosing3$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing3$tax, font=3)

ordiplot(vares_cca4_veg,type="n",main="KRAKEN2", cex.main=1)
text(vares_cca4_veg, dis="cn", scaling="sites", cex=0.6, col="red")
#points(vares_cca4_veg,display="sites",cex=0.5,pc=19,  col = pal2)
ordiellipse(vares_cca4_veg, groups =meta4$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.5, cex=0.7)
#text(x = vars_choosing4$CAP1,y = vars_choosing4$CAP2, 
#    cex=0.8, col="blue", labels=vars_choosing4$tax, font=3)



dev.off()


