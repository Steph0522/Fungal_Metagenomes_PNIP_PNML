library(readxl)
library(cowplot)  
library(agricolae)


metadata_lluvias<- read_excel("/home/yendi/Downloads/Metadatos.xlsx", sheet = "lluvias-septiembre")
metadata_secas<- read_excel("/home/yendi/Downloads/Metadatos.xlsx", sheet = "secas-marzo")

fq_lluvias<- read_excel("Data/fisicoq.xlsx", sheet = "lluvia")
fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
fq_secas2<- read.csv("/home/yendi/Downloads/fisicoq-la.csv")


library(tidyverse)
meta_fq_lluvias<- metadata_lluvias %>% full_join(fq_lluvias) %>% select(SampleID:id_new,pH, MO, N, P) %>% mutate(Season="Rainy")
meta_fq_secas<- metadata_secas %>% full_join(fq_secas) %>%  select(SampleID:id_new,pH, MO, N, P) %>% mutate(Season="Dry")

meta_season<- rbind(meta_fq_lluvias, meta_fq_secas)

shapiro.test(meta_season$pH)#normal
shapiro.test(meta_season$MO)#no normal
shapiro.test(meta_season$P)#no normal
shapiro.test(meta_season$N)#no normal

library(ggpubr)
  a1<-meta_season %>% ggboxplot(., x = "Sitio", y = "pH", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
  a2<-meta_season %>% ggboxplot(., x = "Sitio", y = "N", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )#+stat_compare_means(aes(group = Season), label = "p.signif")
  a3<-meta_season %>% ggboxplot(., x = "Sitio", y = "MO", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )#+stat_compare_means(aes(group = Season), label = "p.signif")
  a4<-meta_season %>% ggboxplot(., x = "Sitio", y = "P", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme(  text = element_text(size=20) )#+stat_compare_means(aes(group = Season), label = "p.signif")
  legends1<- get_legend(a1)
  a<-plot_grid(legends1,a1+theme(legend.position = "none"), a2+theme(legend.position = "none"),a3+theme(legend.position = "none"),a4+theme(legend.position = "none"), nrow = 5, rel_heights = c(0.3,1,1,1,1.4))
  
  a
  b1<-meta_season %>% ggboxplot(., x = "Poligono", y = "pH", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1, scales = "free_x")+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )
  b2<-meta_season %>% ggboxplot(., x = "Poligono", y = "N", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1, scales = "free_x")+theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank() , text = element_text(size=20) )
  b3<-meta_season %>% ggboxplot(., x = "Poligono", y = "MO", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1, scales = "free_x")+theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20)  )
  b4<-meta_season %>% ggboxplot(., x = "Poligono", y = "P", fill = "Season", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1, scales = "free_x")+theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20)  )
  legends<- get_legend(b1)
b<-plot_grid(legends,b1+theme(legend.position = "none"), b2,b3,b4, nrow = 5, rel_heights = c(0.1,1,1,1,1))


c1<-meta_season %>% ggboxplot(., x = "Poligono", y = "pH", fill = "Poligono", facet.by = "Season")+facet_wrap(.~Season, nrow = 1, scales = "free_x")+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )
c2<-meta_season %>% ggboxplot(., x = "Poligono", y = "N", fill = "Poligono", facet.by = "Season")+facet_wrap(.~Season, nrow = 1, scales = "free_x")+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )
c3<-meta_season %>% ggboxplot(., x = "Poligono", y = "MO", fill = "Poligono", facet.by = "Season")+facet_wrap(.~Season, nrow = 1, scales = "free_x")+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )
c4<-meta_season %>% ggboxplot(., x = "Poligono", y = "P", fill = "Poligono", facet.by = "Season")+facet_wrap(.~Season, nrow = 1, scales = "free_x")+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20) )
legends2<- get_legend(c1)
c<-plot_grid(legends2,c1+theme(legend.position = "none"), c2+theme(legend.position = "none"),c3+theme(legend.position = "none"),c4+theme(legend.position = "none"), nrow = 5, rel_heights = c(0.3,1,1,1,1))

ggsave("season_final.png", width =8 , height = 10, dpi = 300, plot = c)
ggsave("poligono-season_final.png", width =18 , height = 10, dpi = 300, plot = b)
ggsave("poligono-season-sitio_final.png", width =18 , height = 10, dpi = 300, plot = a)


#los que solo están en seca
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")
d1<-meta_fq_secas_all %>% ggboxplot(., x = "Poligono", y = "K", fill = "Poligono")+theme_linedraw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=16), legend.position = "top" )
d2<-meta_fq_secas_all %>% ggboxplot(., x = "Poligono", y = "Ca", fill = "Poligono")+theme_linedraw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=16), legend.position = "none" )
d3<-meta_fq_secas_all %>% ggboxplot(., x = "Poligono", y = "Mg", fill = "Poligono")+theme_linedraw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=16) , legend.position = "none")
d4<-meta_fq_secas_all %>% ggboxplot(., x = "Poligono", y = "Fe", fill = "Poligono")+theme_linedraw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=16) , legend.position = "none")
d5<-meta_fq_secas_all %>% ggboxplot(., x = "Poligono", y = "Cu", fill = "Poligono")+theme_linedraw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=16), legend.position = "none" )
d6<-meta_fq_secas_all %>% ggboxplot(., x = "Poligono", y = "Mn", fill = "Poligono")+theme_linedraw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=16) , legend.position = "none")
legends3<- get_legend(d1)
d<-plot_grid(legends3, NULL, d1+theme(legend.position = "none"), d2, d3, d4, d5, d6, ncol = 2, nrow = 4, rel_heights = c(0.3,1,1,1))
ggsave("poligono-metales_final.png", width =14 , height = 10, dpi = 300, plot = d)


e1<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "K", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
e2<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "Ca", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
e3<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "Mg", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
e4<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "Fe", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
e5<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "Cu", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
e6<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "Mn", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
legends4<- get_legend(e1)
e<-plot_grid( e1+theme(legend.position = "none"), e2, e3, e4, e5, e6, ncol = 2, nrow = 3, rel_heights = c(1,1,1.2))
ggsave("poligono-sitio-metales_final.png", width =16 , height = 8, dpi = 300, plot = e)

meta_fq_secas_all<- meta_fq_secas_all %>% unite("intera", c("Poligono", "Sitio"), remove = F)

#shapiro.test(meta_fq_secas_all$ARENA)
meta_fq_secas_all$Poligono<- as.factor(meta_fq_secas_all$Poligono)
meta_fq_secas_all$Sitio<- as.factor(meta_fq_secas_all$Sitio)

#moisture
mod<- lm(moisture ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(moisture ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
moist<-print(out$groups) %>% rownames_to_column(var = "ids")

#pH2
mod<- lm(pH2 ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(pH2 ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
pH2<-print(out$groups) %>% rownames_to_column(var = "ids")

#WHC
mod<- lm(WHC ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(WHC ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
whc<-print(out$groups) %>% rownames_to_column(var = "ids")

#condu
mod<- lm(CONDUC ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(CONDUC ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
ec<-print(out$groups) %>% rownames_to_column(var = "ids")

#silt
mod<- lm(ARCILLA ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(ARCILLA ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
ARC<-print(out$groups) %>% rownames_to_column(var = "ids")

#sand
mod<- lm(ARENA ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(ARENA ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
ARE<-print(out$groups) %>% rownames_to_column(var = "ids")

#LIME
mod<- lm(LIMO ~ Poligono*Sitio, data = meta_fq_secas_all)
mod2<- lm(LIMO ~ intera, data = meta_fq_secas_all)

aov_mod<-aov(mod)
summary(aov_mod)
#plot(mod)
out<-HSD.test(mod2,"intera")
LIM<-print(out$groups) %>% rownames_to_column(var = "ids")

#joinning
all<- moist %>% full_join(pH2, by = "ids") %>% full_join(
  ec, by = "ids") %>% full_join(whc, by = "ids") %>% full_join(
    ARC, by = "ids" ) %>% full_join(ARE, by = "ids") %>% full_join(LIM, by = "ids")
    
write_csv(all, "fisic_lab.csv")


f1<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "moisture", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
#f2<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "Moisture", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
f3<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "CONDUC", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )+ylab("EC")#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
f4<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "WHC", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
f5<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "ARCILLA", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), text = element_text(size=20), legend.position = "none" )+ylab("Silt")#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
f5.1<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "ARENA", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( text = element_text(size=20), legend.position = "none" )+ylab("Sand")#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")
f6<-meta_fq_secas_all %>% mutate_at(c("Poligono"), as.factor) %>% ggboxplot(., x = "Sitio", y = "LIMO", fill="Poligono", facet.by = "Poligono")+facet_wrap(.~Poligono, nrow = 1)+theme_bw()+theme( text = element_text(size=20), legend.position = "none" )+ylab("Lime")#+stat_compare_means(aes(group = Season), method = "t.test", label = "p.signif")

legends4<- get_legend(e1)
f<-plot_grid( f1+theme(legend.position = "none"),  f3, f4, f5,f5.1, f6, ncol = 2, nrow = 3, rel_heights = c(1,1,1.2))
f
ggsave("poligono-sitio-fiscoq-laboratorio.png", width =16 , height = 8, dpi = 300, plot = f)

meta_fq_secas_all<- meta_fq_secas_all %>% unite("intera", c("Poligono", "Sitio"), remove = F)



#pca
#pca secas
meta<- meta_fq_secas_all %>% select(SampleID, Poligono, Sitio, Transecto)
file_pca<- meta_fq_secas_all %>% dplyr::select(
  SampleID, pH:moisture, WHC, EC=CONDUC, 
  Silt=ARCILLA, Lime=LIMO, Sand=ARENA) %>% column_to_rownames(var = "SampleID")

pca_fq<- prcomp(file_pca, center = TRUE, scale. = TRUE)
biplot(pca_fq)
vars<-pca_fq$rotation %>% as.data.frame() %>% rownames_to_column(var = "vars")
ind<- pca_fq$x %>% as.data.frame() %>% rownames_to_column(
  var = "SampleID") %>% inner_join(meta)

ind$Sitio<- as.factor(ind$Sitio)
ind$Poligono<- as.factor(ind$Poligono)

pf<-ggplot()+
  geom_point(data = ind, aes(x = PC1, y = PC2, color=Poligono, shape=Sitio), size=4)+ 
  stat_ellipse(data = ind, aes(x = PC1, y = PC2,color=Poligono))+
  geom_segment(data=vars%>%
                 mutate(PC1=PC1*10, PC2=PC2*10), aes(x=0, y = 0, xend=PC1, yend=PC2),
               arrow=arrow(length=unit(0.15,"cm")),
               alpha = 0.75, color = 'black', size= 0.6)+
  ggrepel::geom_label_repel(data = vars%>%
                              mutate(PC1=PC1*10, PC2=PC2*10), aes(x=PC1, y=PC2, label= vars),
                            segment.colour = NA, col = 'black', fill= "#EEEEEE",

                                                        fontface="bold",  box.padding = 0.2, size=4)+theme_linedraw()
pf
ggsave('pca_fisicoq.png',
    width = 6, height = 4, dpi = 300, plot =pf)

#pca lluvias y secas
meta<- meta_fq_lluvias %>% select(SampleID, Poligono, Sitio, Transecto)
file_pca<- meta_fq_lluvias %>% dplyr::select(
  SampleID, pH:P) %>% column_to_rownames(var = "SampleID")
file_pca<- meta_fq_secas_all %>% dplyr::select(
  SampleID, pH:P) %>% column_to_rownames(var = "SampleID")



pca_fq<- prcomp(file_pca, center = TRUE, scale. = TRUE)
biplot(pca_fq)
vars<-pca_fq$rotation %>% as.data.frame() %>% rownames_to_column(var = "vars")
ind<- pca_fq$x %>% as.data.frame() %>% rownames_to_column(
  var = "SampleID") %>% inner_join(meta)

ind$Sitio<- as.factor(ind$Sitio)
ind$Poligono<- as.factor(ind$Poligono)

pf2<-ggplot()+
  geom_point(data = ind, aes(x = PC1, y = PC2, color=Poligono, shape=Sitio), size=4)+ 
  stat_ellipse(data = ind, aes(x = PC1, y = PC2,color=Poligono))+
  geom_segment(data=vars%>%
                 mutate(PC1=PC1*3, PC2=PC2*3), aes(x=0, y = 0, xend=PC1, yend=PC2),
               arrow=arrow(length=unit(0.15,"cm")),
               alpha = 0.75, color = 'black', size= 0.6)+
  ggrepel::geom_label_repel(data = vars%>%
                              mutate(PC1=PC1*3, PC2=PC2*3), aes(x=PC1, y=PC2, label= vars),
                            segment.colour = NA, col = 'black', fill= "#EEEEEE",
                            
                            fontface="bold",  box.padding = 0.2, size=4)+theme_linedraw()
pf2

library(cowplot)
#pf3<-plot_grid(pf, pf2, labels = c("A. Rainy", "B. Dry" ), hjust = -1)

ggsave('pca_fisicoq_secas_y_lluvias.png',
       width = 10, height = 6, dpi = 300, plot =pf3)


#correlaciones?
View(meta_fq_secas_all)
cors<- meta_fq_secas_all %>% dplyr::select(SampleID, pH:Mn, Humidity=moisture, WHC, Ec=CONDUC, Silt=ARCILLA,
                                           Lime=LIMO, Sand=ARENA) %>% column_to_rownames(var = "SampleID")

corss<- cor(cors, method = "pearson")
library(corrplot)
par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
corrplot(corss, type="upper")
#corrplot(corss, method="circle")
#corrplot(corss, method="number")

library(Hmisc)
corr <- rcorr(as.matrix(cors), type=c("pearson"))
#print(corr)
#cor.out <-corr$r
#cor.out2 <-corr$P

#corrplot(corr$r, p.mat=corr$P, method="color", addCoef.col="black")
corrplot(corr$r, p.mat=corr$P, method="color")
