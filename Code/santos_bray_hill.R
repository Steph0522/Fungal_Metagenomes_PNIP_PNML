#STEP 1: filter otutable and reformat data with calculation of diversity and overlap

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(broom)
library(geosphere)
library(reshape2)
source("Code/functions_beta.R")


#Calculating distances between points
coords<- read_csv("/home/yendi/Documents/corredor_scripts/coord.csv")

coords_mat<- coords%>% mutate(P=paste0("P",pol),
                              S=paste0("S", Sitio),
                              T=paste0("T", Transecto)) %>% unite(
  "SampleID", 
  P:T, sep = "") %>% select(
    SampleID, Longitude, Latitude) %>% column_to_rownames(
      var = "SampleID") %>% as.matrix() 


distance<- distm(coords_mat)/1000
colnames(distance)<- rownames(coords_mat)
rownames(distance)<- rownames(coords_mat)

distance_complete<- distance

distance[upper.tri(distance)] <- NA 

distance_dm<-melt(as.matrix(distance), varnames = c(
  "SampleID.x", "SampleID.y")) %>% drop_na() %>% filter(!value==0)
distance_euclidian<- as.matrix(dist(distance_complete))


#Loading data
library(qiime2R)
library(readxl)
map<- read.csv("/home/yendi/Documents/corredor_scripts/coord.csv") %>% mutate_at(
    c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)

table_single_micop<- read.delim("Data/table_micop_single.txt") 
table_paired_micop<- read.delim("Data/table_micop_paired.txt")
table_qiime2<- data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
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

table_fungi<- read.delim("Data/table_fungi_again.txt", 
                         skip = 1, row.names = 1, check.names = F) %>% dplyr::select_if(
                           is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                             var = "id_sequence") %>% separate(
                               ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                               sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                                 metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names,  -id_fisicoq) %>% column_to_rownames(
                                   var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)




otu <- list(table_single_micop, table_paired_micop, table_qiime2, table_fungi)
otu_match<- lapply(otu, otu.match) # matching to map
otu_single<- lapply(otu_match, otu.single) #remove singletons
otu_norm<- lapply(otu_single, otu.norm)#Normalize
bc.dist<- lapply(otu_norm, beta_div_dist)#Calculate Bray-Curtis dissimilarities
bc.dist2<- lapply(otu_norm, beta_div_dist_hill, q=1)

bc.pcoa2<- lapply(bc.dist2, pcoa_all)

bc.pcoa.axes <- lapply(bc.dist2, pcoa_axes) #obtain pcoa axes
bc.pcoa.eigval <- lapply(bc.dist2, pcoa_eigval) #obtain pcoa eigavlues

#STEP 2 FIGURES:PCOA AND DISTANCES
#Filter BC so that only pairwise comparisons within time points are considered. Transform dissimilarities to similarities.
#bc.dist.tidy.filt<- lapply(bc.dist, bc.dist.tidy.filter)
bc.dist.tidy.filt2<- lapply(bc.dist2, bc.dist.tidy.filter.hill)
#Calculate the range of BC similarities for plotting
#max.sim <- max( b.dist.filt$Similarity)#
#min.sim <- min( b.dist.filt$Similarity)#

max.sim<- 1
min.sim<- 0

#pcoas
qiime2_pcoa <- pcoa_plot(bc.pcoa2[[1]])+
  xlab(paste("PCo1 (", bc.pcoa.eigval[[1]]$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", bc.pcoa.eigval[[1]]$Eigval[2], "%)", sep = "")) +
  ggtitle("QIIME2")

single_pcoa <- pcoa_plot(bc.pcoa2[[2]])+
  xlab(paste("PCo1 (", bc.pcoa.eigval[[2]]$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", bc.pcoa.eigval[[2]]$Eigval[2], "%)", sep = "")) +
  ggtitle("SINGLE MICOP")

paired_pcoa <- pcoa_plot(bc.pcoa2[[3]])+
  xlab(paste("PCo1 (", bc.pcoa.eigval[[3]]$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", bc.pcoa.eigval[[3]]$Eigval[2], "%)", sep = "")) +
  ggtitle("PAIRED MICOP")
  
kraken_pcoa <- pcoa_plot(bc.pcoa2[[4]])+
  xlab(paste("PCo1 (", bc.pcoa.eigval[[4]]$Eigval[2], "%)", sep = "")) +
  ylab(paste("PCo2 (", bc.pcoa.eigval[[4]]$Eigval[3], "%)", sep = "")) +
  ggtitle("KRAKEN2")
    
plot_grid(qiime2_pcoa, single_pcoa, paired_pcoa, kraken_pcoa, nrow = 1, align = "hv")


#Correaltion test
#Perform Pearson correlation test and regression to get stats

cor_test<- lapply(bc.dist.tidy.filt2, cor.b)
lm_test<- lapply(bc.dist.tidy.filt2, lm.b)

stats_qiime2 <- data.frame(label = paste("r = ", signif(cor_test[[1]]$estimate,3), 
                                         "p-value = ", signif(cor_test[[1]]$p.value, 3),
                                         "\nslope = ", signif(lm_test[[1]]$estimate, 3)))
stats_single <- data.frame(label = paste("r = ", signif(cor_test[[2]]$estimate,3), 
                                         "p-value = ", signif(cor_test[[2]]$p.value, 3),
                                         "\nslope = ", signif(lm_test[[2]]$estimate, 3)))
stats_paired <- data.frame(label = paste("r = ", signif(cor_test[[3]]$estimate,3), 
                                         "p-value = ", signif(cor_test[[3]]$p.value, 3),
                                         "\nslope = ", signif(lm_test[[3]]$estimate, 3)))
stats_kraken <- data.frame(label = paste("r = ", signif(cor_test[[4]]$estimate,3), 
                                         "p-value = ", signif(cor_test[[4]]$p.value, 3),
                                         "\nslope = ", signif(lm_test[[4]]$estimate, 3)))




second<- plot_grid(first, leg, nrow = 1, ncol = 2, rel_widths = c(1,0.1),
                   align = "hv", axis = "r" )

plot_distance<- lapply(bc.dist.tidy.filt2, distance.plot0)
leg<- get_legend(qiime2_pcoa)
first<-plot_grid(qiime2_pcoa+theme(legend.position = "none"), 
                 single_pcoa+theme(legend.position = "none"), 
                 paired_pcoa+theme(legend.position = "none"),  
                 kraken_pcoa+theme(legend.position = "none"),   
                 nrow = 1, rel_widths = c(1,1,1,1),
                 align = "hv", axis = "r" )
third<- plot_grid(
  plot_distance[[1]]+ ggtitle(stats_qiime2$label)+theme(plot.title = element_text(size = 10)), 
  plot_distance[[2]]+ ggtitle(stats_single$label)+theme(plot.title = element_text(size = 10)), 
  plot_distance[[3]]+ ggtitle(stats_paired$label)+theme(plot.title = element_text(size = 10)), 
  plot_distance[[4]]+  ggtitle(stats_kraken$label)+theme(plot.title = element_text(size = 10)), 
  nrow = 1)
#plot_grid(get_legend(b), NA, b,c, d,e, nrow = 3, rel_widths = c(5,5), labels = c("a", NA, "b","c", "d", "e"), label_size = 15)

all<-plot_grid(leg, first, third,  nrow = 3, axis = "lr",align = "hv",
          rel_heights = c(0.07,1.5,1))

all
ggsave("Figures_final/Fig1.jaccard_distance.png",width = 12, height = 7, dpi = 300, plot = all, device = "png")

library(ecodist)
library(vegan)
mantel_qiime2<-vegan::mantel(bc.dist2[[1]], distance_complete)
mantel_single<-vegan::mantel(bc.dist2[[2]], distance_complete)
mantel_paired<-vegan::mantel(bc.dist2[[3]], distance_complete)
mantel_kraken<-vegan::mantel(bc.dist2[[4]], distance_complete)

mantels<- data.frame(
  Method=c("QIIME2", "SINGLE MICOP", "PAIRED MICOP", "KRAKEN2"),
  r = c(mantel_qiime2[[3]], mantel_single[[3]], 
        mantel_paired[[3]], mantel_kraken[[3]]),
  p = c(mantel_qiime2[[4]], mantel_single[[4]], 
        mantel_paired[[4]], mantel_kraken[[4]])) %>% mutate_at(c("r"), funs(round(.,digits = 2)))
library(ggpubr)
mantel_q1<-mantels %>% ggtexttable(rows = NULL, theme = ttheme("blank"))%>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 2)

mantel_q2<-mantels %>% ggtexttable(rows = NULL, theme = ttheme("blank"))%>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 2)

library(cowplot)
man<-plot_grid(mantel_bc,  mantel_q2, 
          labels = c("a) Bray curtis", "b)  Horn"), nrow = 1)

man
ggsave("Figures_final/man.png",width = 9, height = 2, dpi = 300, plot = man, device = "png")
