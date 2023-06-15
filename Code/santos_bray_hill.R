#STEP 1: filter otutable and reformat data with calculation of diversity and overlap

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(broom)
library(geosphere)
library(reshape2)
library(readxl)
source("Code/functions_beta.R")
source("Code/funciones_compositional.R")

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

metadata_secas<- read_excel("Data/Metadatos.xlsx", sheet = "secas-marzo")
fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
fq_secas2<- read.csv("Data/fisicoq-la.csv")

meta_fq_secas<- metadata_secas %>% full_join(fq_secas) #%>%  select(Sites:id_new,pH, MO, N, P) %>% mutate(Season="Dry")
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")

env.16S=meta_fq_secas_all
df<- env.16S%>% dplyr::select(SampleID, pH:Mn, moisture, WHC:CONDUC, ARCILLA:ARENA) 
dfs=data.frame(df[1],scale(df[,2:17], center = T, scale = T)) %>% dplyr::select(
  SampleID,P,K,Ca,Mg,moisture,WHC, LIMO) %>% 
  column_to_rownames(var = "SampleID")

veg.16s=read_tsv("Data/vegetacion_data_norm.txt")
dfs2=veg.16s %>% dplyr::select(SampleID,varschose)

#pca environmental

pca_env<- prcomp(dfs, center = F, scale. = F)

#distancias euclidianas

dfs_dist<- dist(dfs, method = "euclidean")
nut.dist<- dfs_dist %>% as.matrix()
nut.dist[upper.tri(nut.dist)] <- NA 

nut.dist.tidy <- nut.dist %>% 
  melt(as.matrix(distance), varnames = c(
    "SampleID.x", "SampleID.y"), value.name = "EucDistenv") %>% 
  drop_na() %>% filter(!EucDistenv==0) %>% 
  filter(!is.na(EucDistenv)) %>% 
  filter(SampleID.x != SampleID.y) 


#Calculating distances between points
coords<- read_csv("Data/coord.csv")

coords_mat<- coords%>% mutate(P=paste0("P",pol),
                              S=paste0("S", Sitio),
                              T=paste0("T", Transecto)) %>% unite(
  "SampleID", 
  P:T, sep = "") %>% dplyr::select(
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
map<- read.csv("Data/coord.csv") %>% mutate_at(
    c(1,2,3,7), as.factor) %>% mutate(SampleID= paste0("P",pol, "S", Sitio,"T", Transecto ))
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
map<-metadata
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

#table_qiime2=table_qiime2[, match(colnames(table_fungi), colnames(table_qiime2))] %>% as.data.frame()


otu <- list(table_qiime2, table_single_micop, table_paired_micop,  table_fungi)
otu_match<- lapply(otu, otu.match) # matching to map# machear para que queden el mismo orden que la metadata y todos iguales
otu_single<- lapply(otu_match, otu.single) #remove singletons
otu_norm<- lapply(otu_single, otu.norm)#Normalize
bc.dist<- lapply(otu_norm, beta_div_dist)#Calculate Bray-Curtis dissimilarities con la data normalizada
bc.dist2<- lapply(otu_single, beta_div_dist_hill, q=1)#Calculate Horn dissimilarities sin la data normalizada

bc.pcoa<- lapply(bc.dist, pcoa_all) #pcoa de bray
bc.pcoa2<- lapply(bc.dist2, pcoa_all) #pcoa de Horn

bc.pcoa.axes <- lapply(bc.dist2, pcoa_axes) #obtain pcoa axes
bc.pcoa.eigval <- lapply(bc.dist2, pcoa_eigval) #obtain pcoa eigavlues

#STEP 2 FIGURES:PCOA AND DISTANCES
#Filter BC so that only pairwise comparisons within time points are considered. Transform dissimilarities to similarities.
bc.dist.tidy.filt<- lapply(bc.dist, bc.dist.tidy.filter)
bc.dist.tidy.filt2<- lapply(bc.dist2, bc.dist.tidy.filter.hill)

#Calculate the range of BC similarities for plotting
#max.sim <- max( b.dist.filt$Similarity)#
#min.sim <- min( b.dist.filt$Similarity)#

max.sim<- 1
min.sim<- 0

#permanova's
set.seed(124)
perm_qiime2<- permanova_beta(bc.dist2[[1]], metadata = metadata)
perm_micop_single<- permanova_beta(bc.dist2[[2]], metadata = metadata)
perm_micop_paired<- permanova_beta(bc.dist2[[3]], metadata = metadata)
perm_fungi<- permanova_beta(bc.dist2[[4]], metadata = metadata)

perm_qiime22<- permanova_beta(bc.dist[[1]], metadata = metadata)
perm_micop_single2<- permanova_beta(bc.dist[[2]], metadata = metadata)
perm_micop_paired2<- permanova_beta(bc.dist[[3]], metadata = metadata)
perm_fungi2<- permanova_beta(bc.dist[[4]], metadata = metadata)

#betadisper
set.seed(125)
permd_qiime2<- permdisp_beta(bc.dist2[[1]], metadata = metadata)
permd_micop_single<- permdisp_beta(bc.dist2[[2]], metadata = metadata)
permd_micop_paired<- permdisp_beta(bc.dist2[[3]], metadata = metadata)
permd_fungi<- permdisp_beta(bc.dist2[[4]], metadata = metadata)

permd_qiime22<- permdisp_beta(bc.dist[[1]], metadata = metadata)
permd_micop_single2<- permdisp_beta(bc.dist[[2]], metadata = metadata)
permd_micop_paired2<- permdisp_beta(bc.dist[[3]], metadata = metadata)
permd_fungi2<- permdisp_beta(bc.dist[[4]], metadata = metadata)

#pcoas
qiime2_pcoa <- pcoa_plot(bc.pcoa2[[1]])+
  #xlab(paste("PCo1 (", bc.pcoa.eigval[[1]]$Eigval[1], "%)", sep = "")) +
  #ylab(paste("PCo2 (", bc.pcoa.eigval[[1]]$Eigval[2], "%)", sep = "")) +
  ggtitle("QIIME2")+  
  xlab("DIM1")+
  ylab("DIM2")+
  labs(title = paste("adonis: F = ", signif(perm_qiime2$F[1], 3), ",",
                     "p-value = ",signif(perm_qiime2$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_qiime2$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_qiime2$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(1, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))

qiime2_pcoa2 <- pcoa_plot(bc.pcoa[[1]])+
  #xlab(paste("PCo1 (", bc.pcoa.eigval[[1]]$Eigval[1], "%)", sep = "")) +
  #ylab(paste("PCo2 (", bc.pcoa.eigval[[1]]$Eigval[2], "%)", sep = "")) +
  ggtitle("QIIME2")+  
  xlab("DIM1")+
  ylab("DIM2")+
  labs(title = paste("adonis: F = ", signif(perm_qiime22$F[1], 3), ",",
                     "p-value = ",signif(perm_qiime22$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_qiime22$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_qiime22$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(1, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))


single_pcoa <- pcoa_plot(bc.pcoa2[[2]])+
#  xlab(paste("PCo1 (", bc.pcoa.eigval[[2]]$Eigval[2], "%)", sep = "")) +
#  ylab(paste("PCo2 (", bc.pcoa.eigval[[2]]$Eigval[3], "%)", sep = "")) +
  ggtitle("SINGLE MICOP")+
  labs(title = paste("adonis: F = ", signif(perm_micop_single$F[1], 3), ",",
                     "p-value = ",signif(perm_micop_single$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_micop_single$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_micop_single$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(2, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))

single_pcoa2 <- pcoa_plot(bc.pcoa[[2]])+
  #  xlab(paste("PCo1 (", bc.pcoa.eigval[[2]]$Eigval[2], "%)", sep = "")) +
  #  ylab(paste("PCo2 (", bc.pcoa.eigval[[2]]$Eigval[3], "%)", sep = "")) +
  ggtitle("SINGLE MICOP")+
  labs(title = paste("adonis: F = ", signif(perm_micop_single2$F[1], 3), ",",
                     "p-value = ",signif(perm_micop_single2$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_micop_single2$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_micop_single2$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(2, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))

paired_pcoa <- pcoa_plot(bc.pcoa2[[3]])+
 # xlab(paste("PCo1 (", bc.pcoa.eigval[[3]]$Eigval[1], "%)", sep = "")) +
#  ylab(paste("PCo2 (", bc.pcoa.eigval[[3]]$Eigval[2], "%)", sep = "")) +
  xlab("DIM1")+
  ylab("DIM2")+
  ggtitle("PAIRED MICOP")+
  labs(title = paste("adonis: F = ", signif(perm_micop_paired$F[1], 3), ",",
                     "p-value = ",signif(perm_micop_paired$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_micop_paired$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_micop_paired$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(2, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))

  
paired_pcoa2 <- pcoa_plot(bc.pcoa[[3]])+
  # xlab(paste("PCo1 (", bc.pcoa.eigval[[3]]$Eigval[1], "%)", sep = "")) +
  #  ylab(paste("PCo2 (", bc.pcoa.eigval[[3]]$Eigval[2], "%)", sep = "")) +
  xlab("DIM1")+
  ylab("DIM2")+
  ggtitle("PAIRED MICOP")+
  labs(title = paste("adonis: F = ", signif(perm_micop_paired2$F[1], 3), ",",
                     "p-value = ",signif(perm_micop_paired2$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_micop_paired2$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_micop_paired2$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(2, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))




kraken_pcoa <- pcoa_plot(bc.pcoa2[[4]])+
#  scale_x_continuous(limits = c(-0.01,0.01))+
 # scale_y_continuous(limits = c(-0.01,0.01))+
 # xlab(paste("PCo1 (", bc.pcoa.eigval[[4]]$Eigval[2], "%)", sep = "")) +
#  ylab(paste("PCo2 (", bc.pcoa.eigval[[4]]$Eigval[3], "%)", sep = "")) +
  ggtitle("KRAKEN2")+
  labs(title = paste("adonis: F = ", signif(perm_fungi$F[1], 3), ",",
                     "p-value = ",signif(perm_fungi$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_fungi$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_fungi$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(2, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))

kraken_pcoa2 <- pcoa_plot(bc.pcoa[[4]])+
  #  scale_x_continuous(limits = c(-0.01,0.01))+
  # scale_y_continuous(limits = c(-0.01,0.01))+
  # xlab(paste("PCo1 (", bc.pcoa.eigval[[4]]$Eigval[2], "%)", sep = "")) +
  #  ylab(paste("PCo2 (", bc.pcoa.eigval[[4]]$Eigval[3], "%)", sep = "")) +
  ggtitle("KRAKEN2")+
  labs(title = paste("adonis: F = ", signif(perm_fungi2$F[1], 3), ",",
                     "p-value = ",signif(perm_fungi2$`Pr(>F)`[1], 5), 
                     "\n betadisper: F = ", signif(permd_fungi2$tab$F[1], 2), ",",
                     "p-value = ",signif(permd_fungi2$tab$`Pr(>F)`[1], 5)))+theme(
                       legend.text = element_text(size = 20), legend.title = element_text(size = 20),
                       axis.title = element_text(size = 16),
                       plot.title = element_text(hjust = 1, size = 12),
                       legend.key.size = unit(2, "cm"))+ guides(
                         fill = guide_legend(override.aes = list(size = 8)))
    

leg<- get_legend(qiime2_pcoa)

first<-plot_grid(qiime2_pcoa+theme(legend.position = "none")+theme(aspect.ratio =6/10),  
                 single_pcoa+theme(legend.position = "none")+theme(aspect.ratio =6/10), 
                 paired_pcoa+theme(legend.position = "none")+theme(aspect.ratio =6/10),   
                 kraken_pcoa+theme(legend.position = "none")+theme(aspect.ratio =6/10), 
                 ncol = 2, nrow = 2, rel_widths = c(1,1,1,1),
                 align = "v",
                 labels = c("A) QIIME2", "B) SINGLE MICOP",
                            "C) PAIRED MICOP", "D) KRAKEN2"),hjust = 0)
second<-plot_grid(qiime2_pcoa2+theme(legend.position = "none")+theme(aspect.ratio =6/10),  
                 single_pcoa2+theme(legend.position = "none")+theme(aspect.ratio =6/10)+
                   xlab("DIM1")+ylab("DIM2") ,
                 paired_pcoa2+theme(legend.position = "none")+theme(aspect.ratio =6/10),   
                 kraken_pcoa2+theme(legend.position = "none")+theme(aspect.ratio =6/10), 
                 ncol = 2, nrow = 2, rel_widths = c(1,1,1,1),
                 align = "v",
                 labels = c("A) QIIME2", "B) SINGLE MICOP",
                            "C) PAIRED MICOP", "D) KRAKEN2"),hjust = 0)
first

pcoas_plot <- plot_grid(first,leg,ncol = 2, rel_widths = c(.9,.1), align = "hv")
pcoas_plot <- plot_grid(second,leg,ncol = 2, rel_widths = c(.9,.1), align = "hv")

pcoas_plot
ggsave("pcoas_plot_final_bray.png",width =9, height =7, dpi = 300, plot = pcoas_plot, device = "png")
#Correaltion test

library(ecodist)
library(vegan)
mantel_qiime2<-vegan::mantel(bc.dist2[[1]], distance_complete)
mantel_single<-vegan::mantel(bc.dist2[[2]], distance_complete)
mantel_paired<-vegan::mantel(bc.dist2[[3]], distance_complete)
mantel_kraken<-vegan::mantel(bc.dist2[[4]], distance_complete)


mantel_qiime2_e<-vegan::mantel(bc.dist2[[1]], dfs_dist)
mantel_single_e<-vegan::mantel(bc.dist2[[2]], dfs_dist)
mantel_paired_e<-vegan::mantel(bc.dist2[[3]], dfs_dist)
mantel_kraken_e<-vegan::mantel(bc.dist2[[4]], dfs_dist)

cor_test<- lapply(bc.dist.tidy.filt2, cor.b)
lm_test<- lapply(bc.dist.tidy.filt2, lm.b)

cor_test_e<- lapply(bc.dist.tidy.filt2, cor.e)
lm_test_e<- lapply(bc.dist.tidy.filt2, lm.e)

stats_qiime2 <- data.frame(label = paste("Pearson: r = ", signif(cor_test[[1]]$estimate,3), 
                                         ", p-value = ", signif(cor_test[[1]]$p.value, 3),
                                         "\nRegression: slope = ", signif(lm_test[[1]]$estimate, 3),
                           "\nMantel: r = ", signif(mantel_qiime2$statistic,3), 
                           " ,p-value = ", signif(mantel_qiime2$signif, 3)))
stats_single <- data.frame(label = paste("Pearson: r = ", signif(cor_test[[2]]$estimate,3), 
                                         ", p-value = ", signif(cor_test[[2]]$p.value, 3),
                                         "\nRegression: slope = ", signif(lm_test[[2]]$estimate, 3),
                                         "\nMantel: r = ", signif(mantel_single$statistic,3), 
                                         ", p-value = ", signif(mantel_single$signif, 3)))
stats_paired <- data.frame(label = paste("Pearson: r = ", signif(cor_test[[3]]$estimate,3), 
                                         ", p-value = ", signif(cor_test[[3]]$p.value, 3),
                                         "\nRegression: slope = ", signif(lm_test[[3]]$estimate, 3),
                                         "\nMantel: r = ", signif(mantel_paired$statistic,3), 
                                         ", p-value = ", signif(mantel_paired$signif, 3)))
stats_kraken <- data.frame(label = paste("Pearson: r = ", signif(cor_test[[4]]$estimate,3), 
                                         ", p-value = ", signif(cor_test[[4]]$p.value, 3),
                                         "\nRegression: slope = ", signif(lm_test[[4]]$estimate, 3),
                                         "\nMantel: r = ", signif(mantel_kraken$statistic,3), 
                                         ", p-value = ", signif(mantel_kraken$signif, 3)))

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



second<- plot_grid(first, leg, nrow = 1, ncol = 2, rel_widths = c(1,0.1),
                   align = "hv", axis = "r" )

plot_distance<- lapply(bc.dist.tidy.filt2, distance.plot0)
plot_distance_env<- lapply(bc.dist.tidy.filt2, distance.plot1)

third<- plot_grid(
  plot_distance[[1]]+ ggtitle(stats_qiime2$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10),   
  plot_distance[[2]]+ ggtitle(stats_single$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10),  
  plot_distance[[3]]+ ggtitle(stats_paired$label)+theme(plot.title =  element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10), 
  plot_distance[[4]]+  ggtitle(stats_kraken$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10), 
  nrow = 1)

blank_plot<-ggplot() + theme_void()+ggtitle("Environmental distance")
leg1<- get_title(blank_plot)
leg2<- plot_grid(NULL, NULL,leg1, NULL, ncol = 4)

fourth<- plot_grid(
  plot_distance_env[[1]]+  ggtitle(stats_qiime2_e$label)+theme(plot.title= element_text(size = 10,hjust = 0.1,vjust = -100), axis.title.x = element_blank()),
  plot_distance_env[[2]]+ ggtitle(stats_single_e$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),
  plot_distance_env[[3]]+ ggtitle(stats_paired_e$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title.x = element_blank()),  
  plot_distance_env[[4]]+  ggtitle(stats_kraken_e$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),  
  nrow = 2,ncol = 2, labels = c("A) QIIME2", "B) SINGLE MICOP",
                       "C) PAIRED MICOP", "D) KRAKEN2"))

four<- plot_grid(fourth, leg2, nrow = 2, rel_heights = c(1,0.05))
four
#plot_grid(get_legend(b), NA, b,c, d,e, nrow = 3, rel_widths = c(5,5), labels = c("a", NA, "b","c", "d", "e"), label_size = 15)

all<-plot_grid(  fourth,third ,nrow = 2, axis = "lr",align = "hv",
          rel_heights = c(1, 1))

all
#ggsave("Fig1.bray_distance_env.tiff",width = 8, height = 8.2, dpi = 300, plot = four, device = "tiff")



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
man<-plot_grid(mantel_q1,  mantel_q2, 
          labels = c("a) Bray curtis", "b)  Horn"), nrow = 1)


ver<-plot_grid(fourth,fourth, nrow = 2)
man
#ggsave("ver.png",width = 6, height = 10, dpi = 300, plot = ver, device = "png")


#plot of distances

qiime2_data<- bc.dist.tidy.filt2[[1]] %>% mutate(Method="QIIME2") %>% pivot_longer(cols = c("SpatialDistance", "EucDist"), names_to = "Distances", values_to = "val_distance")
single_data<- bc.dist.tidy.filt2[[2]] %>% mutate(Method="MICOP SINGLE")%>% pivot_longer(cols = c("SpatialDistance", "EucDist"), names_to = "Distances", values_to = "val_distance")
paired_data<- bc.dist.tidy.filt2[[3]] %>% mutate(Method="MICOP PAIRED")%>% pivot_longer(cols = c("SpatialDistance", "EucDist"), names_to = "Distances", values_to = "val_distance")
kraken_data<- bc.dist.tidy.filt2[[4]] %>% mutate(Method="KRAKEN2")%>% pivot_longer(cols = c("SpatialDistance", "EucDist"), names_to = "Distances", values_to = "val_distance")


joined_data<- rbind(qiime2_data, single_data, paired_data, kraken_data )
joined_data$Method<- factor(joined_data$Method, levels = c(
  "QIIME2", "MICOP SINGLE", "MICOP PAIRED", "KRAKEN2"))

joined_data$Distances<- factor(joined_data$Distances,
                               levels = c("SpatialDistance","EucDist"), 
                               labels = c("Spatial Distance (km)", "Environmental Distance"))

ann_text<-data.frame(val_distance=c( 30,30,30,30,4.2,4.2,4.2,4.2),
                     Similarity=c(0.8,0.3,0.3,0.3, 0.8,0.3,0.3,0.3),
                     Distances=c("Spatial Distance (km)", 
                                 "Spatial Distance (km)", 
                                 "Spatial Distance (km)",
                                 "Spatial Distance (km)",
                                 "Environmental Distance", 
                                 "Environmental Distance", 
                                 "Environmental Distance",
                                 "Environmental Distance"),
                     Method=c("QIIME2", "MICOP SINGLE", "MICOP PAIRED", "KRAKEN2",
                              "QIIME2", "MICOP SINGLE", "MICOP PAIRED", "KRAKEN2"),
                     label=c(stats_qiime2$label, stats_single$label, stats_paired$label, stats_kraken$label,
                             stats_qiime2_e$label, stats_single_e$label, stats_paired_e$label, stats_kraken_e$label)) 
ann_text$Method<- factor(ann_text$Method, levels = c(
  "QIIME2", "MICOP SINGLE", "MICOP PAIRED", "KRAKEN2"))
ann_text$Distances<- factor(ann_text$Distances,
                            levels = c("Spatial Distance (km)", "Environmental Distance"))


jo<-joined_data %>% ggplot(aes(x = val_distance, y =Similarity ))+
  facet_grid(vars(Method), vars(Distances), scales = "free_x")+
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = F) +
  ylab("Horn similarity") +
  ylim(.2, max.sim) +
  theme_linedraw()+theme(legend.position = "none", 
                         axis.text = element_text(size = 12),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8")) +  
  geom_text(data = ann_text,label=ann_text$label, size=3)+xlab("Distances")
jo
ggsave("Fig.Bray_joined_distances_final.png",width = 6.5, height = 8, dpi = 300, plot = jo, device = "png")

