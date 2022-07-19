
#EXPLORATION

library(tidyverse)
library(data.table)
library(tidyverse)
library(readxl)
library(qiime2R)
source("Code/functions_compositional.R")

#load and format files
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
table_single_micop<- read.delim("table_micop_single.txt") 
table_paired_micop<- read.delim("Data/table_micop_paired.txt")
table_qiime2<- data.frame(read_qza("Data/clustered_table_filter.qza")$data, check.names = F) %>% t() %>% as.data.frame() %>% rownames_to_column(
                                var = "SampleID") %>% separate(
                                  SampleID, c(
                                    "id_metagenome", "R", "unmap", "Paired"), 
                                  sep = "_")%>% inner_join(meta) %>% dplyr::select(
                                    -id_metagenome:-Paired, -id_sequence:-id_fisicoq) %>% column_to_rownames(
                                      var="SampleID") %>% t() %>% as.data.frame()
taxonomy_qiime2<- data.frame(read_qza("Data/taxonomy_blast_dfc_0.98.qza")$data, check.names = F) %>% dplyr::select(Feature.ID,Taxon)
taxonomy_single_micop<- read.delim("table_micop_single.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_paired_micop<- read.delim("Data/table_micop_paired.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)

#barplots
barplot_qiime2<- barplot_genus(table_qiime2, taxonomy_qiime2, metadata)+ggtitle("QIIME2")
barplot_micop_single<- barplot_genus2(table_single_micop, taxonomy_single_micop, metadata)+ggtitle("MICOP SINGLE")
barplot_micop_paired<- barplot_genus2(table_paired_micop, taxonomy_paired_micop, metadata)+ggtitle("MICOP PAIRED")

#permanova's
perm_qiime2<- permanova_compo(table_qiime2)
perm_micop_single<- permanova_compo(table_single_micop)
perm_micop_paired<- permanova_compo(table_paired_micop)

#pca's
pca_qiime2<- pca_compositional(table_qiime2)
pca_micop_single<- pca_compositional(table_single_micop)
pca_micop_paired<- pca_compositional(table_paired_micop)

pc1_qiime2<- PC1.f(pca_qiime2)
pc2_qiime2<- PC2.f(pca_qiime2)

pc1_single<- PC1.f(pca_micop_single)
pc2_single<- PC2.f(pca_micop_single)

pc1_paired<- PC1.f(pca_micop_paired)
pc2_paired<- PC2.f(pca_micop_paired)

pca_fig_qiime2<- pca_plot(pca_qiime2, 200, taxonomy_qiime2, metadata)+
  labs(title = paste("Polygon = ",signif(perm_qiime2$`Pr(>F)`[1], 5), ",",
                     "Site = ",signif(perm_qiime2$`Pr(>F)`[2], 5), ",",
                     "Polygon*Site = ",signif(perm_qiime2$`Pr(>F)`[3], 5)))+
  xlab(pc1_qiime2)+ylab(pc2_qiime2)+theme(legend.position = "none", title = element_text(size = 14),axis.title = element_text(size = 16))
pca_fig_micop_single<- pca_plot(
  pca_micop_single, 50, taxonomy_single_micop, metadata)+
  labs(title = paste("Polygon = ",signif(perm_micop_single$`Pr(>F)`[1], 5), ",",
                     "Site = ",signif(perm_micop_single$`Pr(>F)`[2], 5), ",",
                     "Polygon*Site = ",signif(perm_micop_single$`Pr(>F)`[3], 5)))+
  xlab(pc1_single)+ylab(pc2_single)+theme(legend.position = "none", title = element_text(size=14),axis.title = element_text(size = 16))
pca_fig_micop_paired<- pca_plot(
  pca_micop_paired, 50, taxonomy_paired_micop, metadata)+
  labs(title = paste("Polygon = ",signif(perm_micop_paired$`Pr(>F)`[1], 5), ",",
                     "Site = ",signif(perm_micop_paired$`Pr(>F)`[2], 5), ",",
                     "Polygon*Site = ",signif(perm_micop_paired$`Pr(>F)`[3], 5)))+
  xlab(pc1_paired)+ylab(pc2_paired)+theme(title = element_text(size = 14))+theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20),axis.title = element_text(size = 16),
                                                                        legend.key.size = unit(2, "cm"))+ guides(fill = guide_legend(override.aes = list(size = 8)))

library(cowplot)
first<-plot_grid(barplot_qiime2, barplot_micop_single, barplot_micop_paired, ncol = 3)
legds<- get_legend(pca_fig_micop_paired)
second<- plot_grid(pca_fig_qiime2, pca_fig_micop_single,
                   pca_fig_micop_paired+theme(legend.position = "none"), legds, 
                   ncol = 4, rel_widths = c(1,1,1,0.2))
#third<- plot_grid(perm_qiime2, perm_micop_single, perm_micop_paired, ncol = 3)
all<- plot_grid(first, second,  nrow = 2)
ggsave("exploratory_beta.pdf",width = 24, height = 14, dpi = 300, plot = all, device = "pdf")

