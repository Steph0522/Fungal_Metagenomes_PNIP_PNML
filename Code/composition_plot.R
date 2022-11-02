
#EXPLORATION

library(tidyverse)
library(data.table)
library(tidyverse)
library(readxl)
library(qiime2R)
source("Code/functions_compositional.R")
library(RColorBrewer)
n <- 30
taxones_color<- read_csv("taxones_color.csv")
#load and format files
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(
  is.numeric, as.factor) %>% mutate_if(is.character, as.factor) %>% mutate(
    Pol= case_when(
      Poligono==1~ "POL 1",
      Poligono==2~ "POL 2",
      Poligono==3~ "POL 3",
      Poligono==4~ "POL 4",
      Poligono==5~ "POL 5",
      Poligono==6~ "POL 6"),
    Site= case_when(
      Sites==1~ "S1",
      Sites==2~ "S2",
      Sites==3~ "S3",
      Sites==4~ "S4",
      Sites==5~ "S5",
      Sites==6~ "S6",
      Sites==7~ "S7",
      Sites==8~ "S8",
      Sites==9~ "S9",
      Sites==10~ "S10",
      Sites==11~ "S11",
      Sites==12~ "S12"))
metadata$Site <- factor(
  metadata$Site, 
  levels = c("S1","S2","S3","S4",
             "S5","S6","S7","S8",
             "S9","S10","S11","S12"))
meta<- metadata
table_single_micop<- read.delim("Data/table_micop_single.txt") 
table_paired_micop<- read.delim("Data/table_micop_paired.txt")
table_qiime2<- data.frame(
  read_qza("Data/clustered_table_filter.qza")$data,
  check.names = F) %>% t() %>% as.data.frame(
      ) %>% rownames_to_column(
        var = "SampleID") %>% separate(
          SampleID, c(
            "id_metagenome", "R", "unmap", "Paired"), 
          sep = "_")%>% inner_join(
            meta) %>% dplyr::select(
                                    -id_metagenome:-Paired, -id_sequence:-id_fisicoq, -Sites, -Names, -Pol, -Site) %>% column_to_rownames(
                                      var="SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)
table_fungi<- read.delim("Data/table_fungi_again.txt", 
                         skip = 1, row.names = 1, check.names = F) %>% dplyr::select_if(
                           is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "id_sequence") %>% separate(
    ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
    sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
      metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new) %>% column_to_rownames(
        var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)


taxonomy_qiime2<- data.frame(read_qza("Data/taxonomy_blast_dfc_0.98.qza")$data, check.names = F) %>% dplyr::select(Feature.ID,Taxon)
taxonomy_single_micop<- read.delim("Data/table_micop_single.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_paired_micop<- read.delim("Data/table_micop_paired.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_fungi<- read.delim("Data/table_fungi_again.txt", 
                            skip = 1, row.names = 1, check.names = F) %>% select_if(is.character) %>% rownames_to_column(
  var = "#OTU ID") %>% dplyr::select(Feature.ID = "#OTU ID", Taxon= taxonomy)


#barplots
barplot_qiime2<- barplot_genus(table_qiime2, taxonomy_qiime2, metadata)+ggtitle("QIIME2")
barplot_micop_single<- barplot_genus2(table_single_micop, taxonomy_single_micop, metadata)+ggtitle("MICOP SINGLE")
barplot_micop_paired<- barplot_genus2(table_paired_micop, taxonomy_paired_micop, metadata)+ggtitle("MICOP PAIRED")
barplot_fungi<- barplot_genus(table_fungi, taxonomy_fungi, metadata)+ggtitle("QIIME2")


#transform clr data
clr_qiime2<- transform_clr(table_qiime2)
clr_micop_single<- transform_clr(table_single_micop)
clr_micop_paired<- transform_clr(table_paired_micop)
clr_fungi<- transform_clr(table_fungi)

#permanova's
perm_qiime2<- permanova_compo(clr_qiime2, meta)
perm_micop_single<- permanova_compo(clr_micop_single, meta)
perm_micop_paired<- permanova_compo(clr_micop_paired, meta)
perm_fungi<- permanova_compo(clr_fungi, meta)


#pca's
#pca_qiime21<- pca_compositional(table_qiime2)
#pca_micop_single<- pca_compositional(table_single_micop)

pca_qiime2<- pca_compo(clr_qiime2)
pca_micop_single<- pca_compo(clr_micop_single)
pca_micop_paired<- pca_compo(clr_micop_paired)
pca_fungi<- pca_compo(clr_fungi)

pc1_qiime2<- PC1.f(pca_qiime2)
pc2_qiime2<- PC2.f(pca_qiime2)

pc1_single<- PC1.f(pca_micop_single)
pc2_single<- PC2.f(pca_micop_single)

pc1_paired<- PC1.f(pca_micop_paired)
pc2_paired<- PC2.f(pca_micop_paired)

pc1_fungi<- PC1.f(pca_fungi)
pc2_fungi<- PC2.f(pca_fungi)

pca_fig_qiime2<- pca_compositional_sites(pca_qiime2, scales = 250, taxonomy_qiime2)+
  labs(title = paste("F = ",signif(perm_qiime2$F[1], 3), ",",
                     "p-value = ",signif(perm_qiime2$`Pr(>F)`[1], 5)))+
  xlab(pc1_qiime2)+ylab(pc2_qiime2)+theme(legend.position = "none", title = element_text(size = 14),axis.title = element_text(size = 16))

pca_fig_micop_single<- pca_compositional_sites(
  pca_micop_single, 50, taxonomy_single_micop)+
  labs(title = paste("F = ",signif(perm_micop_single$F[1], 3), ",",
                     "p-value = ",signif(perm_micop_single$`Pr(>F)`[1], 5)))+
  xlab(pc1_single)+ylab(pc2_single)+theme(legend.position = "none", title = element_text(size=14),axis.title = element_text(size = 16))
pca_fig_micop_paired<- pca_compositional_sites(
  pca_micop_paired, 50, taxonomy_paired_micop)+
  labs(title = paste("F = ",signif(perm_micop_paired$F[1], 3), ",",
                     "p-value = ",signif(perm_micop_paired$`Pr(>F)`[1], 5)))+
    xlab(pc1_paired)+ylab(pc2_paired)+theme(title = element_text(size = 14))+theme(
      legend.text = element_text(size = 20), legend.title = element_text(size = 20),
 axis.title = element_text(size = 16))+guides(
   fill = guide_legend(verride.aes = list(size = 8), 
   shape = guide_legend(override.aes = list(size = 8))))
pca_fig_fungi<- pca_compositional_sites(
  pca_fungi, 3, taxonomy_fungi)+
  labs(title = paste("F = ",signif(perm_fungi$F[1], 3), ",",
                     "p-value = ",signif(perm_fungi$`Pr(>F)`[1], 5)))+
  xlab(pc1_fungi)+ylab(pc2_fungi)+theme(title = element_text(size = 14))+theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20),axis.title = element_text(size = 16),
                                                                                 legend.key.size = unit(2, "cm"))+ guides(fill = guide_legend(override.aes = list(size = 8)))

library(cowplot)
first<-plot_grid(barplot_qiime2, barplot_micop_single, barplot_micop_paired,barplot_fungi ,ncol = 4 )
legds<- get_legend(pca_fig_micop_paired)
second<- plot_grid(pca_fig_qiime2, pca_fig_micop_single,
                   pca_fig_micop_paired+theme(legend.position = "none"),
                   pca_fig_fungi+theme(legend.position = "none"),legds, 
                   ncol = 5, rel_widths = c(1,1,1,1,0.2), align = "hv", axis = "r")
#third<- plot_grid(perm_qiime2, perm_micop_single, perm_micop_paired, ncol = 3)
all<- plot_grid(first, second,  nrow = 2)
all
ggsave("third.pdf",width = 12, height =7, dpi = 300, plot = third, device = "pdf")

first2<-plot_grid(barplot_qiime2, barplot_micop_single, barplot_micop_paired,barplot_fungi ,ncol = 2, nrow = 2 )
second2<- plot_grid(pca_fig_qiime2+theme(aspect.ratio =6/10), pca_fig_micop_single+theme(aspect.ratio =6/10),
                    pca_fig_micop_paired+theme(legend.position = "none")+theme(aspect.ratio =6/10),
                   pca_fig_fungi+theme(legend.position = "none")+theme(aspect.ratio =6/10),
                   ncol = 2, nrow = 2, rel_widths = c(1,1,1,1), align = "hv")
legd<- plot_grid(NULL, legds, NULL, nrow = 3, ncol = 1)
third<- plot_grid(second2, legd, ncol = 2, rel_widths = c(1,0.2), align = "hv")
third
