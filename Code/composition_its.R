
#EXPLORATION

library(tidyverse)
library(data.table)
library(tidyverse)
library(readxl)
library(qiime2R)
source("Code/functions_compositional.R")

#load and format files
metadata<-read_tsv("/home/yendi/Documents/ITS_corredor/maps_its/its_map.txt") %>% mutate_if(is.numeric, as.factor)
table_qiime2<- data.frame(read_qza("/home/yendi/Documents/ITS_corredor/table_merged_paired.qza")$data, check.names = F) 
taxonomy_qiime2<- data.frame(read_qza("/home/yendi/Documents/ITS_corredor/taxonomy_paired_hybrid.qza")$data, check.names = F) %>% dplyr::select(Feature.ID,Taxon)

#barplots
  library(RColorBrewer)
  n <- 30
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector[6]<-"#02F0CE"
  col_vector[7]<-"#900C3F"
  library(tidyverse)
  library(ggh4x)
  table_genus<- table_qiime2%>%rownames_to_column(
    var = "Feature.ID") %>%  inner_join(taxonomy_qiime2) %>% separate(
      Taxon, c("k","p","c","o","f","g","s"), sep = ";" ) %>% mutate_at(
        c("g"), ~str_replace(., "g__", "")) %>% mutate_if(
          is.character, ~replace_na(., "Unassigned")) %>% group_by(
            g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
              var = "g") %>%  mutate(
                all= rowSums(.)) %>% arrange(-all) %>% relabunda(.) %>% rownames_to_column(
                  var = "Taxon")%>% filter(!Taxon=="unidentified" ,
                                           !Taxon=="Unassigned") %>% slice(
                                             c(1:30))  %>% pivot_longer(
                                               ., cols = -Taxon, names_to ="SampleID", 
                                               values_to = "relab" ) %>% filter(!SampleID=="all")
  
  barplot_genus<- table_genus %>% inner_join(metadata) %>% mutate(
    Pol= case_when(
      Poligono=="P1"~ "POL 1",
      Poligono=="P2"~ "POL 2",
      Poligono=="P3"~ "POL 3",
      Poligono=="P4"~ "POL 4",
      Poligono=="P5"~ "POL 5",
      Poligono=="P6"~ "POL 6"))%>% ggplot(
        aes(Pol, relab, fill=Taxon)) +geom_col() +facet_nested(
          .~Pol+Season, scales = "free_x")+scale_fill_manual(
            values=col_vector)+theme_linedraw()+scale_x_discrete(
              labels=rep(1:3,12))+ylab("Relative abundance (%)")+
    xlab("")+theme(legend.text = element_text(face = "italic"))+
    guides(fill = guide_legend(nrow = 30))+
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 12), 
          legend.key.size = unit(0.6, 'cm'), #change legend key size
          legend.key.height = unit(0.45, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'),
          legend.box = "vertical")
  print(barplot_genus)
barplot_genus

  library(ALDEx2)
  set.seed(123)
  aldex.clr.transform <- aldex.clr(table_qiime2, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  otu_pca<- prcomp(aldex.clr.transform.data)
tab<- otu_pca
  pca_p<- ggplot() +
    geom_segment(data=data.frame(tab$rotation) %>%   #arrows
                   rownames_to_column(var = "Feature.ID")%>%  
                   mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                   top_n(8, a) %>% #keep 10 furthest away points
                   mutate(PC1=PC1*500, PC2=PC2*500),
                 aes(x=0, xend=PC1, y=0, yend=PC2),
                 arrow = arrow(length = unit(0.3,"cm")))+
    geom_point(data=data.frame(tab$x) %>% #individuals
                 rownames_to_column(var = "SampleID")%>%
                 left_join(metadata, by = "SampleID"),
               aes(x=PC1, y=PC2, fill=Poligono,shape= Season), size=5) +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Polygon")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_polygon(data=data.frame(tab$x) %>% #individuals
                   rownames_to_column(var = "SampleID")%>%
                   left_join(metadata, by = "SampleID")%>%
                   drop_na() %>%
                   group_by(Poligono) %>% 
                   slice(chull(PC1, PC2)),
                 aes(x=PC1, y=PC2, fill=Poligono, color=Poligono),
                 alpha = 0.3,
                 show.legend = FALSE)+
    ggrepel::geom_text_repel(data=data.frame(tab$rotation) %>%   #arrows
                               rownames_to_column(var = "Feature.ID")%>%  
                               mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                               top_n(8, a) %>% #keep 10 furthest away points
                               mutate(PC1=PC1*500, PC2=PC2*500)%>% left_join(
                                 taxonomy_qiime2)%>% dplyr::select(
                                   Taxon, PC1, PC2, Feature.ID)%>%
                               mutate_at(
                                 c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% mutate(
                                   tax= str_extract(Taxon, "[^__]+$")) %>%
                               mutate_at(c("tax"), funs(tax = case_when(
                                 tax=="Fungi" ~ "Unidentified",
                                 tax=="sajor-caju" ~ "Lentinus",
                                 TRUE~as.character(tax)))),
                             aes(x=PC1, y=PC2, label= tax),
                             segment.colour = NA, col = 'black', fill= "#EEEEEE",
                             fontface="italic",  box.padding = 0.2, size=6)
pca_p
PC1.f <- function(pcx){paste("PC1 : ", round(pcx$sdev[1]^2/sum(pcx$sdev^2),3)*100, "%",sep="")}
PC2.f <- function(pcx){paste("PC2 : ", round(pcx$sdev[2]^2/sum(pcx$sdev^2),3)*100, "%",sep="")}
  set.seed(123)
  library(vegan)
  library(ALDEx2)
  aldex.clr.transform <- aldex.clr(table_qiime2, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1)) 
  meta_just <- data.frame(
    aldex.clr.transform.data, check.names = F) %>% rownames_to_column(
      var = "SampleID") %>% inner_join(
        metadata) 
  perm <- how(nperm = 999)
  permanova_ma <- adonis2(aldex.clr.transform.data~Poligono*Sitio*Season,
                          data = meta_just,
                          method = "euclidian",
                          permutations =perm) 
  print(permanova_ma)
  library(ggpubr)
  Permanova_table <- data.frame(permanova_ma,
                                check.names = F) %>%
    round(., digits = 3) %>%replace(is.na(.), "-")%>% rownames_to_column(
      var="Factor") %>% ggtexttable(., rows = NULL, theme = ttheme("blank", base_size = 8)) %>%
    tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
    tab_add_hline(at.row = c(10), row.side = "bottom", linewidth = 3, linetype = 1)
  print(Permanova_table)
