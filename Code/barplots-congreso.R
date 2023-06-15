
#EXPLORATION

library(tidyverse)
library(data.table)
library(tidyverse)
library(readxl)
library(qiime2R)
source("Code/functions_beta.R")
source("Code/functions_compositional.R")

#source("Code/selecting_vars.R")
library(RColorBrewer)
n <- 30
taxones_color<- read_csv("Data/taxones_color.csv")
set.seed(123)
#load and format files
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(
  is.numeric, as.factor) %>% mutate_if(is.character, as.factor) %>% mutate(
    Pol= case_when(
      Poligono==1~ "P1",
      Poligono==2~ "P2",
      Poligono==3~ "P3",
      Poligono==4~ "P4",
      Poligono==5~ "P5",
      Poligono==6~ "P6"),
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
                                 metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names, -Pol, -Site, -id_fisicoq) %>% column_to_rownames(
                                   var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)

table_qiime2_s<- otu.single(table_qiime2)
table_fungi_s<- otu.single(table_fungi)
table_single_micop_s<- otu.single(table_single_micop)
table_paired_micop_s<- otu.single(table_paired_micop)

taxonomy_qiime2<- data.frame(read_qza("Data/taxonomy_blast_dfc_0.98.qza")$data, check.names = F) %>% dplyr::select(Feature.ID,Taxon)
taxonomy_single_micop<- read.delim("Data/table_micop_single.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_paired_micop<- read.delim("Data/table_micop_paired.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_fungi<- read.delim("Data/table_fungi_again.txt", 
                            skip = 1, row.names = 1, check.names = F) %>% select_if(is.character) %>% rownames_to_column(
                              var = "#OTU ID") %>% dplyr::select(Feature.ID = "#OTU ID", Taxon= taxonomy)
#qiime2
table_genus_qiime2<- table_qiime2 %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy_qiime2) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = ";" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", ""))%>% 
  dplyr::mutate(g = stringr::str_trim(g, side = "both")) %>% mutate_if(
    is.character, ~replace_na(., "Unassigned")) %>% group_by(
      g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
        var = "g") %>%  mutate(
          all= rowSums(.)) %>% dplyr::arrange(-all) %>% relabunda(.) %>% rownames_to_column(
            var = "Taxon")%>% filter(!Taxon=="unidentified" ,
                                     !Taxon=="Unassigned") %>% slice(
                                       c(1:30))  %>% pivot_longer(
                                         ., cols = -Taxon, names_to ="SampleID", 
                                         values_to = "relab" ) %>% filter(!SampleID=="all")
cols<- table_genus_qiime2 %>% inner_join(taxones_color) %>% arrange(Taxon)
col <- as.character(cols$color)
names(col) <- as.character(cols$Taxon)

barplot_genus_qiime2<- table_genus_qiime2 %>% inner_join(metadata) %>% ggplot(
  aes(Pol, relab, fill=Taxon)) +geom_col(width = 0.5) +facet_grid(
    .~Pol, scales = "free_x")+
  theme_linedraw()+theme(panel.spacing = unit(0.001, "cm"))+scale_x_discrete(
    labels=rep(1:3,12))+ylab("Relative abundance (%)")+
  xlab("")+theme(legend.text = element_text(face = "italic"))+
  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        legend.box = "vertical")+
  scale_fill_manual(name = "Taxon",values =col )#+scale

#single

table_genus_single<- table_single_micop %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy_single_micop) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = "__" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", ""))%>% 
  dplyr::mutate(g = stringr::str_trim(g, side = "both")) %>% mutate_if(
    is.character, ~replace_na(., "Unassigned")) %>% group_by(
      g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
        var = "g") %>%  mutate(
          all= rowSums(.)) %>% dplyr::arrange(-all) %>% relabunda(.) %>% rownames_to_column(
            var = "Taxon")%>% filter(!Taxon=="unidentified" ,
                                     !Taxon=="Unassigned") %>% slice(
                                       c(1:30))  %>% pivot_longer(
                                         ., cols = -Taxon, names_to ="SampleID", 
                                         values_to = "relab" ) %>% filter(!SampleID=="all")
cols<- table_genus_single %>% inner_join(taxones_color) %>% arrange(Taxon)
col <- as.character(cols$color)
names(col) <- as.character(cols$Taxon)

barplot_genus_single<- table_genus_single %>% inner_join(metadata) %>% ggplot(
  aes(Pol, relab, fill=Taxon)) +geom_col(width = 0.5) +facet_grid(
    .~Pol, scales = "free_x")+
  theme_linedraw()+theme(panel.spacing = unit(0.001, "cm"))+scale_x_discrete(
    labels=rep(1:3,12))+ylab("Relative abundance (%)")+
  xlab("")+theme(legend.text = element_text(face = "italic"))+
  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        legend.box = "vertical")+
  scale_fill_manual(name = "Taxon",values =col )#+scale

#paired

table_genus_paired<- table_paired_micop %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy_paired_micop) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = "__" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", ""))%>% 
  dplyr::mutate(g = stringr::str_trim(g, side = "both")) %>% mutate_if(
    is.character, ~replace_na(., "Unassigned")) %>% group_by(
      g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
        var = "g") %>%  mutate(
          all= rowSums(.)) %>% dplyr::arrange(-all) %>% relabunda(.) %>% rownames_to_column(
            var = "Taxon")%>% filter(!Taxon=="unidentified" ,
                                     !Taxon=="Unassigned") %>% slice(
                                       c(1:30))  %>% pivot_longer(
                                         ., cols = -Taxon, names_to ="SampleID", 
                                         values_to = "relab" ) %>% filter(!SampleID=="all")
cols<- table_genus_paired %>% inner_join(taxones_color) %>% arrange(Taxon)
col <- as.character(cols$color)
names(col) <- as.character(cols$Taxon)

barplot_genus_paired<- table_genus_paired %>% inner_join(metadata) %>% ggplot(
  aes(Pol, relab, fill=Taxon)) +geom_col(width = 0.5) +facet_grid(
    .~Pol, scales = "free_x")+
  theme_linedraw()+theme(panel.spacing = unit(0.001, "cm"))+scale_x_discrete(
    labels=rep(1:3,12))+ylab("Relative abundance (%)")+
  xlab("")+theme(legend.text = element_text(face = "italic"))+
  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        legend.box = "vertical")+
  scale_fill_manual(name = "Taxon",values =col )#+scale

#kraken

table_genus_fungi<- table_fungi %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy_fungi) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = ";" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", ""))%>% 
  dplyr::mutate(g = stringr::str_trim(g, side = "both")) %>% mutate_if(
    is.character, ~replace_na(., "Unassigned")) %>% group_by(
      g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
        var = "g") %>%  mutate(
          all= rowSums(.)) %>% dplyr::arrange(-all) %>% relabunda(.) %>% rownames_to_column(
            var = "Taxon")%>% filter(!Taxon=="unidentified" ,
                                     !Taxon=="Unassigned") %>% slice(
                                       c(1:30))  %>% pivot_longer(
                                         ., cols = -Taxon, names_to ="SampleID", 
                                         values_to = "relab" ) %>% filter(!SampleID=="all")
cols<- table_genus_fungi %>% inner_join(taxones_color) %>% arrange(Taxon)
col <- as.character(cols$color)
names(col) <- as.character(cols$Taxon)

barplot_genus_fungi<- table_genus_fungi %>% inner_join(metadata) %>% ggplot(
  aes(Pol, relab, fill=Taxon)) +geom_col(width = 0.5) +facet_grid(
    .~Pol, scales = "free_x")+
  theme_linedraw()+theme(panel.spacing = unit(0.001, "cm"))+scale_x_discrete(
    labels=rep(1:3,12))+ylab("Relative abundance (%)")+
  xlab("")+theme(legend.text = element_text(face = "italic"))+
  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        legend.box = "vertical")+
  scale_fill_manual(name = "Taxon",values =col )#+scale

bars<-cowplot::plot_grid(barplot_genus_qiime2+ggtitle("QIIME2")+ 
                     theme(strip.text.x = element_text(size = 14, face = "bold")),
                   barplot_genus_single+ggtitle("MICOP SINGLE")+
                     theme(strip.text.x = element_text(size = 14, face = "bold"))+theme(axis.title.y = element_blank()),
                   barplot_genus_paired+ggtitle("MICOP PAIRED")+
                     theme(strip.text.x = element_text(size = 14, face = "bold")),
                   barplot_genus_fungi+ggtitle("KRAKEN2")+
                     theme(strip.text.x = element_text(size = 14, face = "bold"))+theme(axis.title.y = element_blank()), 
                   ncol = 2, nrow = 2, rel_widths = c(1,1,2,1) ,
                   labels = c("A)", "B)", "C)", "D)"), hjust = 0)

bars
ggsave("bars.png",width = 14, height =12, plot = bars, device = "png")
