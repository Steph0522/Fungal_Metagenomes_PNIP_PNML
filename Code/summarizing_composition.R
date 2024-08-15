relabunda<- function(x){(t(t(x)/colSums(x)))*100}
library(tidyverse)
library(readxl)
library(qiime2R)
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
                                 metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names, -Pol, -Site, -id_fisicoq) %>% column_to_rownames(
                                   var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)


taxonomy_qiime2<- data.frame(read_qza("Data/taxonomy_blast_dfc_0.98.qza")$data, check.names = F) %>% dplyr::select(Feature.ID,Taxon)
taxonomy_single_micop<- read.delim("Data/table_micop_single.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_paired_micop<- read.delim("Data/table_micop_paired.txt") %>% rownames_to_column(var = "Feature.ID") %>% dplyr::select(Feature.ID) %>% mutate(Taxon=Feature.ID)
taxonomy_fungi<- read.delim("Data/table_fungi_again.txt", 
                            skip = 1, row.names = 1, check.names = F) %>% select_if(is.character) %>% rownames_to_column(
                              var = "#OTU ID") %>% dplyr::select(Feature.ID = "#OTU ID", Taxon= taxonomy)
#summaricing

sum(table_paired_micop)
sum(table_single_micop)
sum(table_qiime2)
sum(table_fungi)

#number of features
nrow(table_paired_micop)
nrow(table_qiime2)
nrow(table_single_micop)
nrow(table_fungi)

#min and max
min(colSums(table_qiime2))
max(colSums(table_qiime2))
mean(colSums(table_qiime2))
median(colSums(table_qiime2))


min(colSums(table_single_micop))
max(colSums(table_single_micop))
mean(colSums(table_single_micop))
median(colSums(table_single_micop))

min(colSums(table_paired_micop))
max(colSums(table_paired_micop))
mean(colSums(table_paired_micop))
median(colSums(table_paired_micop))


min(colSums(table_fungi))
max(colSums(table_fungi))
mean(colSums(table_fungi))
median(colSums(table_fungi))



table_genus<-table_paired_micop %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy_paired_micop) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = "__" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", "")) %>% mutate_at(
        c("g"), ~str_replace(., "g__", "")) %>% mutate_if(
          is.character, ~replace_na(., "Unassigned")) %>% group_by(
            g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
              var = "g") %>%  mutate(
                all= rowSums(.)) %>% arrange(-all) %>% relabunda(.) %>% as.data.frame() %>% 
  rownames_to_column(var = "Taxon")%>% filter(!Taxon=="unidentified" , !Taxon=="Unassigned") %>% slice(
                    c(1:30))  %>% pivot_longer(
                      ., cols = -Taxon, names_to ="SampleID", 
                      values_to = "relab" ) %>% filter(!SampleID=="all")
#cols<- table_genus %>% inner_join(taxones_color) %>% arrange(Taxon)
#col <- as.character(cols$color)
#names(col) <- as.character(cols$Taxon)}
#qiime2<-table_genus(table_qiime2, taxonomy_qiime2)
table_genus %>% group_by(Taxon) %>% summarise_if(is.numeric, mean ) %>% arrange(-relab)
table_genus %>% group_by(Taxon) %>% summarise_if(is.numeric, sd ) %>% arrange(-relab)
