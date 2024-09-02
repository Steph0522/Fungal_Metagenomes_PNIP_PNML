#Alpha diversity
library(tidyverse)
library(readxl)
library(qiime2R)
library(ggpubr)

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
table_fungi<- read.delim("Data/table_kraken.txt", 
                          row.names = 1, check.names = F) %>% dplyr::select_if(
                           is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                             var = "id_sequence") %>% separate(
                               ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                               sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                                 metadata) %>% dplyr::select(-id_sequence:-Transecto, -id_metagenome, -Sites, -id_new, -Names, -Pol, -Site, -id_fisicoq) %>% column_to_rownames(
                                   var = "SampleID") %>% t() %>% as.data.frame() %>% mutate_all(as.numeric)



#let's summarice the total of features and counts in each method
library(kableExtra)
tables_sum<- list(table_paired_micop,  table_single_micop, table_qiime2, table_fungi)
names_rows<- c("Paired micop",  "Single Micop", "QIIME2", "Kraken2")

summarice_freq <- function(x){
  list(sum(colSums(x)),nrow(x), min(colSums(x)), 
       max(colSums(x)), mean(colSums(x)), median(colSums(x)))}
summaricing_freq<-sapply(tables_sum, summarice_freq )
dd  <-  as.data.frame(matrix(unlist(summaricing_freq), nrow=6)) 
table_summaricing<-dd %>% t() %>% as.data.frame() %>% dplyr::rename( "Total frequency" = "V1", 
                                                                     "Number of features"="V2",
                                                                     "Min frequency of sample"=V3,
                                                                     "Max frequency of sample"=V4,
                                                                     "Mean frequency of sample"=V5, 
                                                                     "Median frequency of sample"=V6) %>% round(., digits = 1)
rownames(table_summaricing)<- names_rows
table_summaricing %>% round()

library(readxl)
writexl::write_xlsx(table_summaricing, "Table1. Summarizing.xls")
