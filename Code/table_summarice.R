#Alpha diversity
library(tidyverse)
library(readxl)
library(qiime2R)
library(ggpubr)

metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
table_single_micop<- read.delim("Data/table_micop_single.txt") 
table_paired_micop<- read.delim("Data/table_micop_paired.txt")
table_qiime2<- data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
                          check.names = F) %>% t() %>% as.data.frame() %>% rownames_to_column(
                            var = "SampleID") %>% separate(
                              SampleID, c(
                                "id_metagenome", "R", "unmap", "Paired"), 
                              sep = "_")%>% inner_join(metadata) %>% dplyr::select(
                                -id_metagenome:-Paired, -id_sequence:-id_fisicoq) %>% column_to_rownames(
                                  var="SampleID") %>% t() %>% as.data.frame()

#let's summarice the total of features and counts in each method
library(kableExtra)
tables_sum<- list(table_paired_micop,  table_single_micop, table_qiime2)
names_rows<- c("Paired micop",  "Single Micop", "QIIME2")

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

library(readxl)
writexl::write_xlsx(table_summaricing, "Table1. Summarizing.xls")
