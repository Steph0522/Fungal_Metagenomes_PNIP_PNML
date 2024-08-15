
library(tidyverse)

table_all<- read.delim("/home/yendi/Downloads/table_all.txt", check.names = F,row.names = 1) %>% dplyr::select(-taxonomy)
metadata<- read.delim("/home/yendi/Downloads/FINALMAP18S", check.names = F) %>% mutate(
  Treatments=case_when(
    Treatment==1 ~ "Wet",
    Treatment==2 ~ "Dry",
    Treatment==3 ~ "Extreme-dry"))
table_ro<- table_all %>% dplyr::select_at(vars(contains("NR")))
#table_ri<- table_all %>% dplyr::select_at(vars(contains("RI")))
#table_nr<- table_all %>% dplyr::select_at(vars(contains("NR")))

table_ro_wet<- table_ro %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "#SampleID") %>% inner_join(metadata) %>% filter(
    Treatments=="Wet") %>% mutate(sampleids=paste0(`#SampleID`, "_wet")) %>% dplyr::select(
      -BarcodeSequence:-Treatments, -"#SampleID") %>% column_to_rownames(
        var = "sampleids") %>% rownames_to_column(var = "sampleids")

table_ro_dry<- table_ro %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "#SampleID") %>% inner_join(metadata) %>% filter(
    Treatments=="Dry") %>% mutate(sampleids=paste0(`#SampleID`, "_dry")) %>% dplyr::select(
      -BarcodeSequence:-Treatments, -"#SampleID") %>% column_to_rownames(
        var = "sampleids") %>% rownames_to_column(var = "sampleids")


table_ro_exdry<- table_ro %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "#SampleID") %>% inner_join(metadata) %>% filter(
    Treatments=="Extreme-dry") %>% mutate(sampleids=paste0(`#SampleID`, "_exdry")) %>% dplyr::select(
      -BarcodeSequence:-Treatments, -"#SampleID") %>% column_to_rownames(
        var = "sampleids") %>% rownames_to_column(var = "sampleids")

write_tsv(table_ro_wet,"/home/yendi/Desktop/table_wet_nr.txt")
write_tsv(table_ro_dry,"/home/yendi/Desktop/table_dry_nr.txt")
write_tsv(table_ro_exdry,"/home/yendi/Desktop/table_exdry_nr.txt")

