fungi_brack_table<- read.delim("/home/yendi//Downloads/kraken_fungi_300/fungi/report/bracken_reports/table_fungi_300.txt", skip = 1, check.names = F)
fungi_brack_table<- read.delim("/home/yendi//Downloads/kraken_fungi_300/plus_pf/reports/table_pf_300_fungi.csv",  check.names = F)


otu<- fungi_brack_table %>% dplyr::select_if(is.numeric) %>% column_to_rownames(var = "#OTU ID")
taxa<- fungi_brack_table %>% dplyr::select(Feature.ID="#OTU ID", Taxon="taxonomy")
parse_tax<- parse_taxonomy(taxa) %>% mutate_if(
  is.character, ~replace(., is.na(.), "Unassigned")) %>% rownames_to_column(
    var = "ids") 
parse_tax$ids <- str_replace_all(parse_tax$ids, fixed(" "), "")
parse_tax<- parse_tax %>% column_to_rownames(var = "ids")
gen<- qiime2R::summarize_taxa(otu, parse_tax)$Genus
gen2<- qiime2R::summarize_taxa(otu, parse_tax)$Genus
library(qiime2R)
library(tidyverse)

fungi_bracken<-taxa_barplot(gen)
pluspf_bracken<-taxa_barplot(gen2)
pluspf_bracken
fungi_bracken
