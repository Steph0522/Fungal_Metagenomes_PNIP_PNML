library(tidyverse)
glomus_seqs<- read.delim("/home/yendi/Downloads/list_glom_seqs_edit.csv")
taxonomy<- read.delim("/home/yendi/Downloads/ncbi_taxa/taxonomy.tsv")

tax<- glomus_seqs %>% left_join(taxonomy) %>% mutate(
  taxonomy=case_when(
    specie=="Ramaria_curta"~"k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Gomphales; f__Gomphaceae; g__Phaeoclavulina; s__curta",
    specie=="Ramaria_stricta" ~"k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Gomphales; f__Gomphaceae; g__Ramaria; s__stricta",
    specie== "Ramaria_aurea"~"k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Gomphales; f__Gomphaceae; g__Ramaria; s__aurea",
    specie== "Ramaria_gracilis"~"k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Gomphales; f__Gomphaceae; g__Ramaria; s__gracilis",
    specie=="Lentaria_byssiseda"~"k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Gomphales; f__Gomphaceae; g__Lentaria; s__byssiseda",
    specie== "Lentaria_aff_micheneri"~"k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Gomphales; f__Gomphaceae; g__Lentaria; s__aff_micheneri",
    TRUE ~ as.character(Taxon)
  )
)
taxonomy_glom<- tax %>% dplyr::select(Feature.ID=ids, Taxon=taxonomy)
write.table(taxonomy_glom,"/home/yendi/Downloads/taxonomy_glom.txt", sep = "\t", row.names = F)
