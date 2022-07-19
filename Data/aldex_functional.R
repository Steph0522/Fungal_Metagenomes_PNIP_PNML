EC_predicted <- read.delim("/home/yendi/Downloads/pathways_table_all.txt",  check.names = F) %>% dplyr::select(-"#OTU ID") %>% column_to_rownames(var = "taxonomy")
Alpha.t_asv_table<- read.delim("/home/yendi/Downloads/alpha_all.tsv")


table<- EC_predicted %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sampleid")
table_meta<- table %>% inner_join(Alpha.t_asv_table) %>% group_by(
  Type, Treatment, Plot) %>% summarise_if(is.numeric, sum) %>% unite(
    "ids",Type:Plot,  sep = "_") %>% column_to_rownames(var = "ids") %>% dplyr::select(-Loc:-value) %>% t() %>% as.data.frame() %>% mutate_if(
      is.numeric, as.integer)

table_nr<- table_meta %>% dplyr::select_at(vars(contains("Bulk")))
table_ri<- table_meta %>% dplyr::select_at(vars(contains("Rhizosphere")))
table_ro<- table_meta %>% dplyr::select_at(vars(contains("Root")))

table_nr_ri<- cbind(table_nr, table_ri)
table_nr_ro<- cbind(table_nr, table_ro)
table_ro_ri<- cbind(table_ro, table_ri)

conds1<- c(rep("bs", 9), rep("ri", 9))
conds2<- c(rep("bs", 9), rep("ro", 9))
conds3<- c(rep("ro", 9), rep("ri", 9))

library(ALDEx2)
aldex_nr_ri<-aldex(table_nr_ri, conditions = conds1, mc.samples = 1000, test = "t", effect = T, denom = "all")
aldex_nr_ro<-aldex(table_nr_ro, conditions = conds2, mc.samples = 1000, test = "t", effect = T, denom = "all")
aldex_ro_ri<-aldex(table_nr_ri, conditions = conds3, mc.samples = 1000, test = "t", effect = T, denom = "all")

write.table(aldex_nr_ri, "aldex_nr_ri_funct.txt", sep = "\t")
write.table(aldex_nr_ro, "aldex_nr_ro_funct.txt", sep = "\t")
write.table(aldex_ro_ri, "aldex_ro_ri_funct.txt", sep = "\t")

