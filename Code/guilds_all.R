guild_qiime2<- read.delim(
  "../Data/tabla_qiime2_blast_OTUS.guilds.txt", 
  check.names = F, row.names = 1) %>% rownames_to_column(var = "OTUID") %>% dplyr::select(OTUID, taxonomy:Notes)


guild_single<- read.delim("../Data/table_single.guilds.txt", 
                          check.names = F)%>% dplyr::select(OTUID, taxonomy:Notes)
guild_paired<- read.delim("../Data/table_paired.guilds.txt", 
                          check.names = F)%>% dplyr::select(OTUID, taxonomy:Notes)




guild_fungi<- read.delim("../Data/table_fungi_again.guilds.txt", 
                         check.names = F, row.names = 1) %>% rownames_to_column(var = "OTUID") %>% dplyr::select(OTUID, taxonomy:Notes)


guilds_all<- guild_qiime2 %>% full_join(guild_single) %>% full_join(guild_paired) %>% full_join(guild_fungi) %>% filter(!Taxon=="-") %>% filter(!Guild=="NULL")

write_tsv(guilds_all, "../Data/guilds_all.tsv")
