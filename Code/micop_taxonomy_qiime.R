da<- read.delim("/home/yendi/Documents/FUNGuild/table_micop_single.txt")
dat<- da %>% separate(taxonomy, c("k", "p", "c", "o", "f", "g", "s")) %>% 
  mutate("k1"="k_","p1"="p_","c1"="c_","o1"="o_","f1"="f_","g1"="g_","s1"="s_" ) %>% unite(
., "king",c("k1", "k"))%>% unite(
  ., "phyl",c("p1", "p"))%>% unite(
    ., "cla",c("c1", "c"))%>% unite(
      ., "ord",c("o1", "o"))%>% unite(
        ., "dam",c("f1", "f"))%>% unite(
          ., "gen",c("g1", "g"))%>% unite(
            ., "spe",c("s1", "s")) %>% unite(., "taxonomy",king:spe, sep = ";")


write_tsv(dat,"/home/yendi/Documents/FUNGuild/table_micop_single.txt")
write.table(dat,"/home/yendi/Documents/FUNGuild/table_micop_paired.txt", sep = "\t", row.names = F)
