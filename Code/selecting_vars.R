

#vars<-table_qiime2 %>% rownames_to_column(
 # var = "Feature.ID") %>% inner_join(taxonomy_qiime2) %>% mutate(
  #  sums = rowSums(table_qiime2[1:32]) ) %>% arrange(-sums)



qiime2_vars<- c("153dcc0b36958edaffccd9da88294bac5f230e3e",
                "19d9e2e72ae907ce2dbf1fc56e57e360e2159e7d",
                "e84be7057e98c927410b8f5e36ff202bedfef793",
                "693a6eddf026a7afc3cffc3b64e4dd1114767c0f")

#vars<-table_single_micop %>% rownames_to_column(
 # var = "Feature.ID") %>% inner_join(taxonomy_single_micop) %>% mutate(
  #  sums = rowSums(table_single_micop[1:32]) ) %>% arrange(-sums)

single_vars<- c("Eukaryota__Basidiomycota__Agaricomycetes__Boletales__Serpulaceae__Serpula__Serpula lacrymans" ,
                "Eukaryota__Basidiomycota__Agaricomycetes__Russulales__Bondarzewiaceae__Heterobasidion__Heterobasidion irregulare",
                "Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Penicillium__Penicillium arizonense"  ,
                "Eukaryota__Ascomycota__Eurotiomycetes__Eurotiales__Aspergillaceae__Penicillium__Penicillium arizonense",
                "Eukaryota__Ascomycota__Leotiomycetes__Unclassified__Pseudeurotiaceae__Pseudogymnoascus__Pseudogymnoascus destructans")


#vars<-table_paired_micop %>% rownames_to_column(
 # var = "Feature.ID") %>% inner_join(taxonomy_paired_micop) %>% mutate(
  #  sums = rowSums(table_paired_micop[1:32]) ) %>% arrange(-sums)

paired_vars<- c("Eukaryota__Ascomycota__Dothideomycetes__Botryosphaeriales__Botryosphaeriaceae__Diplodia__Diplodia corticola" ,
                "Eukaryota__Unclassified__Unclassified__Unclassified__Enterocytozoonidae__Enterocytozoon__Enterocytozoon bieneusi"   ,
                 "Eukaryota__Basidiomycota__Agaricomycetes__Agaricales__Tricholomataceae__Laccaria__Laccaria bicolor",
                "Eukaryota__Basidiomycota__Agaricomycetes__Agaricales__Tricholomataceae__Laccaria__Laccaria bicolor" )


vars<-table_fungi[1:72,]%>% rownames_to_column(
  var = "Feature.ID") %>% inner_join(taxonomy_fungi) %>% mutate(
    sums = rowSums(table_fungi[1:72,]) ) %>% arrange(-sums)

fungi_vars<-c("101028", "2587410", "80884", "98403")
