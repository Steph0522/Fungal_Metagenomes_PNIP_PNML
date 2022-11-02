    table_genus<-table_single_micop %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy_single_micop) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = "__" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", "")) %>% mutate_at(
        c("g"), ~str_replace(., "g__", "")) %>% mutate_if(
          is.character, ~replace_na(., "Unassigned")) %>% group_by(
            g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
              var = "g") %>%  mutate(
                all= rowSums(.)) %>% arrange(-all) %>% relabunda(.) %>% rownames_to_column(
                  var = "Taxon")%>% filter(!Taxon=="unidentified" , !Taxon=="Unassigned") %>% slice(
                    c(1:30))  %>% pivot_longer(
                      ., cols = -Taxon, names_to ="SampleID", 
                      values_to = "relab" ) %>% filter(!SampleID=="all")
#cols<- table_genus %>% inner_join(taxones_color) %>% arrange(Taxon)
#col <- as.character(cols$color)
#names(col) <- as.character(cols$Taxon)}
#qiime2<-table_genus(table_qiime2, taxonomy_qiime2)
ver<-table_genus %>% group_by(Taxon) %>% summarise_if(is.numeric, mean ) %>% arrange(-relab)
