barplot_class<- function(tab, taxonomy, metadata){
  library(RColorBrewer)
  n <- 30
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector[6]<-"#02F0CE"
  col_vector[7]<-"#900C3F"
  library(tidyverse)
  library(ggh4x)
  table_class<- tab %>%rownames_to_column(
    var = "Feature.ID") %>%  inner_join(taxonomy) %>% separate(
      Taxon, c("k","p","c","o","f","g","s"), sep = ";" ) %>% mutate_at(
        c("c"), ~str_replace(., "c__", "")) %>% mutate_if(
          is.character, ~replace_na(., "Unassigned")) %>% group_by(
            c) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
              var = "c") %>%  mutate(
                all= rowSums(.)) %>% arrange(-all) %>% relabunda(.) %>% rownames_to_column(
                  var = "Taxon")%>% filter(!Taxon=="unidentified" , !Taxon=="Unassigned") %>% slice(
                    c(1:30))  %>% pivot_longer(
                      ., cols = -Taxon, names_to ="SampleID", 
                      values_to = "relab" ) %>% filter(!SampleID=="all")
  
  barplot_class<- table_class %>% inner_join(metadata) %>% mutate(
    Pol= case_when(
      Poligono==1~ "POL 1",
      Poligono==2~ "POL 2",
      Poligono==3~ "POL 3",
      Poligono==4~ "POL 4",
      Poligono==5~ "POL 5",
      Poligono==6~ "POL 6"),
    Site= case_when(
      Sitio==1~ "S1",
      Sitio==2~ "S2"))%>% ggplot(
        aes(SampleID, relab, fill=Taxon)) +geom_col() +facet_nested(
          .~Pol+Site, scales = "free_x")+scale_fill_manual(values=col_vector)+theme_linedraw()
  print(barplot_class)}
