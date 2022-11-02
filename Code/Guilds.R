library(readxl)
library(tidyverse)
library(ggpubr)
#Guilds
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
guild_colors<- read_tsv("colors_guild")
trophic_colors<- read_tsv("colors_trophic")

  
names_col<- read.delim(
  "Data/tabla_qiime2_blast_OTUS.guilds.txt", 
  check.names = F, row.names = 1)%>% t() %>% 
  as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% separate(
      SampleID, c(
        "id_metagenome", "R", "unmap", "Paired"), 
      sep = "_")%>% inner_join(meta) %>% dplyr::select(
        SampleID) %>% t() %>% as.vector()
guild_qiime2<- read.delim(
  "Data/tabla_qiime2_blast_OTUS.guilds.txt", 
  check.names = F, row.names = 1) %>% filter(
  !Guild=="-",! `Confidence Ranking`=="Possible" )
colnames(guild_qiime2)[1:36]<- names_col


guild_single<- read.delim("Data/table_micop_single.guilds.txt", 
                          check.names = F)%>% filter(
  !Guild=="-",! `Confidence Ranking`=="Possible" )
guild_paired<- read.delim("Data/table_micop_paired.guilds.txt", 
                          check.names = F)%>% filter(
                            !Guild=="-",! `Confidence Ranking`=="Possible" )

relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}


barplot_guild<- function(tab,  metadata){
  library(RColorBrewer)
 # n <- 30
  #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #col_vector[6]<-"#02F0CE"
  #col_vector[7]<-"#900C3F"
  #col_vector[4]<-"#0000FF"
    library(tidyverse)
  library(ggh4x)
  table_guild<- tab %>%dplyr::select(
    P3S1T3:P5S1T2, Guild) %>% group_by(
            Guild) %>% summarise_if(
              is.numeric, sum) %>% column_to_rownames(
              var = "Guild") %>%  mutate(
                all= rowSums(.)) %>% arrange(-all) %>%
    relabunda(.) %>% rownames_to_column(
                  var = "Guild")%>% filter(
                    ., !grepl('-', Guild)) %>% filter(
                      !Guild=="NULL") %>% # %>% slice(
                                      #      c(1:30))  %>% 
    pivot_longer(., cols = -Guild, names_to ="SampleID", 
                 values_to = "relab" ) %>% filter(
                   !SampleID=="all") 
  cols_guild<- table_guild %>% inner_join(guild_colors) %>% arrange(Guild)
  col_guild <- as.character(cols_guild$col)
  names(col_guild) <- as.character(cols_guild$Guild)
  
  barplot_guild<- table_guild %>% inner_join(metadata) %>% mutate(
    Pol= case_when(
      Poligono==1~ "POL 1",
      Poligono==2~ "POL 2",
      Poligono==3~ "POL 3",
      Poligono==4~ "POL 4",
      Poligono==5~ "POL 5",
      Poligono==6~ "POL 6"),
    Site= case_when(
      Sitio==1~ "S1",
      Sitio==2~ "S2"))%>% ggbarplot(
        x = "Pol", y = "relab",add = "mean",
        facet.by = "Poligono", fill="Guild", 
        position = position_fill()) +facet_nested(
          .~Poligono, scales = "free_x")+scale_fill_manual(
            name = "Guild",values =col_guild )+
theme_linedraw()+ylab("Relative abundance")+
    xlab("")+theme(legend.text = element_text(face = "plain"))+
  #  guides(fill = guide_legend(nrow = 30))+
    theme(legend.title = element_text(size = 8),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 8), 
          axis.text.x = element_blank(),
          legend.key.size = unit(0.6, 'cm'), #change legend key size
          legend.key.height = unit(0.45, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'),
          legend.box = "vertical")
  print(barplot_guild)
}
barplot_Trophic<- function(tab,  metadata){
  library(RColorBrewer)
  n <- 30
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector[6]<-"#02F0CE"
  col_vector[7]<-"#900C3F"
  col_vector[4]<-"#0000FF"
  library(tidyverse)
  library(ggh4x)
  table_Trophic<- tab %>%dplyr::select(P3S1T3:P5S1T2, Trophic=`Trophic Mode`) %>% group_by(
    Trophic) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
      var = "Trophic") %>%  mutate(
        all= rowSums(.)) %>% arrange(-all) %>%
    relabunda(.) %>% rownames_to_column(
      var = "Trophic")%>%# filter(
    #., !grepl('-', Trophic)) %>% 
    filter(
      !Trophic=="NULL", !Trophic=="Pathotroph-Saprotroph-Symbiotroph") %>% # %>% slice(
    #      c(1:30))  %>% 
    pivot_longer(., cols = -Trophic, names_to ="SampleID", 
                 values_to = "relab" ) %>% filter(!SampleID=="all")
  cols_trophic<- table_Trophic %>% inner_join(trophic_colors) %>% arrange(Trophic)
  col_trophic <- as.character(cols_trophic$col)
  names(col_trophic) <- as.character(cols_trophic$Trophic)
  
  barplot_Trophic<- table_Trophic %>% inner_join(
    metadata) %>% mutate(
    Pol= case_when(
      Poligono==1~ "POL 1",
      Poligono==2~ "POL 2",
      Poligono==3~ "POL 3",
      Poligono==4~ "POL 4",
      Poligono==5~ "POL 5",
      Poligono==6~ "POL 6"),
    Site= case_when(
      Sitio==1~ "S1",
      Sitio==2~ "S2"))%>% ggbarplot(
        x = "Pol", y = "relab",add = "mean",
        facet.by = "Poligono", fill="Trophic",
        position = position_fill()) +
    facet_nested(
          .~Poligono, scales = "free_x")+scale_fill_manual(
            values=col_trophic)+theme_linedraw()+ylab("Relative abundance")+
    xlab("")+theme(legend.text = element_text(face = "plain"))+
    #guides(fill = guide_legend(nrow = 30))+
    theme(legend.title = element_text(size = 8),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 8), 
          axis.text.x = element_blank(),
          legend.key.size = unit(0.6, 'cm'), #change legend key size
          legend.key.height = unit(0.45, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'),
          legend.box = "vertical")
  print(barplot_Trophic)
}

barplot_guild_qiime2<- barplot_guild(guild_qiime2, metadata)+ggtitle("QIIME2")
barplot_guild_single<- barplot_guild(guild_single, metadata)+ggtitle("MICOP SINGLE")+theme(axis.title.y = element_blank())
barplot_guild_paired<- barplot_guild(guild_paired, metadata)+ggtitle("MICOP PAIRED")+theme(axis.title.y = element_blank())

barplot_trophic_qiime2<- barplot_Trophic(guild_qiime2, metadata)
barplot_trophic_single<- barplot_Trophic(guild_single, metadata)+theme(axis.title.y = element_blank())
barplot_trophic_paired<- barplot_Trophic(guild_paired, metadata)+theme(axis.title.y = element_blank())


guilds_tro<-cowplot::plot_grid(barplot_guild_qiime2, barplot_guild_single, barplot_guild_paired,
                   barplot_trophic_qiime2, barplot_trophic_single, barplot_trophic_paired,
                   ncol = 3, nrow = 2)
ggsave("exploratory_guilds.png",width = 12, height = 7, dpi = 300, plot = guilds_tro, device = "png")

