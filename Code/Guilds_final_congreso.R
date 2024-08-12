library(readxl)
library(tidyverse)
library(ggpubr)
## Guilds
metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
guild_colors<- read.csv("Data/guilds_color")
trophic_colors<- read.csv("Data/trophic_color")

# fisrt -------------------------------------------------------------------

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
  !Guild=="-" ) %>% mutate(Method="QIIME2")
colnames(guild_qiime2)[1:36]<- names_col
guild_qiime2<- guild_qiime2 %>% rownames_to_column(var = "otu")

guild_single<- read.delim("Data/table_micop_single.guilds.txt", 
                          check.names = F)%>% filter(
  !Guild=="-" ) %>% mutate(Method="MICOP SINGLE")

guild_paired<- read.delim("Data/table_micop_paired.guilds.txt", 
                          check.names = F)%>% filter(
                            !Guild=="-" ) %>% mutate(Method="MICOP PAIRED")
guilds_fungi<- read.delim("Data/table_fungi_again.guilds.txt", 
                          check.names = F, row.names = 1)%>% filter(
                            !Guild=="-" ) %>% mutate(
                            otu=paste0("otu", rownames(.))) %>% select_at(vars(!starts_with("kraken")))

guild_fungi<- read.delim("Data/table_fungi_again.guilds.txt", 
                          check.names = F, row.names = 1) %>% filter(
                            !Guild=="-" ) %>% mutate(
                              otu=paste0("otu", rownames(.))) %>% column_to_rownames(var = "otu")%>% select_if(is.numeric)%>% t() %>% as.data.frame() %>% rownames_to_column(
                                var = "id_sequence") %>% separate(
                                  ., "id_sequence", c("kraken", "fungi", "id_metagenome", "report", "bracken"), 
                                  sep = "_") %>% dplyr::select(-kraken, -fungi, -report, -bracken) %>% full_join(
                                    metadata) %>% dplyr::select(-id_sequence:-Names, -id_metagenome) %>% column_to_rownames(
                                      var = "SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column(var = "otu") %>% inner_join(guilds_fungi)%>% mutate(Method="KRAKEN2")


relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}


guild_all<- rbind(guild_qiime2, guild_paired, guild_single, guild_fungi)

guild_all_qiime2<- guild_all %>% filter(Method=="QIIME2") %>% 
  group_by(Guild) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Guild") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  filter(all>3.8) %>% 
 rownames_to_column(var = "Guild") %>% 
  mutate(Guild= case_when(
      Guild=="Bryophyte Parasite-Dung Saprotroph-Ectomycorrhizal-Fungal Parasite-Leaf Saprotroph-Plant Parasite-Undefined Saprotroph-Wood Saprotroph"~ "Sap.Indefinido-Ectomi-Parásito",
 Guild=="Endophyte-Insect Pathogen"~"Endofito-Patógeno",
 Guild=="Lichenized" ~"Liquenizado",
 Guild=="Wood Saprotroph" ~"Saprótrofo de madera",
  Guild=="Endophyte-Litter Saprotroph-Soil Saprotroph-Undefined Saprotroph"  ~"Endofito-Sap.Idenfinido",
       .default = as.character(Guild))  ) %>% dplyr::select(Guild,QIIME2=all) 
  
guild_all_single<- guild_all %>% filter(Method=="MICOP SINGLE") %>% 
  group_by(Guild) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Guild") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  filter(all>6.5) %>% 
  rownames_to_column(var = "Guild") %>% 
  mutate(Guild= case_when(
    Guild=="Dung Saprotroph-Wood Saprotroph" ~ "Sap. Estiércol y Madera",
    Guild=="Dung Saprotroph-Undefined Saprotroph-Wood Saprotroph"~"Sap.Indefinido",
    Guild=="Ectomycorrhizal"   ~"Ectomicorriza",
    Guild=="Wood Saprotroph" ~"Saprótrofo de madera",
    Guild=="Plant Pathogen-Wood Saprotroph"   ~"Patógeno-Sap.Planta",
    .default = as.character(Guild))  )%>% dplyr::select(Guild,"SINGLE MICOP"=all) 


guild_all_paired<- guild_all %>% filter(Method=="MICOP PAIRED") %>% 
  group_by(Guild) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Guild") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  filter(all>1) %>% 
  rownames_to_column(var = "Guild") %>% 
  mutate(Guild= case_when(
    Guild=="Epiphyte-Plant Pathogen-Wood Saprotroph"~"Epífita",
    Guild=="Animal Pathogen-Endophyte-Endosymbiont-Epiphyte-Soil Saprotroph-Undefined Saprotroph"
    ~"Patógeno Animal-Sap.Indefinido",
    Guild=="Ectomycorrhizal"   ~"Ectomicorriza",
    Guild=="Wood Saprotroph" ~"Saprótrofo de madera",
    Guild=="Plant Pathogen-Wood Saprotroph"   ~"Patógeno-Sap.Planta",
    .default = as.character(Guild))  ) %>% dplyr::select(Guild,"PAIRED MICOP"=all) 



guild_all_fungi<- guild_all %>% filter(Method=="KRAKEN2") %>% 
  group_by(Guild) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Guild") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  filter(all>4) %>% 
  rownames_to_column(var = "Guild") %>% 
  mutate(Guild= case_when(
    Guild== "Plant Pathogen-Undefined Saprotroph"  ~"Endófito-Patógeno planta",
    Guild=="Animal Pathogen-Endophyte-Endosymbiont-Epiphyte-Soil Saprotroph-Undefined Saprotroph"
    ~"Patógeno Animal-Sap.Indefinido",
    Guild== "Clavicipitaceous Endophyte-Plant Pathogen"      ~"Endófito-Patógeno planta",
    Guild=="Endophyte-Plant Pathogen" ~"Endófito-Patógeno planta",
    Guild=="Animal Pathogen-Endophyte-Plant Pathogen-Soil Saprotroph-Wood Saprotroph"
    ~"Patógeno-Sap.Planta",
    Guild=="Endophyte-Plant Pathogen-Undefined Saprotroph"~"Endófito-Patógeno-Saprótrofo",
    .default = as.character(Guild))  ) %>% dplyr::select(Guild,"KRAKEN2"=all) 



guild_all_qiime2t<- guild_all %>% filter(Method=="QIIME2") %>% 
  group_by(`Trophic Mode`) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Trophic Mode") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  #filter(all>3.8) %>% 
  rownames_to_column(var = "Trophic") %>% 
  mutate(Trophic= case_when(
    Trophic=="Pathotroph-Saprotroph-Symbiotroph"~ "Patótrofo/Sapótrofo/Simbiótrofo",
    Trophic=="Saprotroph"~"Saprótrofo",
    Trophic=="Saprotroph-Symbiotroph"  ~"Saprótrofo/Simbiótrofo",
    Trophic=="Pathotroph"  ~"Patrótrofo",
    Trophic=="Symbiotroph"   ~"Simbiótrofo",
    Trophic=="Pathotroph-Symbiotroph" ~"Patrótrofo/Simbiótrofo",
    Trophic=="Pathotroph-Saprotroph" ~"Patrótrofo/Saprótrofo",
    .default = as.character(Trophic))  ) %>% dplyr::select(Trophic,QIIME2=all) 


guild_all_singlet<- guild_all %>% filter(Method=="MICOP SINGLE") %>% 
  group_by(`Trophic Mode`) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Trophic Mode") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  #filter(all>3.8) %>% 
  rownames_to_column(var = "Trophic") %>% 
  mutate(Trophic= case_when(
    Trophic=="Pathotroph-Saprotroph-Symbiotroph"~ "Patótrofo/Sapótrofo/Simbiótrofo",
    Trophic=="Saprotroph"~"Saprótrofo",
    Trophic=="Saprotroph-Symbiotroph"  ~"Saprótrofo/Simbiótrofo",
    Trophic=="Pathotroph"  ~"Patrótrofo",
    Trophic=="Symbiotroph"   ~"Simbiótrofo",
    Trophic=="Pathotroph-Symbiotroph" ~"Patrótrofo/Simbiótrofo",
    Trophic=="Pathotroph-Saprotroph" ~"Patrótrofo/Saprótrofo",
    .default = as.character(Trophic))  ) %>% dplyr::select(Trophic,"MICOP SINGLE"=all) 



guild_all_pairedt<- guild_all %>% filter(Method=="MICOP PAIRED") %>% 
  group_by(`Trophic Mode`) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Trophic Mode") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  #filter(all>3.8) %>% 
  rownames_to_column(var = "Trophic") %>% 
  mutate(Trophic= case_when(
    Trophic=="Pathotroph-Saprotroph-Symbiotroph"~ "Patótrofo/Sapótrofo/Simbiótrofo",
    Trophic=="Saprotroph"~"Saprótrofo",
    Trophic=="Saprotroph-Symbiotroph"  ~"Saprótrofo/Simbiótrofo",
    Trophic=="Pathotroph"  ~"Patrótrofo",
    Trophic=="Symbiotroph"   ~"Simbiótrofo",
    Trophic=="Pathotroph-Symbiotroph" ~"Patrótrofo/Simbiótrofo",
    Trophic=="Pathotroph-Saprotroph" ~"Patrótrofo/Saprótrofo",
    .default = as.character(Trophic))  ) %>% dplyr::select(Trophic,"MICOP PAIRED"=all) 


guild_all_fungit<- guild_all %>% filter(Method=="KRAKEN2") %>% 
  group_by(`Trophic Mode`) %>% summarise_if(is.numeric, mean)  %>% 
  column_to_rownames(var = "Trophic Mode") %>% relabunda() %>% 
  mutate(all=rowMeans(.)) %>% arrange(-all) %>% 
  #filter(all>3.8) %>% 
  rownames_to_column(var = "Trophic") %>% 
  mutate(Trophic= case_when(
    Trophic=="Pathotroph-Saprotroph-Symbiotroph"~ "Patótrofo/Sapótrofo/Simbiótrofo",
    Trophic=="Saprotroph"~"Saprótrofo",
    Trophic=="Saprotroph-Symbiotroph"  ~"Saprótrofo/Simbiótrofo",
    Trophic=="Pathotroph"  ~"Patrótrofo",
    Trophic=="Symbiotroph"   ~"Simbiótrofo",
    Trophic=="Pathotroph-Symbiotroph" ~"Patrótrofo/Simbiótrofo",
    Trophic=="Pathotroph-Saprotroph" ~"Patrótrofo/Saprótrofo",
    .default = as.character(Trophic))  ) %>% dplyr::select(Trophic,"KRAKEN2"=all) 


all_tro<- guild_all_qiime2t %>% full_join(guild_all_singlet) %>% 
  full_join(guild_all_pairedt) %>% full_join(guild_all_fungit) %>% mutate_if(is.numeric, ~replace_na(., 0))


all_guild<- guild_all_qiime2 %>% full_join(guild_all_single) %>% 
  full_join(guild_all_paired) %>% full_join(guild_all_fungi) %>% mutate_if(is.numeric, ~replace_na(., 0))

write.csv(all_tro, "Data/all_tro.csv")
write.csv(all_guild, "Data/all_guild.csv")


don1<-all_tro %>% dplyr::select(Trophic,QIIME2) %>% mutate_at(c("QIIME2"),~round(.,digits = 0))
don1$fraction <- don1$QIIME2 / sum(don1$QIIME2)
don1$ymax <- cumsum(don1$fraction)
don1$ymin <- c(0, head(don1$ymax, n=-1))
don1$labelPosition <- (don1$ymax + data$ymin) /2
don1$label <- paste0(don1$Trophic, "\n ", don1$QIIME2, "%")
#fix(don1)

ggplot(don1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Trophic)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + # x here controls label position (inner / outer)
  #scale_fill_brewer() +
  # scale_color_brewer(palette="black") +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#barplot_guild<- function(tab,  metadata){
  library(RColorBrewer)
 # n <- 30
  #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #col_vector[6]<-"#02F0CE"
  #col_vector[7]<-"#900C3F"
  #col_vector[4]<-"#0000FF"
library(tidyverse)
library(ggh4x)
table_guild<- guild_all %>%dplyr::select(
    P5S1T3:P6S2T1, Guild, Method) %>% group_by(
            Guild, Method) %>% summarise_if(
              is.numeric, sum) %>% unite(
                "name", Guild:Method) %>%  column_to_rownames(
              var = "name") %>%  mutate(
                all= rowSums(.)) %>% arrange(-all) %>%
    relabunda(.) %>% rownames_to_column(
                  var = "name")%>% #filter(
                    #., !grepl('-', name)) %>% 
  filter(  !name=="NULL_QIIME2") %>% # %>% slice(
                                      #      c(1:30))  %>% 
    pivot_longer(., cols = -name, names_to ="SampleID", 
                 values_to = "relab" ) %>% filter(
                   !SampleID=="all") %>% separate(name, c("Guild", "Method"), sep = "_") %>%  mutate_at(c("Method"), as.factor)
cols_guild<- table_guild %>% inner_join(guild_colors) %>% arrange(Guild)
col_guild <- as.character(cols_guild$col)
names(col_guild) <- as.character(cols_guild$Guild)
  
table_guild$Methods <- factor(table_guild$Method,
                                levels = c("QIIME2","MICOP SINGLE","MICOP PAIRED" , "KRAKEN2"))
  
  
table_guilds<- table_guild %>% inner_join(metadata) %>% group_by(Guild,Method) %>% 
  summarise_if(is.numeric, mean) %>% 
  mutate_at(c("relab"), ~round(., digits = 3)) %>% 
  pivot_wider(names_from = Method, values_from = relab) %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  filter(QIIME2)


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
        x = "Sites", y = "relab",add = "mean",
        facet.by = "Methods", fill="Guild", 
        position = position_fill()) +facet_nested(
          .~Methods, scales = "free_x")+scale_fill_manual(
            name = "Guild",values =col_guild )+
theme_linedraw()+ylab("Relative abundance")+
    xlab("")+theme(legend.text = element_text(face = "plain"))+
  #  guides(fill = guide_legend(nrow = 30))+
    theme(legend.title = element_text(size = 9),
          #axis.ticks = element_blank(),
          legend.text = element_text(size = 9), 
          axis.text.x = element_text(size = 10),
          legend.key.size = unit(0.6, 'cm'), #change legend key size
          legend.key.height = unit(0.45, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'),
          strip.text.x = element_text(size = 16),
          legend.box = "vertical")
  print(barplot_guild)
#}
#barplot_Trophic<- function(tab,  metadata){
  library(RColorBrewer)
  n <- 30
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector[6]<-"#02F0CE"
  col_vector[7]<-"#900C3F"
  col_vector[4]<-"#0000FF"
  library(tidyverse)
  library(ggh4x)
  table_Trophic<- guild_all %>%dplyr::select(P3S1T3:P5S1T2, Trophic=`Trophic Mode`, Method) %>% group_by(
    Trophic, Method) %>% summarise_if(is.numeric, sum) %>% unite(
      "name", Trophic:Method) %>%  column_to_rownames(
        var = "name")%>%  mutate(
        all= rowSums(.)) %>% arrange(-all) %>%
    relabunda(.) %>% rownames_to_column(
      var = "name")%>%# filter(
    #., !grepl('-', Trophic)) %>% 
       filter(!str_detect(name, "^Pathotroph-Saprotroph-Symbiotroph")) %>% # %>% slice(
    #      c(1:30))  %>% 
    pivot_longer(., cols = -name, names_to ="SampleID", 
                 values_to = "relab" ) %>% filter(!SampleID=="all")%>% separate(name, c("Trophic", "Method"), sep = "_") %>%  mutate_at(c("Method"), as.factor)
  cols_trophic<- table_Trophic %>% inner_join(trophic_colors) %>% arrange(Trophic)
  col_trophic <- as.character(cols_trophic$col)
  names(col_trophic) <- as.character(cols_trophic$Trophic)
  table_Trophic$Methods <- factor(table_Trophic$Method,
                                levels = c("QIIME2","MICOP SINGLE","MICOP PAIRED" , "KRAKEN2"))
  
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
        x = "Sites", y = "relab",add = "mean",
        facet.by = "Methods", fill="Trophic",
        position = position_fill()) +
    facet_nested(
          .~Methods, scales = "free_x")+scale_fill_manual(
            values=col_trophic)+theme_linedraw()+ylab("Relative abundance")+
    xlab("")+theme(legend.text = element_text(face = "plain"))+
    #guides(fill = guide_legend(nrow = 30))+
    theme(legend.title = element_text(size = 9),
         # axis.ticks = element_blank(),
          legend.text = element_text(size = 9), 
          axis.text.x = element_text(size = 10),
          legend.key.size = unit(0.6, 'cm'), #change legend key size
          legend.key.height = unit(0.45, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'),
          strip.text.x = element_text(size = 16),
          legend.box = "vertical")+xlab("Sites")
  print(barplot_Trophic)
#}

#barplot_guild_qiime2<- barplot_guild(guild_qiime2, metadata)+ggtitle("QIIME2")
#barplot_guild_single<- barplot_guild(guild_single, metadata)+ggtitle("MICOP SINGLE")+theme(axis.title.y = element_blank())
#barplot_guild_paired<- barplot_guild(guild_paired, metadata)+ggtitle("MICOP PAIRED")+theme(axis.title.y = element_blank())

#barplot_trophic_qiime2<- barplot_Trophic(guild_qiime2, metadata)
#barplot_trophic_single<- barplot_Trophic(guild_single, metadata)+theme(axis.title.y = element_blank())
#barplot_trophic_paired<- barplot_Trophic(guild_paired, metadata)+theme(axis.title.y = element_blank())


#guilds_tro<-cowplot::plot_grid(barplot_guild_qiime2, barplot_guild_single, barplot_guild_paired,
 #                  barplot_trophic_qiime2, barplot_trophic_single, barplot_trophic_paired,
  #                 ncol = 3, nrow = 2)
  library(cowplot)
all<- plot_grid(barplot_guild, barplot_Trophic, nrow = 2)
all
ggsave("exploratory_guilds_final.png",width = 10, height = 7, dpi = 300, plot = all, device = "png")

