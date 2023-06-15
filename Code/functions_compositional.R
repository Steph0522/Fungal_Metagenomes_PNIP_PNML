#functions
taxones_color<- read.delim("Data/taxones_color.csv",sep = ",")


relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}

barplot_genus<- function(tab, taxonomy, metadata){

  library(tidyverse)
  library(ggh4x)
  table_genus<- tab %>%rownames_to_column(
    var = "Feature.ID") %>%  inner_join(taxonomy) %>% separate(
      Taxon, c("k","p","c","o","f","g","s"), sep = ";" ) %>% mutate_at(
        c("g"), ~str_replace(., "g__", ""))%>% 
    dplyr::mutate(g = stringr::str_trim(g, side = "both")) %>% mutate_if(
          is.character, ~replace_na(., "Unassigned")) %>% group_by(
            g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
              var = "g") %>%  mutate(
              all= rowSums(.)) %>% dplyr::arrange(-all) %>% relabunda(.) %>% rownames_to_column(
                  var = "Taxon")%>% filter(!Taxon=="unidentified" ,
                                           !Taxon=="Unassigned") %>% slice(
                    c(1:30))  %>% pivot_longer(
                    ., cols = -Taxon, names_to ="SampleID", 
                    values_to = "relab" ) %>% filter(!SampleID=="all")
  cols<- table_genus %>% inner_join(taxones_color) %>% arrange(Taxon)
  col <- as.character(cols$color)
  names(col) <- as.character(cols$Taxon)
  
barplot_genus<- table_genus %>% inner_join(metadata) %>% ggplot(
      aes(SampleID, relab, fill=Taxon)) +geom_col() +facet_nested(
        .~Pol, scales = "free_x")+#scale_fill_manual(
          #values=col_vector)+
  theme_linedraw()+scale_x_discrete(
            labels=rep(1:3,12))+ylab("Relative abundance (%)")+
  xlab("")+theme(legend.text = element_text(face = "italic"))+
  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        legend.box = "vertical")+
  scale_fill_manual(name = "Taxon",values =col )#+scale_y_continuous(limits = c(0,100))
print(barplot_genus)
}

barplot_genus2<- function(tab, taxonomy, metadata){
  
  #library(RColorBrewer)
  #n <- 50
  #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #col_vector[6]<-"#02F0CE"
  #col_vector[7]<-"#900C3F"
  library(tidyverse)
  library(ggh4x)
table_genus<- tab %>%rownames_to_column(
  var = "Feature.ID") %>%  inner_join(taxonomy) %>% separate(
    Taxon, c("k","p","c","o","f","g","s"), sep = "__" ) %>% mutate_at(
      c("g"), ~str_replace(., "g__", "")) %>% mutate_at(
        c("g"), ~str_replace(., "g__", "")) %>% mutate_if(
        is.character, ~replace_na(., "Unassigned")) %>% group_by(
          g) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
            var = "g") %>%  mutate(
              all= rowSums(.)) %>% dplyr::arrange(-all) %>% relabunda(.) %>% rownames_to_column(
                var = "Taxon")%>% filter(!Taxon=="unidentified" , !Taxon=="Unassigned") %>% slice(
                  c(1:30))  %>% pivot_longer(
                    ., cols = -Taxon, names_to ="SampleID", 
                    values_to = "relab" ) %>% filter(!SampleID=="all")
cols<- table_genus %>% inner_join(taxones_color) %>% arrange(Taxon)
col <- as.character(cols$color)
names(col) <- as.character(cols$Taxon)

barplot_genus<- table_genus %>% inner_join(metadata) %>% inner_join(taxones_color) %>% ggplot(
      aes(SampleID, relab, fill=Taxon)) +geom_col() +facet_nested(
        .~Pol, scales = "free_x")+#scale_fill_manual(
          #values=col_vector)+
  theme_linedraw()+scale_x_discrete(
            labels=rep(1:3,12))+ylab("Relative abundance (%)")+
  xlab("")+theme(legend.text = element_text(face = "italic"))+
   guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        legend.box = "vertical")+
  scale_fill_manual(name = "Taxon",values =col )
print(barplot_genus)
}


transform_clr<- function(x){
  library(ALDEx2)
  set.seed(123)
  aldex.clr.transform <- aldex.clr(x, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  return(aldex.clr.transform.data)}


pca_compo<- function(x){
  otu_pca<- prcomp(x)}

pca_compositional<- function(x){
  library(ALDEx2)
  set.seed(123)
  aldex.clr.transform <- aldex.clr(x, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  otu_pca<- prcomp(aldex.clr.transform.data)}

PC1.f <- function(pcx){paste("PC1 : ", round(pcx$sdev[1]^2/sum(pcx$sdev^2),3)*100, "%",sep="")}
PC2.f <- function(pcx){paste("PC2 : ", round(pcx$sdev[2]^2/sum(pcx$sdev^2),3)*100, "%",sep="")}

pca_plot<- function(tab, scales, taxonomys, feature){ggplot() +
    geom_segment(data=data.frame(tab$rotation) %>%   #arrows
                   rownames_to_column(var = "Feature.ID")%>%  
                   mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                   top_n(5, a) %>% #keep 10 furthest away points
                   mutate(PC1=PC1*scales, PC2=PC2*scales),
                 aes(x=0, xend=PC1, y=0, yend=PC2),
                 arrow = arrow(length = unit(0.3,"cm")))+
    geom_point(data=data.frame(tab$x) %>% #individuals
                 rownames_to_column(var = "SampleID")%>%
                 left_join(metadata, by = "SampleID"),
               aes(x=PC1, y=PC2, fill=Sites),shape=21, size=4) +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo" )+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical") +
    geom_polygon(data=data.frame(tab$x) %>% #individuals
                   rownames_to_column(var = "SampleID")%>%
                   left_join(metadata, by = "SampleID")%>%
                   drop_na() %>%
                   group_by(Sites) %>% 
                   slice(chull(PC1, PC2)),
                 aes(x=PC1, y=PC2, fill=Sites, color=Sites),
                 alpha = 0.3,
                 show.legend = FALSE)+
    ggrepel::geom_label_repel(data=data.frame(tab$rotation) %>%   #arrows
                                rownames_to_column(var = "Feature.ID")%>%  
                                mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                                top_n(5, a) %>% #keep 10 furthest away points
                                mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
                                  taxonomys)%>% dplyr::select(
                                    Taxon, PC1, PC2, Feature.ID)%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% mutate(
                                    tax= str_extract(Taxon, "[^__]+$")) %>%
                                mutate_at(c("tax"), funs(tax = case_when(
                                  tax=="Fungi" ~ "Unidentified",
                                  tax=="sajor_caju" ~ "Lentinus",
                                  TRUE~as.character(tax)))),
                              aes(x=PC1, y=PC2, label= tax),
                              segment.colour = NA, col = 'black', fill= "#EEEEEE",
                              fontface="italic",  box.padding = 0.2, size=4)}



permanova_compo<- function(table_transformed, metadata){
  tab<-table_transformed %>% as.data.frame() %>% rownames_to_column(
        var = "SampleID") %>% inner_join(metadata)
    library(vegan)
  perm<- how(nperm = 999)
  
perm<-adonis2(table_transformed~Site, data = tab, method = 
                "euclidian", permutations =perm)
  print(perm)}

permdisp_compo<- function(table_transformed, metadata){
  tab<-table_transformed %>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  mat<- dist(table_transformed, method = "euclidean")
  permdisp<-betadisper(mat, tab$Site)
  permdisp2<- permutest( permdisp, permutations = 999)

  print(permdisp2)
}


pca_compositional_sites<-function(pca){
  metadata1<- as.data.frame(pca$x) %>% rownames_to_column(var = "SampleID") %>% 
    inner_join(metadata)
    y<-ggordiplots::gg_ordiplot(pca, metadata1$Sites, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())#+
    #geom_text(data=data.frame(pca$rotation) %>%   #arrows
     #           rownames_to_column(var = "Feature.ID")%>%  
      #          mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
       #         mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
        #          taxonomy)%>% dplyr::select(
         #           Taxon, PC1, PC2, Feature.ID,a)%>%
          #      mutate_at(
           #       c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% filter(
            #        Feature.ID %in% vars) %>% separate(Taxon, c(letters[1:7]),sep = separ),  #keep 10 furthest away points
             # aes(x=PC1, y=PC2, label= f),
              #fontface="italic",  box.padding = 0.5, size=4)
  print(a)}


pca_new<-function(pca, scales, taxonomys, feature){
    metadata1<- as.data.frame(pca$x) %>% rownames_to_column(var = "SampleID") %>% 
    inner_join(metadata)
  y<-ggordiplots::gg_ordiplot(pca, metadata1$Sites, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggrepel::geom_label_repel(data=data.frame(pca$rotation) %>%   #arrows
                                rownames_to_column(var = "Feature.ID")%>%  
                                mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                                top_n(5, a) %>% #keep 10 furthest away points
                                mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
                                  taxonomys)%>% dplyr::select(
                                    Taxon, PC1, PC2, Feature.ID)%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(
                                    .,";s__unidentified", ""))%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(
                                    .,";g__unidentified", "")) %>% mutate(
                                    tax= str_extract(Taxon, "[^__]+$")) %>%
                                mutate_at(c("tax"), funs(tax = case_when(
                                  tax=="Fungi" ~ "Unidentified",
                                  tax=="pseudograminearum"~"Fusarium",
                                  tax=="oryzae"~ "Aspergillus oryzae",
                                  tax=="oreades"~ "",
                                  tax=="solani"~"Rhizoctonia solani",
                                  tax=="romaleae"~"Encephalitozoon romaleae",
                                  tax=="Pseudogymnoascus verrucosus"~"",
                                                                    tax=="sajor_caju" ~ "Lentinus",
                                  TRUE~as.character(tax)))),
                              aes(x=PC1, y=PC2, label= tax),
                              segment.colour = NA, col = 'black', fill= "#EEEEEE",
                              fontface="italic",  box.padding = 0.1, size=4, label.padding = 0.1)
  print(a)}

