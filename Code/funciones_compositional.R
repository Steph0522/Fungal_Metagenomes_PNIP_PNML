#functions

pca_compositional<- function(x){
  library(ALDEx2)
  aldex.clr.transform <- aldex.clr(x, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  otu_pca<- prcomp(aldex.clr.transform.data)}

pca_plot<- function(tab, scales, taxonomys, feature){ggplot() +
      geom_segment(data=data.frame(tab$rotation) %>%   #arrows
                   rownames_to_column(var = "Feature.ID")%>%  
                   mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                   top_n(8, a) %>% #keep 10 furthest away points
                   mutate(PC1=PC1*scales, PC2=PC2*scales),
                 aes(x=0, xend=PC1, y=0, yend=PC2),
                 arrow = arrow(length = unit(0.3,"cm")))+
    geom_point(data=data.frame(tab$x) %>% #individuals
                 rownames_to_column(var = "SampleID")%>%
                 left_join(metadata, by = "SampleID"),
               aes(x=PC1, y=PC2, fill=Poligono),shape=21, size=4) +
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
                   group_by(Poligono) %>% 
                   slice(chull(PC1, PC2)),
                 aes(x=PC1, y=PC2, fill=Poligono, color=Poligono),
                 alpha = 0.3,
                 show.legend = FALSE)+
    ggrepel::geom_label_repel(data=data.frame(tab$rotation) %>%   #arrows
    rownames_to_column(var = "Feature.ID")%>%  
    mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
    top_n(8, a) %>% #keep 10 furthest away points
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
