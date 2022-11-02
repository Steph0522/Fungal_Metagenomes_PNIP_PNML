# Functions of beta diversity

rel_ab <- function(otu, total = 100) {
  t(t(otu)/colSums(otu)) * 100
}

log_norm <- function(otu) {
  log(otu + 1)
}
beta_div_dist <- function(otu, method = "bray") {
  require(phyloseq)
  physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE))
  as.matrix(phyloseq::distance(physeq, method))
} 

beta_div_dist_hill<- function(otu, q=NULL){
  require(hillR)
otu<-otu  
  if (q==0) {
b<- hilldiv::pair_dis(otu, qvalue = 0, metric = "U")
return(b[["L1_UqN"]])} 
if (q==1) {
b<- hilldiv::pair_dis(otu, qvalue = 1, metric = "U")
return(b[["L1_UqN"]])}
if (q==2) {
  b<- hilldiv::pair_dis(otu, qvalue = 2, metric = "C")
return(b[["L1_CqN"]])

      }
  
}

pcoa_all <- function(dist) {
  require(ape)
  require(vegan)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- wcmdscale(as.dist(dist), eig=TRUE)
  return(pcoa)
}

pcoa_axes <- function(dist) {
  require(ape)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  as.data.frame(pcoa$vectors) %>%
    mutate(SampleID = rownames(.)) %>%
    inner_join(map, by = "SampleID")
}
pcoa_eigval <- function (dist) {
  require(ape)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
  data.frame( PC = 1:length(eigval), Eigval = eigval, CumEigval = cumsum(eigval))
}

bc.dist.tidy.filter <-function(bc.dist){ 
  require(reshape2)
  bc.dist[upper.tri(bc.dist)] <- NA 
  bc.dist.tidy<-bc.dist %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
  
  b.dist.filt <- bc.dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value) %>% 
    mutate(Similarity = 1 - Distance)
  return(b.dist.filt)
}


bc.dist.tidy.filter.hill <-function(bc.dist){ 
  require(reshape2)
  bc.dist.tidy<-bc.dist %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
  
  b.dist.filt <- bc.dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value) %>% 
    mutate(Similarity = 1 - Distance)
  return(b.dist.filt)
}

overlap_function<- function(otu){
  #Generate presence/absence OTU table
  bi.otu <- otu > 0
  #Find the number of OTUs overlapping between pairs of samples. 
  #%*% indicates that this should be treated as a matrix multiplication in 
  #linear algebra
  overlap <- t(bi.otu) %*% (bi.otu)
  #The universe of vOTUs shared by a pair of samples is the  sum of the 
  #overlap and the number of vOTus unique to each of the samples.  
  total <- (colSums(bi.otu) - overlap) + t(colSums(bi.otu) - overlap) + overlap
  #Calculate the percent of vOTUs shared
  pct.shared <- overlap / total
  #Remove redundant values
  pct.shared[upper.tri(pct.shared)] <- NA
  #Generate data frame
  pct.tidy <- pct.shared %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y"), value.name = "PercentShared")  %>% 
    filter(!is.na(PercentShared)) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) 
}

pcoa_plot<-function(pca){
  y<-ggordiplots::gg_ordiplot(pca, map$Site, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites", nrow = 1))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "top", 
          legend.box = "horizontal",
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


otu.match <- function(x){x[,match(map$SampleID, colnames(x))]}
otu.single <- function(otu){otu[rowSums(otu>0)>1,]}
otu.norm <- function(otu){otu %>% rel_ab() %>% log_norm()}
#bc.dist <-funct beta_div_dist(otu.norm)
cor.b <- function(x){
  cor.test(x$SpatialDistance,x$Similarity, method= "spearman", alternative = "two.sided") %>% tidy()
}
lm.b <- function(x){lm(Similarity ~ SpatialDistance, data = x) %>% tidy() %>% filter(term == "SpatialDistance")}

distance.plot <- function(x){
  x %>% 
    ggplot(aes(SpatialDistance, Similarity)) +
    geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
    geom_smooth(method = "lm", color = "black", se = F) +
    xlab("Distance between plots (km)") +
    ylab("Bray-Curtis similarity") +
    ylim(.2, max.sim) +
    xlim(0,60)+
    theme_linedraw()+theme(legend.position = "none", 
                           axis.text = element_text(size = 12),
                           strip.text = element_text(size = 12, face = "bold.italic"),
                           panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                           colour = "#E5E8E8"), 
                           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                           colour = "#E5E8E8"))}

distance.plot0 <- function(x){
  x %>% 
    ggplot(aes(SpatialDistance, Similarity)) +
    geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
    geom_smooth(method = "lm", color = "black", se = F) +
    xlab("Distance between plots (km)") +
    ylab("Horn similarity") +
    ylim(.2, max.sim) +
    xlim(0,60)+
    theme_linedraw()+theme(legend.position = "none", 
                           axis.text = element_text(size = 12),
                           strip.text = element_text(size = 12, face = "bold.italic"),
                           panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                           colour = "#E5E8E8"), 
                           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                           colour = "#E5E8E8"))}


