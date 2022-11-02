library(tidyverse)
coords<- read_csv("/home/yendi/Documents/corredor_scripts/coord.csv")

coords_mat<- coords%>% mutate(P=paste0("P",pol),
                              S=paste0("S", Sitio),
                              T=paste0("T", Transecto)
                              ) %>% unite(
                                "SampleID", 
                                P:T, sep = "") %>% select(
 SampleID, Longitude, Latitude) %>% column_to_rownames(var = "SampleID") %>% as.matrix() 

library(geosphere)

distance<- distm(coords_mat)
colnames(distance)<- rownames(coords_mat)
rownames(distance)<- rownames(coords_mat)

library(reshape2)
distance[upper.tri(distance)] <- NA 

distance_dm<-melt(as.matrix(distance), varnames = c(
  "SampleID.x", "SampleID.y")) %>% drop_na() %>% filter(!value==0)
ver<- distance_dm %>% arrange(value)

write_tsv(distance_dm, "Data/distance_points.tsv")
