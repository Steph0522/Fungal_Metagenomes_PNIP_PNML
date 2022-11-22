
library(broom)
library(cowplot)
library(tidyverse)
library(readxl)
library(reshape2)
#Load soil nutrient data
metadata_secas<- read_excel("Data/Metadatos.xlsx", sheet = "secas-marzo")
fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
fq_secas2<- read.csv("Data/fisicoq-la.csv")

meta_fq_secas<- metadata_secas %>% full_join(fq_secas) #%>%  select(Sites:id_new,pH, MO, N, P) %>% mutate(Season="Dry")
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")

env.16S=meta_fq_secas_all
df<- env.16S%>% dplyr::select(SampleID, pH:Mn, moisture, WHC:CONDUC, ARCILLA:ARENA) 
dfs=data.frame(df[1],scale(df[,2:17], center = T, scale = T)) %>% dplyr::select(
  SampleID,P,K,Ca,Mg,moisture,WHC, LIMO) %>% 
 column_to_rownames(var = "SampleID")
dfs_dist<- dist(dfs, method = "euclidean")

coords<- read_csv("Data//coord.csv")

coords_mat<- coords%>% mutate(P=paste0("P",pol),
                              S=paste0("S", Sitio),
                              T=paste0("T", Transecto)
) %>% unite(
  "SampleID", 
  P:T, sep = "") %>% select(
    SampleID, Longitude, Latitude) %>% column_to_rownames(var = "SampleID") %>% as.matrix() 

library(geosphere)

distance<- distm(coords_mat)/1000
colnames(distance)<- rownames(coords_mat)
rownames(distance)<- rownames(coords_mat)
distance_complete<- distance
metadata_secas<- read_excel("Data/Metadatos.xlsx", sheet = "secas-marzo")

fq_secas<- read_excel("Data/fisicoq.xlsx", sheet = "seca")
distance_dm<- read.delim("Data/distance_points.tsv") %>% mutate(value=value/1000)

fq_secas2<- read.csv("Data/fisicoq-la.csv")
meta_fq_secas<- metadata_secas %>% full_join(fq_secas) #%>%  select(Sites:id_new,pH, MO, N, P) %>% mutate(Season="Dry")
meta_fq_secas_all<- metadata_secas %>% full_join(fq_secas) %>% full_join(fq_secas2, by = "SampleID")

ward <- meta_fq_secas_all

#Load mapping file
map <-metadata_secas

#Environmental distance analysis
#Remove variables with no variation, z-transform each variable, 
#and format it as a matrix
#nut.mtx <- ward %>% 
  #select(-LIMO, -Fe, -WHC, -id_sequence:-id_fisicoq, -type, -type2, -pH2, -Names, -Sites) %>% 
 # select(SampleID, Ca, Mg, K,P, moisture,WHC,LIMO) %>% 
  #select(SampleID,pH:ARENA, -type, -type2, -pH2) %>% 
  #  gather(key = "Variable", value = "Value", -SampleID) %>% 
  #group_by(Variable) %>% 
  #mutate(zValue = (Value - mean(Value))/sd(Value)) %>% 
  #select(SampleID, Variable, zValue) %>% 
  #spread(key = Variable, value = zValue) %>% 
  #as.data.frame()
#row.names(nut.mtx) <- nut.mtx$SampleID
#nut.mtx <- nut.mtx[,-1]
#nut.mtx <- as.matrix(nut.mtx)

#Calculate the environmental distance and filter redundant values
#nut.dist2 <- as.matrix(dist(nut.mtx, method = "euclidean"))
nut.dist<- dfs_dist %>% as.matrix()

nut.dist[upper.tri(nut.dist)] <- NA 


#Create a long data frame and remove pairwise comparisons between the same sample
nut.dist.tidy <- nut.dist %>% 
  melt(as.matrix(distance), varnames = c(
    "SampleID.x", "SampleID.y"), value.name = "EucDist") %>% drop_na() %>% filter(
      !EucDist==0) %>% 
  filter(!is.na(EucDist)) %>% 
  filter(SampleID.x != SampleID.y) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID"))  %>% 
  full_join(distance_dm) %>% dplyr::rename(SpatialDistance = "value")


#Perform correlation analysis and regression 
cor <- cor.test(nut.dist.tidy$SpatialDistance, nut.dist.tidy$EucDist, method = "pearson", alternative = "two.sided") %>% tidy()
lm <- lm(EucDist ~ SpatialDistance, data = nut.dist.tidy) %>% tidy() %>% filter(term == "SpatialDistance")
#mant<- vegan::mantel(distance_complete,dist(nut.mtx))
mant2<- vegan::mantel(distance_complete, dfs_dist)
#mant;mant2
dist.stats <- data.frame(label = paste("correlation test: r = ", signif(cor$estimate,2), ",",
                                       "\nlinear regression: slope = ", signif(lm$estimate, 3), ",",
                                       "p-value = ", signif(cor$p.value, 3),
                                       "\nmantel test: r =",  signif(mant2$statistic, 2),
                                       ", p-value = ",  signif(mant2$signif, 3)))






#Plot
environmental.p <- nut.dist.tidy %>% 
  ggplot(aes(SpatialDistance, EucDist)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(color = "black", se = F, method = "lm") +
  scale_color_brewer(name = "Block", palette = "Set1", direction = -1) +
  xlab("Distance between plots (Km)") +
  ylab("Environmental distance") +
 # scale_x_continuous(breaks = seq(0, 18, by = 3)) +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "top")  + theme_linedraw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold.italic"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',   colour = "#E5E8E8"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "#E5E8E8"))+
  ggtitle(dist.stats$label)+theme(plot.title = element_text(size = 14))


environmental.p
ggsave("Figures/Fig.env_distance_geom.png",width = 6, height = 4, dpi = 300, plot = environmental.p, device = "png")
ggsave("Figures/Fig.variables_distance_geom.png",width = 9, height = 7, dpi = 300, plot = variables.p, device = "png")

#Individual variable analysis
#Remove variables with no variation and format it as a matrix
nut.mtx2 <- ward %>% 
  select(SampleID,-P,K,Ca,Mg,moisture,WHC, Silt=LIMO) %>% 
  #select(SampleID,pH:ARENA, -type, -type2, -pH) %>% 
  #select(SampleID,pH:ARENA, -type, -type2, -pH, -pH2, -ARENA, -Cu, -CONDUC) %>% 
    as.data.frame()
row.names(nut.mtx2) <- nut.mtx2$SampleID
nut.mtx2 <- nut.mtx2[,-1]
nut.mtx2 <- as.matrix(nut.mtx2)

#Function to calculate the absolute differences between pairs of samples for each parameter
get_differences <- function(x){
  vec <- nut.mtx2[,x]
  vec %*% t(vec)
  var.dist <- dist(as.matrix(data.frame(x = 0, y = vec))) %>% as.matrix()
  var.dist[upper.tri(var.dist)] <- NA 
  var.dist %>% 
    as.data.frame() %>% 
    mutate(SampleID.x = row.names(.)) %>% 
    gather(key = "SampleID.y", value = "VarDist", -SampleID.x) %>% 
    filter(!is.na(VarDist))
}

#Generate a data frame with the absolute differences and spatial distance for each parameter
var.list <- list()

for(i in 1:ncol(nut.mtx2)){
  var.name <- colnames(nut.mtx2)[i]
  var.list[[var.name]] <- get_differences(i) 
}

var.tidy <- plyr::ldply(var.list, function(x) x) %>% 
  filter(SampleID.x != SampleID.y) %>% 
  inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
  inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
  full_join(distance_dm) %>% dplyr::rename(SpatialDistance = "value")
names(var.tidy)[1] <- "Variable"

#Run correlations and regressions
run_cor <- function(x){
  cor.test(x$SpatialDistance, x$VarDist, method = "pearson", alternative = "two.sided") %>% 
    tidy()
}

#Aggregate all the stats
stats <- var.tidy %>% 
  group_by(Variable) %>% 
  nest() %>% 
  mutate(cor = map(data, run_cor)) %>% 
  unnest(cor) %>% 
  ungroup() %>% 
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>% 
  mutate(label = paste("r = ", signif(estimate,2), ", P = ", signif(p.adj,3), sep = "")) #%>% 
  #inner_join(alt.names, by = "Variable")

#Plot
variables.p <- var.tidy %>% 
  inner_join(stats, by = "Variable") %>% 
  ggplot(aes(SpatialDistance, VarDist)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "gray25") +
  #geom_text(data = bc.stats, aes(x = 9, y = 0.4, label = label), hjust = 0) + 
  geom_smooth(color = "black", se = F, method = "lm") +
  scale_color_brewer(name = "Block", palette = "Set1", direction = -1) +
  xlab("Distance between plots (km)") +
  ylab("Difference between samples") +
  #scale_x_continuous(breaks = seq(0, 18, by = 3)) +
  facet_wrap(~ Variable + label, scales = "free", ncol = 3) + 
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "top",
        strip.background =element_rect(fill="gray25"),
        strip.text = element_text(colour = "white"))+theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "top")  + theme_linedraw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold.italic"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',   colour = "#E5E8E8"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "#E5E8E8"))

variables.p

#Put it all together
top <- plot_grid(NA, environmental.p, NA, rel_widths = c(1,2,1), nrow = 1, labels = c(NA, "a", NA), label_size = 15)

plot_grid(top, variables.p, nrow = 2, rel_heights = c(1,3), labels = c(NA, "b"), label_size = 15)

