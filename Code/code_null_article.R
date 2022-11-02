#R script for the manuscript: Using null models to compare bacterial and microeukaryotic metacommunity assembly under shifting environmental conditions"
#These script were used for data preparation, data analyses and data visualization. Created in 2018-06-13

#Version: 1.0.136
#R version 3.5.0 


##Loading datasets:
setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")

#OTU tables:
OTU16S=read.csv("16S_OTUfiltered.csv", head=T, row.names = 1, sep=";")
OTU18S=read.csv("18S_OTUfiltered.csv", head=T, row.names = 1, sep=";")

#Separating OTU table from the Taxonomy:
TAX_16S=OTU16S[,183:188]
OTU16S=OTU16S[,1:182]

TAX_18S=OTU18S[,183:193]
OTU18S=OTU18S[,1:182]

#Environmental data:
ENV=read.csv("Env_table.csv", head=T, row.names = 1, sep=";") #this file available on the DiVA repository
p1 <- 'P'
ENV_16S <- subset(ENV, grepl(p1, ENV$ID_16Snew)) #subsetting for 16S
ENV_16S=ENV_16S[,-1]

ENV_18S <- subset(ENV, grepl(p1, ENV$ID_18S)) #subsetting for 18S
ENV_18S=ENV_18S[,-2]

#Phylogeny:
library(ape)
OTU_tree16S <- read.tree("16Stree.tre")
OTU_tree18S <- read.tree("18Stree.tre")

library("phyloseq")
#creating phyloseq pbject for 16S
OTU = otu_table(OTU16S, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(TAX_16S))
ENV_16S=sample_data(ENV_16S)

physeq16S=merge_phyloseq(phyloseq(OTU, TAX), ENV_16S, OTU_tree16S)
physeq16S

#creating phyloseq pbject for 18S
OTU = otu_table(OTU18S, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(TAX_18S))
ENV_18S=sample_data(ENV_18S)

physeq18S=merge_phyloseq(phyloseq(OTU, TAX),ENV_18S, OTU_tree18S)
physeq18S

##Removing temporary pools from all datasets:
physeq16Sss=subset_samples(physeq16S, Code != "15" & Code !="18" & Code !="19" & Code !="23")
physeq16Sss

physeq18Sss=subset_samples(physeq18S, Code != "15" & Code !="18" & Code !="19" & Code !="23")
physeq18Sss

ENV_16Sss=sample_data(physeq16Sss)
ENV_18Sss=sample_data(physeq18Sss)

write.csv(ENV_16Sss, file="/Volumes/FREECOM HDD/Phd_project/R_copy/ENV_16Ssubset.csv")
write.csv(ENV_18Sss, file="/Volumes/FREECOM HDD/Phd_project/R_copy/ENV_18Ssubset.csv")


##Pruning and subsampling:
physeq_p16S=prune_taxa(taxa_sums(physeq16Sss)>10, physeq16Sss)
physeq_p16S
hist(colSums(otu_table(physeq_p16S)))
library(vegan)
rarecurve(t(otu_table(physeq_p16S)), step=50, cex=0.5)
abline(v=4940)

physeq_p18S=prune_taxa(taxa_sums(physeq18Sss)>10, physeq18Sss)
physeq_p18S

rarecurve(t(otu_table(physeq_p18S)), step=50, cex=0.5)
abline(v=3429)

physeq_p16S_rel = rarefy_even_depth(physeq_p16S, sample.size = min(sample_sums(physeq_p16S)))
physeq_p16S_rel
colSums(otu_table(physeq_p16S_rel))

physeq_p18S_rel = rarefy_even_depth(physeq_p18S, sample.size = min(sample_sums(physeq_p18S)))
physeq_p18S_rel
colSums(otu_table(physeq_p18S_rel))

#the following files are available on the DiVA repository
write.csv(otu_table(physeq_p16S_rel), file="/Volumes/FREECOM HDD/Phd_project/R_copy/OTU16S_rar.csv")
write.csv(otu_table(physeq_p18S_rel), file="/Volumes/FREECOM HDD/Phd_project/R_copy/OTU18S_rar.csv")

write.csv(tax_table(physeq_p16S_rel), file="/Volumes/FREECOM HDD/Phd_project/R_copy/OTU16S_tax.csv")
write.csv(tax_table(physeq_p18S_rel), file="/Volumes/FREECOM HDD/Phd_project/R_copy/OTU18S_tax.csv")


#Subsetting fasta.file based on the filtered/prunned, subsampled OTU table:

library("seqinr")
fasta16S<- read.fasta(file = "/Volumes/FREECOM HDD/Phd_project/16S/16Srep.fasta" ,as.string = TRUE, set.attributes = FALSE)

include_list=row.names(otu_table(physeq_p16S_rel))
fasta16Sss=subset(fasta16S, names(fasta16S) %in% include_list)
names(fasta16Sss)

write.fasta(fasta16Sss, names(fasta16Sss), file.out = "/Volumes/FREECOM HDD/Phd_project/16S/16Srep_ss.fasta",nbchar = 60, as.string = FALSE)

fasta18S<- read.fasta(file = "/Volumes/FREECOM HDD/Phd_project/18S/18Srep.fasta" ,as.string = TRUE, set.attributes = FALSE)

include_list=row.names(otu_table(physeq_p18S_rel))
fasta18Sss=subset(fasta18S, names(fasta18S) %in% include_list)
names(fasta18Sss)

write.fasta(fasta18Sss, names(fasta18Sss), file.out = "/Volumes/FREECOM HDD/Phd_project/18S/18Srep_ss.fasta",nbchar = 60, as.string = FALSE)


##################
# Data analyses #
#################


# Assessing environmental conditions and selecting structuring variables
library(reshape2)
library(picante)
env.16S=read.csv("ENV_table.csv", sep=";", row.names = 1)
df <- env.16S[,c("Date","Water_temp","Cond","Depth","Water_col","Chl_a","TN","TP","Daphnia","Copepod","Chironomidae")]
dfs=data.frame(scale(df[,2:11], center = T, scale = T)) #standardize environmental data

dfs$Date<-NA
dfs$Date<- df$Date
melted_df=melt(dfs, id="Date")
mean_table=summarySE(melted_df, measurevar = "value", groupvars=c("Date","variable"), na.rm = T)

# checking collinearity
cor.table(na.omit(dfs[,1:10]))

# RDA separately for 16S and 18S data
library(vegan)
setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")
spp.16S=read.csv("OTU16S_rar.csv", row.names = 1)
env.16S=read.csv("ENV_16Ssubset.csv", sep=";", row.names = 1)

spp.16S=data.frame(t(spp.16S))

mm=cbind(spp.16S,env.16S)
mm=na.omit(mm)
dim(mm)
head(mm)

env.16S=mm[,4588:4607]
spp.16S=mm[,1:4587]

dca=decorana(spp.16S)
plot(dca) # RDA is okay to use

#select variables that did not show collinearity
env.16S=env.16S[,c("Water_temp","Cond","Depth","TN","TP","Daphnia","Copepod","Chironomidae")]

spp.16S_hell=decostand(spp.16S, "hell") # Hellinger transformation
env.16S_st=data.frame(scale(env.16S, scale=T, center=F)) # standardize env. data

PC1=rda(spp.16S_hell~., env.16S_st)
plot(PC1, display = "sites" ,type = c("text"))
plot(PC1, type="n")
text(PC1, dis="cn", lwd=1, cex=1)
text(PC1, dis="sites", type=c("text"), cex=0.8)
points(PC1, dis="sites", lwd=1, cex=0.5)

quartz.save("RDA0_16S.pdf", type="pdf")

# Forward selection procedure
PC1=rda(spp.16S~., env.16S)
plot(PC1, type="n")
text(PC1, dis="cn", lwd=2.5, cex=1.2)
text(PC1, dis="sites", type=c("text"), cex=0.8)
points(PC1, dis="sites", lwd=2.5)

cap.env=capscale(spp.16S_hell~.,env.16S_st,distance = "bray")
cap.env

mod0.env=capscale(spp.16S_hell~1,env.16S_st,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
step.env$anova # These are the selected structuring variables that are included in the subsequent analyses for 16S


# Let's repeat the procedure but with the microeukaryotic dataset
setwd("/Volumes/FREECOM HDD/Phd_project/R_copy/")
spp.18S=read.csv("OTU18S_rar.csv", row.names = 1)
env.18S=read.csv("ENV_18Ssubset.csv", sep=",", row.names = 1)
spp.18S=data.frame(t(spp.18S))

mm=cbind(spp.18S,env.18S)
mm=na.omit(mm)
dim(mm)
head(mm)

env.18S=mm[,1337:1356]
spp.18S=mm[,1:1336]

dca=decorana(spp.18S)
plot(dca) # RDA is okay to use

#select variables that did not show collinearity
env.18S=env.18S[,c("Water_temp","Cond","Depth","TN","TP","Daphnia","Copepod","Chironomidae")]

env.18S_st=data.frame(scale(env.18S, scale = T, center = F))
spp.18S_hell=decostand(spp.18S, "hell")

PC1=rda(spp.18S_hell~., env.18S_st)
plot(PC1, display = "sites" ,type = c("text"))
plot(PC1, type="n")
text(PC1, dis="cn", lwd=1, cex=1)
text(PC1, dis="sites", type=c("text"), cex=0.8)
points(PC1, dis="sites", lwd=1, cex=0.5)

quartz.save("RDA0_18S.pdf", type="pdf")

# Forward selection procedure
PC1=rda(spp.18S~., env.18S)
plot(PC1, type="n")
text(PC1, dis="cn", lwd=2.5, cex=1.2)
text(PC1, dis="sites", type=c("text"), cex=0.8)
points(PC1, dis="sites", lwd=2.5)

cap.env=capscale(spp.18S_hell~.,env.18S_st,distance = "bray")
cap.env

mod0.env=capscale(spp.18S_hell~1,env.18S_st,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
step.env$anova # These are the selected structuring variables that are included in the subsequent analyses for 18S

# Plot the temporal trends of the selected environmental variables over time
# Figure 2
setwd("/Volumes/FREECOM HDD/Phd_project/R_copy/")
env.16S=read.csv("ENV_table.csv", sep=";", row.names = 1)
df <- env.16S[,c("Date","Water_temp","Cond","Depth","TN","TP","Daphnia","Copepod")]
melted_df=melt(df, id="Date")
mean=summarySE(melted_df, measurevar = "value", groupvars=c("Date","variable"), na.rm = T)

mean$variable=factor(mean$variable, levels = c("Water_temp","Depth","Cond","TP","TN","Daphnia","Copepod"))

conservation_status <- c(
  Water_temp = "Water temperature (°C)",
  Depth = "Depth (cm)",
  Cond = "Conductivity (µS/cm)",
  TP = "Total phosphorus (µg/L)",
  TN = "Total nitrogen (µg/L)",
  Daphnia = "Daphnia abundance (ind/L)",
  Copepod = "Copepoda abundance (ind/L)")
library(plyr)
mean$conservation2 <- plyr::revalue(mean$variable, conservation_status)

library(ggplot2)
break.vec <- c(seq(from = as.Date("2015-08-14"), to = as.Date("2015-09-19"),
                   by = "4 day"))

meanplot=ggplot(data=mean, aes(x=as.Date(Date), y=value, group=variable), size=1.1) + geom_point()+
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_line()+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1)+
  facet_grid(conservation2~., scales = "free_y", labeller = labeller(conservation2=label_wrap_gen(15)))
meanplot=meanplot+theme_bw(base_size = 10)+theme_linedraw()
meanplot=meanplot+theme_linedraw(base_size=10)+labs(x="Sampling occasion",y="Mean±SD")+theme(legend.position="none")
meanplot
scaleFUN <- function(x) sprintf("%.0f", x)
meanplot=meanplot + 
  scale_y_continuous(labels=scaleFUN)+
  theme(axis.text = element_text(size=8),
        axis.title.x = element_text(size=10, colour = "black", vjust = -0.4),panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.title = element_text(size=8, face="bold"),axis.text.x = element_text(angle = 90, vjust = 0.8),
        panel.grid = element_blank(),
        strip.text.y = element_text(margin = margin(0,0,0,0, "cm")))+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")
meanplot+labs(subtitle = "    Dry period        Wet period")+theme(plot.subtitle = element_text(face="bold"))

quartz.save("meanplot.pdf", type = "pdf") 
quartz.save("meanplot.tiff", type = "tiff")

# It is clear from the Figure 2 plot that there are two distinct periods. But are they statistically different?
#### non-parametric test between dry and wet periods
# Table S1 in Supplementary material
library(car)
kruskal.test(env.16S$Cond~env.16S$period)
leveneTest(env.16S$Cond~env.16S$period)

kruskal.test(env.16S$Copepod~env.16S$period)
leveneTest(env.16S$Copepod~env.16S$period)

kruskal.test(env.16S$Daphnia~env.16S$period)
leveneTest(env.16S$Daphnia~env.16S$period)

kruskal.test(env.16S$Depth~env.16S$period)
leveneTest(env.16S$Depth~env.16S$period)

kruskal.test(env.16S$TN~env.16S$period)
leveneTest(env.16S$TN~env.16S$period)

kruskal.test(env.16S$TP~env.16S$period)
leveneTest(env.16S$TP~env.16S$period)

kruskal.test(env.16S$Water_temp~env.16S$period)
leveneTest(env.16S$Water_temp~env.16S$period)


#NMDS; assessing dry-wet periods 
library(ggfortify)
library(cluster)
library(vegan)

# NMDS plot showing the changes in environmental conditions in relation to the dry
# and wet period based on Euclidean distances

setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")
env.16S=read.csv("ENV_table.csv", sep=";", row.names = 1)

env.16S=env.16S[which(env.16S$Code=="10"),]

df <- env.16S[,c("Water_temp","Cond","Depth","TN","TP","Copepod","Daphnia","period")]
df=na.omit(df)

dff=df[,c("Cond","Water_temp","TN","TP","Depth","Daphnia","Copepod")]
dff=scale(dff, scale = T, center = F)
dff <- dff[,colSums(is.na(dff))<nrow(dff)]

dffd = (vegdist(dff, "euc"))
anova(betadisper(dffd, df$period)) 
adonis_location = adonis(dff ~ df$period, method = "euc") #PERMANOVA
adonis_location$aov.tab

dff = as.matrix((vegdist(dff, "euc")))
vare.mds=metaMDS(dff, distance = "bray", autotransform = F, na.rm=T, k=3)
data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$period <- df$period  #  add the grp variable created earlier
data.scores  #look at the data


grp.a <- data.scores[data.scores$period == "dry", ][chull(data.scores[data.scores$period == 
                                                                   "dry", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$period == "wet", ][chull(data.scores[data.scores$period == 
                                                                   "wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp B


hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

# Figure S7 in the Supplementary material
P20=ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=period, fill=period),alpha=0.2) +
  scale_fill_manual(values=c("dry" = "#CC4800", "wet" = "#3078CB")) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=period),size=2) + # add the point markers
  scale_shape_manual(values=c(19,21))+
  theme_bw() + 
  theme(axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=10), # remove x-axis labels
        axis.title.y = element_text(size=10), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "right")

P20 + annotate("text", x = 5, y = -4.5, label = "stress: 0.056")
quartz.save("ENV_DIS.pdf", type="pdf") 


# Separation of environmental conditions between the dry and wet period in each rock pool based on 
# PERMANOVA and PERMDISP using Euclidean distances for the environmental data
# Table S2 and S3 in Supplementary material

setwd("/Volumes/FREECOM HDD/Phd_project/")
env.16S=read.csv("ENV_table.csv", sep=";", row.names = 1)

#subsetting to individual pool; run this line only 
#if you want to check the separation at individual rock pool level
env.16S=env.16S[which(env.16S$Code=="28"),] #change rock pool ID
# otherwise, continue from this line:
df <- env.16S[,c("Water_temp","Cond","Depth","TN","TP","Copepod","Daphnia","period")]
df=na.omit(df)

dff=df[,c("Cond","Water_temp","Depth","TN","TP","Copepod","Daphnia")]
dff=scale(dff, scale = T, center = F)
dff <- dff[,colSums(is.na(dff))<nrow(dff)]

adonis_location = adonis(dff ~ df$period, method = "euc")
adonis_location$aov.tab #PERMANOVA
dffd = (vegdist(dff, "euc"))
anova(betadisper(dffd, df$period)) #PERMDISP


# NMDS plot showing the changes in bacterial or microeukaryotic communities in relation to the dry and wet
# period based on Bray-Curtis dissimilarities 
setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")
otu.16S=read.csv("OTU18S_rar.csv", sep=",", row.names = 1) # Change to "OTU16S_rar.csv" if you want to assess the bacterial dataset
otu.16S=data.frame(t(otu.16S))
otulist=row.names(otu.16S)

otu.16S_hell=decostand(otu.16S, "hell")

#subsetting to individual pool; run this lines only 
#if you want to check the separation at individual rock pool level
df=subset(env.16S, row.names(env.16S) %in% otulist)
df=df[which(df$Code=="28"),] # change rock pool ID
sitelist=row.names(df)
#otherwise run from this line:
otu.16S_hell=subset(otu.16S_hell, row.names(otu.16S_hell) %in% sitelist)


vare.mds=metaMDS(otu.16S_hell, distance = "bray", autotransform = F, k=3)

data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$period <- df$period  #  add the grp variable created earlier
data.scores$code <- df$Code
data.scores$date <- df$Date
data.scores  #look at the data

grp.a <- data.scores[data.scores$period == "dry", ][chull(data.scores[data.scores$period == 
                                                                        "dry", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$period == "wet", ][chull(data.scores[data.scores$period == 
                                                                        "wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

# Figure S8 or S9 depending on the OTU data you loaded
P20=ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=period, fill=period),alpha=0.2) +
  scale_fill_manual(values=c("dry" = "#CC4800", "wet" = "#3078CB")) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=period),size=2) + # add the point markers
  scale_shape_manual(values=c(19,21))+
  theme_bw() + 
  theme(axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=8), # remove x-axis labels
        axis.title.y = element_text(size=8), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "right")

P20+ annotate("text", x = -1.3, y = -1.5, label = "stress: 0.144")

quartz.save("OTU18S_DIS.pdf", type="pdf")

# PERMANOVA and PERMDISP using Bray-Curtis distances for the bacterial or microeukaryotic dataset
# Table S2 and S3 in Supplementary material

adonis_location = adonis(otu.16S_hell ~ df$period, method = "bray")
adonis_location$aov.tab #PERMANOVA
dffd = (vegdist(otu.16S_hell, "bray"))
anova(betadisper(dffd, df$period)) #PERMDISP



############################
## Null model approaches ##
##########################

# Elements of Metacommunity Structure (EMS)

# EMS for 16S
abund_table<-read.csv("OTU16S_rar.csv",row.names=1,check.names=FALSE)
metadata16S<-read.csv("ENV_16Ssubset.csv",header=T,row.names=1, sep=",")

meta_table=metadata16S
meta_table$Date
tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)

abund_table=t(abund_table)
row.names(abund_table)

abund_table814=subset(abund_table, row.names(abund_table) %in% tp814i)
abund_table818=subset(abund_table, row.names(abund_table) %in% tp818i)
abund_table822=subset(abund_table, row.names(abund_table) %in% tp822i)
abund_table826=subset(abund_table, row.names(abund_table) %in% tp826i)
abund_table830=subset(abund_table, row.names(abund_table) %in% tp830i)
abund_table903=subset(abund_table, row.names(abund_table) %in% tp903i)
abund_table907=subset(abund_table, row.names(abund_table) %in% tp907i)
abund_table911=subset(abund_table, row.names(abund_table) %in% tp911i)
abund_table915=subset(abund_table, row.names(abund_table) %in% tp915i)
abund_table919=subset(abund_table, row.names(abund_table) %in% tp919i)

abund_table814 = abund_table814[,which(colSums(abund_table814) != 0)]
abund_table818 = abund_table818[,which(colSums(abund_table818) != 0)]
abund_table822 = abund_table822[,which(colSums(abund_table822) != 0)]
abund_table826 = abund_table826[,which(colSums(abund_table826) != 0)]
abund_table830 = abund_table830[,which(colSums(abund_table830) != 0)]
abund_table903 = abund_table903[,which(colSums(abund_table903) != 0)]
abund_table907 = abund_table907[,which(colSums(abund_table907) != 0)]
abund_table911 = abund_table911[,which(colSums(abund_table911) != 0)]
abund_table915 = abund_table915[,which(colSums(abund_table915) != 0)]
abund_table919 = abund_table919[,which(colSums(abund_table919) != 0)]

library(metacom)

#814
Met814=Metacommunity(abund_table814, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met814$Coherence)
tur1r1=data.frame(Met814$Turnover)
clum1r1=t(Met814$Boundary)
write.csv(coh1r1, "./Coh1_814_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_814_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_814_r1.csv",quote=F)

OM814=OrderMatrix(abund_table814, outputScores=T, binary=F)
write.csv(data.frame(OM814$sitescores), "./OM814site.csv", quote=F)

#818
Met818=Metacommunity(abund_table818, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met818$Coherence)
tur1r1=data.frame(Met818$Turnover)
clum1r1=t(Met818$Boundary)
write.csv(coh1r1, "./Coh1_818_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_818_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_818_r1.csv",quote=F)

OM818=OrderMatrix(abund_table818, outputScores=T, binary=F)
write.csv(data.frame(OM818$sitescores), "./OM818site.csv", quote=F)

#822
Met822=Metacommunity(abund_table822, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met822$Coherence)
tur1r1=data.frame(Met822$Turnover)
clum1r1=t(Met822$Boundary)
write.csv(coh1r1, "./Coh1_822_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_822_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_822_r1.csv",quote=F)

OM822=OrderMatrix(abund_table822, outputScores=T, binary=F)
write.csv(data.frame(OM822$sitescores), "./OM822site.csv", quote=F)

#826
Met826=Metacommunity(abund_table826, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met826$Coherence)
tur1r1=data.frame(Met826$Turnover)
clum1r1=t(Met826$Boundary)
write.csv(coh1r1, "./Coh1_826_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_826_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_826_r1.csv",quote=F)

OM826=OrderMatrix(abund_table826, outputScores=T, binary=F)
write.csv(data.frame(OM826$sitescores), "./OM826site.csv", quote=F)

#830
Met830=Metacommunity(abund_table830, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met830$Coherence)
tur1r1=data.frame(Met830$Turnover)
clum1r1=t(Met830$Boundary)
write.csv(coh1r1, "./Coh1_830_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_830_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_830_r1.csv",quote=F)

OM830=OrderMatrix(abund_table830, outputScores=T, binary=F)
write.csv(data.frame(OM830$sitescores), "./OM830site.csv", quote=F)

#903
Met903=Metacommunity(abund_table903, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met903$Coherence)
tur1r1=data.frame(Met903$Turnover)
clum1r1=t(Met903$Boundary)
write.csv(coh1r1, "./Coh1_903_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_903_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_903_r1.csv",quote=F)

OM903=OrderMatrix(abund_table903, outputScores=T, binary=F)
write.csv(data.frame(OM903$sitescores), "./OM903site.csv", quote=F)

#907
Met907=Metacommunity(abund_table907, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met907$Coherence)
tur1r1=data.frame(Met907$Turnover)
clum1r1=t(Met907$Boundary)
write.csv(coh1r1, "./Coh1_907_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_907_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_907_r1.csv",quote=F)

OM907=OrderMatrix(abund_table907, outputScores=T, binary=F)
write.csv(data.frame(OM907$sitescores), "./OM907site.csv", quote=F)

#911
Met911=Metacommunity(abund_table911, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met911$Coherence)
tur1r1=data.frame(Met911$Turnover)
clum1r1=t(Met911$Boundary)
write.csv(coh1r1, "./Coh1_911_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_911_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_911_r1.csv",quote=F)

OM911=OrderMatrix(abund_table911, outputScores=T, binary=F)
write.csv(data.frame(OM911$sitescores), "./OM911site.csv", quote=F)

#915
Met915=Metacommunity(abund_table915, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met915$Coherence)
tur1r1=data.frame(Met915$Turnover)
clum1r1=t(Met915$Boundary)
write.csv(coh1r1, "./Coh1_915_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_915_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_915_r1.csv",quote=F)

OM915=OrderMatrix(abund_table915, outputScores=T, binary=F)
write.csv(data.frame(OM915$sitescores), "./OM915site.csv", quote=F)

#919
Met919=Metacommunity(abund_table919, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met919$Coherence)
tur1r1=data.frame(Met919$Turnover)
clum1r1=t(Met919$Boundary)
write.csv(coh1r1, "./Coh1_919_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_919_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_919_r1.csv",quote=F)

OM919=OrderMatrix(abund_table919, outputScores=T, binary=F)
write.csv(data.frame(OM919$sitescores), "./OM919site.csv", quote=F)

# EMS for 18S
abund_table<-read.csv("OTU18S_rar.csv",row.names=1,check.names=FALSE)
metadata16S<-read.csv("ENV_18Ssubset.csv",header=T,row.names=1, sep=",")

meta_table=metadata16S
meta_table$Date
tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)

abund_table=t(abund_table)
row.names(abund_table)

abund_table814=subset(abund_table, row.names(abund_table) %in% tp814i)
abund_table818=subset(abund_table, row.names(abund_table) %in% tp818i)
abund_table822=subset(abund_table, row.names(abund_table) %in% tp822i)
abund_table826=subset(abund_table, row.names(abund_table) %in% tp826i)
abund_table830=subset(abund_table, row.names(abund_table) %in% tp830i)
abund_table903=subset(abund_table, row.names(abund_table) %in% tp903i)
abund_table907=subset(abund_table, row.names(abund_table) %in% tp907i)
abund_table911=subset(abund_table, row.names(abund_table) %in% tp911i)
abund_table915=subset(abund_table, row.names(abund_table) %in% tp915i)
abund_table919=subset(abund_table, row.names(abund_table) %in% tp919i)

abund_table814 = abund_table814[,which(colSums(abund_table814) != 0)]
abund_table818 = abund_table818[,which(colSums(abund_table818) != 0)]
abund_table822 = abund_table822[,which(colSums(abund_table822) != 0)]
abund_table826 = abund_table826[,which(colSums(abund_table826) != 0)]
abund_table830 = abund_table830[,which(colSums(abund_table830) != 0)]
abund_table903 = abund_table903[,which(colSums(abund_table903) != 0)]
abund_table907 = abund_table907[,which(colSums(abund_table907) != 0)]
abund_table911 = abund_table911[,which(colSums(abund_table911) != 0)]
abund_table915 = abund_table915[,which(colSums(abund_table915) != 0)]
abund_table919 = abund_table919[,which(colSums(abund_table919) != 0)]

library(metacom)

#814
Met814=Metacommunity(abund_table814, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met814$Coherence)
tur1r1=data.frame(Met814$Turnover)
clum1r1=t(Met814$Boundary)
write.csv(coh1r1, "./Coh1_814_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_814_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_814_r1.csv",quote=F)

OM814=OrderMatrix(abund_table814, outputScores=T, binary=F)
write.csv(data.frame(OM814$sitescores), "./OM814site.csv", quote=F)

#818
Met818=Metacommunity(abund_table818, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met818$Coherence)
tur1r1=data.frame(Met818$Turnover)
clum1r1=t(Met818$Boundary)
write.csv(coh1r1, "./Coh1_818_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_818_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_818_r1.csv",quote=F)

OM818=OrderMatrix(abund_table818, outputScores=T, binary=F)
write.csv(data.frame(OM818$sitescores), "./OM818site.csv", quote=F)

#822
Met822=Metacommunity(abund_table822, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met822$Coherence)
tur1r1=data.frame(Met822$Turnover)
clum1r1=t(Met822$Boundary)
write.csv(coh1r1, "./Coh1_822_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_822_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_822_r1.csv",quote=F)

OM822=OrderMatrix(abund_table822, outputScores=T, binary=F)
write.csv(data.frame(OM822$sitescores), "./OM822site.csv", quote=F)

#826
Met826=Metacommunity(abund_table826, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met826$Coherence)
tur1r1=data.frame(Met826$Turnover)
clum1r1=t(Met826$Boundary)
write.csv(coh1r1, "./Coh1_826_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_826_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_826_r1.csv",quote=F)

OM826=OrderMatrix(abund_table826, outputScores=T, binary=F)
write.csv(data.frame(OM826$sitescores), "./OM826site.csv", quote=F)

#830
Met830=Metacommunity(abund_table830, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met830$Coherence)
tur1r1=data.frame(Met830$Turnover)
clum1r1=t(Met830$Boundary)
write.csv(coh1r1, "./Coh1_830_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_830_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_830_r1.csv",quote=F)

OM830=OrderMatrix(abund_table830, outputScores=T, binary=F)
write.csv(data.frame(OM830$sitescores), "./OM830site.csv", quote=F)

#903
Met903=Metacommunity(abund_table903, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met903$Coherence)
tur1r1=data.frame(Met903$Turnover)
clum1r1=t(Met903$Boundary)
write.csv(coh1r1, "./Coh1_903_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_903_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_903_r1.csv",quote=F)

OM903=OrderMatrix(abund_table903, outputScores=T, binary=F)
write.csv(data.frame(OM903$sitescores), "./OM903site.csv", quote=F)

#907
Met907=Metacommunity(abund_table907, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met907$Coherence)
tur1r1=data.frame(Met907$Turnover)
clum1r1=t(Met907$Boundary)
write.csv(coh1r1, "./Coh1_907_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_907_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_907_r1.csv",quote=F)

OM907=OrderMatrix(abund_table907, outputScores=T, binary=F)
write.csv(data.frame(OM907$sitescores), "./OM907site.csv", quote=F)

#911
Met911=Metacommunity(abund_table911, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met911$Coherence)
tur1r1=data.frame(Met911$Turnover)
clum1r1=t(Met911$Boundary)
write.csv(coh1r1, "./Coh1_911_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_911_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_911_r1.csv",quote=F)

OM911=OrderMatrix(abund_table911, outputScores=T, binary=F)
write.csv(data.frame(OM911$sitescores), "./OM911site.csv", quote=F)

#915
Met915=Metacommunity(abund_table915, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met915$Coherence)
tur1r1=data.frame(Met915$Turnover)
clum1r1=t(Met915$Boundary)
write.csv(coh1r1, "./Coh1_915_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_915_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_915_r1.csv",quote=F)

OM915=OrderMatrix(abund_table915, outputScores=T, binary=F)
write.csv(data.frame(OM915$sitescores), "./OM915site.csv", quote=F)

#919
Met919=Metacommunity(abund_table919, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
coh1r1=data.frame(Met919$Coherence)
tur1r1=data.frame(Met919$Turnover)
clum1r1=t(Met919$Boundary)
write.csv(coh1r1, "./Coh1_919_r1.csv",quote=F)
write.csv(tur1r1, "./Tur1_919_r1.csv",quote=F)
write.csv(clum1r1, "./Clum1_919_r1.csv",quote=F)

OM919=OrderMatrix(abund_table919, outputScores=T, binary=F)
write.csv(data.frame(OM919$sitescores), "./OM919site.csv", quote=F)


# Thereafter I compiled the outcomes of the EMS tests for both 16S and 18S in Excel and load it below for visualization
# loading for the plots

setwd("/Volumes/FREECOM HDD/Phd_project/EMS/")
EMS=read.csv("R_input2.csv", sep=";")

# Testing differences in coherence (z-value) between periods
library(car)
EMSb=EMS[which(EMS$group=="b"),]
kruskal.test(EMSb$coh_z~EMSb$period)
EMSe=EMS[which(EMS$group=="e"),]
kruskal.test(EMSe$coh_z~EMSe$period)

# Plot coherence (z-values)
# Figure 3
scaleFUN <- function(x) sprintf("%.3f", x)
conservation_status <- c(b = "Bacteria",e = "Microeukaryotes")
sp10 <- ggplot(EMS, aes(x=as.Date(Date), y=coh_z)) +
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_hline(yintercept=1.96, linetype=3, alpha=0.8)+
  geom_hline(yintercept=-1.96, linetype=3, alpha=0.8)+
  geom_point(aes(shape=factor(pattern),color=turn_z, size=clump))+
  scale_size(name="Boundary clumping\n(Morisita's index)", breaks = c(1, 1.1, 1.2, 1.3))+
  scale_shape_manual(values=c(15,17,16), labels=c("Checkerboards", "Nested-\nclumped species loss", "Random"), name="Metacommunity type")+
  scale_colour_gradient(low = "black", high = "grey", name="Turnover (z-value)")
sp10=sp10+ facet_wrap(~group, strip.position="right", ncol=1,labeller=labeller(group=conservation_status))+
  theme_bw(base_size = 15)+theme_linedraw()+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")
sp10=sp10+theme_linedraw()
plot10=sp10+theme_bw(base_size = 15)+theme_linedraw(base_size=15)+labs(x="Sampling occasion", y="Coherence (z-value)")
plot10+ theme(plot.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.text = element_text(size=13),
        strip.text = element_text(hjust = 0.5, size=15))+
  guides(shape = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(size = 5)))+
  labs(subtitle = "    Dry period           Wet period")+theme(plot.subtitle = element_text(face="bold"))

quartz.save("EMS_complett_MS.pdf", type="pdf")
quartz.save("EMS_complett_MS.tiff", type="tiff")


# Incidence-based Raup-Crick beta-diversity

setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")
abund_table<-read.csv("OTU18S_rar.csv",row.names=1)
metadata18S<-read.csv("ENV_18Ssubset.csv",header=T,row.names=1, sep=",")

meta_table=metadata18S
meta_table$Date
tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)

abund_table=t(abund_table)
row.names(abund_table)

abund_table814=subset(abund_table, row.names(abund_table) %in% tp814i)
abund_table818=subset(abund_table, row.names(abund_table) %in% tp818i)
abund_table822=subset(abund_table, row.names(abund_table) %in% tp822i)
abund_table826=subset(abund_table, row.names(abund_table) %in% tp826i)
abund_table830=subset(abund_table, row.names(abund_table) %in% tp830i)
abund_table903=subset(abund_table, row.names(abund_table) %in% tp903i)
abund_table907=subset(abund_table, row.names(abund_table) %in% tp907i)
abund_table911=subset(abund_table, row.names(abund_table) %in% tp911i)
abund_table915=subset(abund_table, row.names(abund_table) %in% tp915i)
abund_table919=subset(abund_table, row.names(abund_table) %in% tp919i)

abund_table814 = abund_table814[,which(colSums(abund_table814) != 0)]
abund_table818 = abund_table818[,which(colSums(abund_table818) != 0)]
abund_table822 = abund_table822[,which(colSums(abund_table822) != 0)]
abund_table826 = abund_table826[,which(colSums(abund_table826) != 0)]
abund_table830 = abund_table830[,which(colSums(abund_table830) != 0)]
abund_table903 = abund_table903[,which(colSums(abund_table903) != 0)]
abund_table907 = abund_table907[,which(colSums(abund_table907) != 0)]
abund_table911 = abund_table911[,which(colSums(abund_table911) != 0)]
abund_table915 = abund_table915[,which(colSums(abund_table915) != 0)]
abund_table919 = abund_table919[,which(colSums(abund_table919) != 0)]

rc1=raup_crick(abund_table814, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc2=raup_crick(abund_table818, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3=raup_crick(abund_table822, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4=raup_crick(abund_table826, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5=raup_crick(abund_table830, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6=raup_crick(abund_table903, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc7=raup_crick(abund_table907, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc8=raup_crick(abund_table911, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc9=raup_crick(abund_table915, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc10=raup_crick(abund_table919, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

rcc1=subset(melt(as.matrix(rc1)),value!=0)
rcc2=subset(melt(as.matrix(rc2)),value!=0)
rcc3=subset(melt(as.matrix(rc3)),value!=0)
rcc4=subset(melt(as.matrix(rc4)),value!=0)
rcc5=subset(melt(as.matrix(rc5)),value!=0)
rcc6=subset(melt(as.matrix(rc6)),value!=0)
rcc7=subset(melt(as.matrix(rc7)),value!=0)
rcc8=subset(melt(as.matrix(rc8)),value!=0)
rcc9=subset(melt(as.matrix(rc9)),value!=0)
rcc10=subset(melt(as.matrix(rc10)),value!=0)

rcc1$date=1
rcc2$date=2
rcc3$date=3
rcc4$date=4
rcc5$date=5
rcc6$date=6
rcc7$date=7
rcc8$date=8
rcc9$date=9
rcc10$date=10

rcc_m=rbind(rcc1, rcc2,rcc3,rcc4,rcc5,rcc6,rcc7,rcc8,rcc9,rcc10)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summary_se=summarySE(rcc_m, measurevar = "value", groupvars = "date")
# Let's save the summary_se object
write.csv(summary_se, "./RC_bact.csv")


# same procedure but with the 18S dataset
abund_table<-read.csv("OTU18S_rar.csv",row.names=1)
metadata18S<-read.csv("ENV_18Ssubset.csv",header=T,row.names=1, sep=",")

meta_table=metadata18S
meta_table$Date
tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)

abund_table=t(abund_table)
row.names(abund_table)

abund_table814=subset(abund_table, row.names(abund_table) %in% tp814i)
abund_table818=subset(abund_table, row.names(abund_table) %in% tp818i)
abund_table822=subset(abund_table, row.names(abund_table) %in% tp822i)
abund_table826=subset(abund_table, row.names(abund_table) %in% tp826i)
abund_table830=subset(abund_table, row.names(abund_table) %in% tp830i)
abund_table903=subset(abund_table, row.names(abund_table) %in% tp903i)
abund_table907=subset(abund_table, row.names(abund_table) %in% tp907i)
abund_table911=subset(abund_table, row.names(abund_table) %in% tp911i)
abund_table915=subset(abund_table, row.names(abund_table) %in% tp915i)
abund_table919=subset(abund_table, row.names(abund_table) %in% tp919i)

abund_table814 = abund_table814[,which(colSums(abund_table814) != 0)]
abund_table818 = abund_table818[,which(colSums(abund_table818) != 0)]
abund_table822 = abund_table822[,which(colSums(abund_table822) != 0)]
abund_table826 = abund_table826[,which(colSums(abund_table826) != 0)]
abund_table830 = abund_table830[,which(colSums(abund_table830) != 0)]
abund_table903 = abund_table903[,which(colSums(abund_table903) != 0)]
abund_table907 = abund_table907[,which(colSums(abund_table907) != 0)]
abund_table911 = abund_table911[,which(colSums(abund_table911) != 0)]
abund_table915 = abund_table915[,which(colSums(abund_table915) != 0)]
abund_table919 = abund_table919[,which(colSums(abund_table919) != 0)]

rc1=raup_crick(abund_table814, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc2=raup_crick(abund_table818, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3=raup_crick(abund_table822, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4=raup_crick(abund_table826, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5=raup_crick(abund_table830, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6=raup_crick(abund_table903, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc7=raup_crick(abund_table907, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc8=raup_crick(abund_table911, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc9=raup_crick(abund_table915, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc10=raup_crick(abund_table919, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

rcc1=subset(melt(as.matrix(rc1)),value!=0)
rcc2=subset(melt(as.matrix(rc2)),value!=0)
rcc3=subset(melt(as.matrix(rc3)),value!=0)
rcc4=subset(melt(as.matrix(rc4)),value!=0)
rcc5=subset(melt(as.matrix(rc5)),value!=0)
rcc6=subset(melt(as.matrix(rc6)),value!=0)
rcc7=subset(melt(as.matrix(rc7)),value!=0)
rcc8=subset(melt(as.matrix(rc8)),value!=0)
rcc9=subset(melt(as.matrix(rc9)),value!=0)
rcc10=subset(melt(as.matrix(rc10)),value!=0)

rcc1$date=1
rcc2$date=2
rcc3$date=3
rcc4$date=4
rcc5$date=5
rcc6$date=6
rcc7$date=7
rcc8$date=8
rcc9$date=9
rcc10$date=10

rcc_m=rbind(rcc1, rcc2,rcc3,rcc4,rcc5,rcc6,rcc7,rcc8,rcc9,rcc10)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summary_se=summarySE(rcc_m, measurevar = "value", groupvars = "date")
# Let's save the summary_se object
write.csv(summary_se, "./RC_euk.csv")


# After the calculation of the incidence-based Raup-Crick beta-diversity, I compiled the results together for visualization
# Load the compiled results
setwd("/Volumes/FREECOM HDD/Phd_project/Betadiv/RC/")
betadiv=read.csv("RC_result.csv", sep=";")

# Testing the differences between periods
library(car)
betab=betadiv[which(betadiv$group=="Bacteria"),]
kruskal.test(betab$value~betab$period)
betae=betadiv[which(betadiv$group=="Microeukaryotes"),]
kruskal.test(betae$value~betae$period)

# Plot incidence-based Raup-Crick beta-diversity
# Figure 4
as.numeric(betadiv$date)
bd16=ggplot(betadiv, aes(x=as.Date(date), y=value, group=group, lty=group)) +
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, lty=1)+
  geom_point(aes(shape=group),size=5)+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")
bd16=bd16+facet_grid(group~.)+theme_bw(base_size = 18)+theme_linedraw()
plot16=bd16+theme_linedraw(base_size=18)+labs(x="Sampling occasion",y="Incidence-based beta-diversity (±SE)")+theme(legend.position="none")
plot16
scaleFUN <- function(x) sprintf("%.2f", x)
plot16=plot16 + 
  scale_y_continuous(labels=scaleFUN)+
  theme(axis.text = element_text(size=16),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size=16, angle = 90),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid = element_blank(),
        strip.text = element_text(hjust = 0.5, size=18))+
  labs(subtitle = "   Dry period       Wet period")+theme(plot.subtitle = element_text(face="bold"))
plot16

quartz.save("Betadiv_RC.pdf", type="pdf")
quartz.save("Betadiv_RC.tiff", type="tiff")



# Qantitative Process Estimate (QPE) from Stegen et al. 2013
# This approach requires that phylogenetic distances (PD) among taxa reflect differences in 
# the ecological niches they inhabit, thus, carry a phylogenetic signal. 
# The presence of phylogenetic signals was tested using Mantel correlograms, as described in Stegen et al. 2013
# The necessary R script was provided by Jianjun Wang and used in the study of Langenheder et al. 2017 FEMS Microbiology Ecology

# After performing phylogenetic signal tests, let's use mantel.correlog function and plot the results
# The following script was run separately on 16S and 18S dataset!
library(readr)
library(vegan)
library(picante)
library(ape)
library(dplyr)

Cond_spp.dist.4.correlog <- read.csv("Cond_spp.dist.4.correlog.csv", row.names=1,header=T)
Cond_phylo.dist.4.correlog <- read.csv("Cond_phylo.dist.4.correlog.csv", row.names=1,header=T)

phylo.sig.correlog = mantel.correlog(Cond_spp.dist.4.correlog,Cond_phylo.dist.4.correlog,
                                     nperm=999,
                                     cutoff=FALSE,
                                     n.class=50,
                                     mult="bonferroni");
phylo.sig.correlog[["mantel.res"]]

write.csv(phylo.sig.correlog$mantel.res, "Cond_signal.csv",  quote=F);


TP_spp.dist.4.correlog <- read.csv("TP_spp.dist.4.correlog.csv", row.names=1,header=T)
TP_phylo.dist.4.correlog <- read.csv("TP_phylo.dist.4.correlog.csv", row.names=1,header=T)

phylo.sig.correlog = mantel.correlog(TP_spp.dist.4.correlog,TP_phylo.dist.4.correlog,
                                     nperm=999,
                                     cutoff=FALSE,
                                     n.class=50,
                                     mult="bonferroni");
phylo.sig.correlog[["mantel.res"]]

write.csv(phylo.sig.correlog$mantel.res, "TP_signal.csv",  quote=F);



TN_spp.dist.4.correlog <- read.csv("TN_spp.dist.4.correlog.csv", row.names=1,header=T)
TN_phylo.dist.4.correlog <- read.csv("TN_phylo.dist.4.correlog.csv", row.names=1,header=T)

phylo.sig.correlog = mantel.correlog(TN_spp.dist.4.correlog,TN_phylo.dist.4.correlog,
                                     nperm=999,
                                     cutoff=FALSE,
                                     n.class=50,
                                     mult="bonferroni");
phylo.sig.correlog[["mantel.res"]]

write.csv(phylo.sig.correlog$mantel.res, "TN_signal.csv",  quote=F);

Daphnia_spp.dist.4.correlog <- read.csv("Daphnia_spp.dist.4.correlog.csv", row.names=1,header=T)
Daphnia_phylo.dist.4.correlog <- read.csv("Daphnia_phylo.dist.4.correlog.csv", row.names=1,header=T)

phylo.sig.correlog = mantel.correlog(Daphnia_spp.dist.4.correlog,Daphnia_phylo.dist.4.correlog,
                                     nperm=999,
                                     cutoff=FALSE,
                                     n.class=50,
                                     mult="bonferroni");
phylo.sig.correlog[["mantel.res"]]

write.csv(phylo.sig.correlog$mantel.res, "Daphnia_signal.csv",  quote=F);


Depth_spp.dist.4.correlog <- read.csv("Depth_spp.dist.4.correlog.csv", row.names=1,header=T)
Depth_phylo.dist.4.correlog <- read.csv("Depth_phylo.dist.4.correlog.csv", row.names=1,header=T)

phylo.sig.correlog = mantel.correlog(Depth_spp.dist.4.correlog,Depth_phylo.dist.4.correlog,
                                     nperm=999,
                                     cutoff=FALSE,
                                     n.class=50,
                                     mult="bonferroni");
phylo.sig.correlog[["mantel.res"]]

write.csv(phylo.sig.correlog$mantel.res, "Depth_signal.csv",  quote=F);


Water_temp_spp.dist.4.correlog <- read.csv("Water_temp_spp.dist.4.correlog.csv", row.names=1,header=T)
Water_temp_phylo.dist.4.correlog <- read.csv("Water_temp_phylo.dist.4.correlog.csv", row.names=1,header=T)

phylo.sig.correlog = mantel.correlog(Water_temp_spp.dist.4.correlog,Water_temp_phylo.dist.4.correlog,
                                     nperm=999,
                                     cutoff=FALSE,
                                     n.class=50,
                                     mult="bonferroni");
phylo.sig.correlog[["mantel.res"]]

write.csv(phylo.sig.correlog$mantel.res, "Water_temp_signal.csv",  quote=F);



# Plotting the Mantel correlograms
# Figure S5 and S6

# 16S dataset (Figure S5)
setwd("/Volumes/FREECOM HDD/Phd_project/Phylo_signal/16S/results")

par(mfrow=c(3,2))
phylo.sig<-read.csv("Cond_signal.csv", row.names = 1, header=T, sep=",")
p16Cond=plot(phylo.sig[,c(1,3)],
             xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Conductivity")

phylo.sig<-read.csv("TP_signal.csv", row.names = 1, header=T, sep=",")
p16TP=plot(phylo.sig[,c(1,3)],
           xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Total phosphorus")

phylo.sig<-read.csv("TN_signal.csv", row.names = 1, header=T, sep=",")
p16TN=plot(phylo.sig[,c(1,3)],
           xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Total nitrogen")

phylo.sig<-read.csv("Water_temp_signal.csv", row.names = 1, header=T, sep=",")
p16Water_temp=plot(phylo.sig[,c(1,3)],
                   xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Water temperature")

phylo.sig<-read.csv("Depth_signal.csv", row.names = 1, header=T, sep=",")
p16Depth=plot(phylo.sig[,c(1,3)],
              xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Depth")

phylo.sig<-read.csv("Daphnia_signal.csv", row.names = 1, header=T, sep=",")
p16Daphnia=plot(phylo.sig[,c(1,3)],
                xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Daphnia")

quartz.save("Phylo_signal16S.pdf", type="pdf")


# 18S dataset (Figure S6)
setwd("/Volumes/FREECOM HDD/Phd_project/Phylo_signal/18S/results")

dev.off()
par(mfrow=c(3,3))
phylo.sig<-read.csv("Cond_signal.csv", row.names = 1, header=T, sep=",")
p18Cond=plot(phylo.sig[,c(1,3)],
             xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Conductivity")

phylo.sig<-read.csv("TP_signal.csv", row.names = 1, header=T, sep=",")
p18TP=plot(phylo.sig[,c(1,3)],
           xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Total phosphorus")

phylo.sig<-read.csv("Water_temp_signal.csv", row.names = 1, header=T, sep=",")
p18WT=plot(phylo.sig[,c(1,3)],
           xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Water temperature")

phylo.sig<-read.csv("Daphnia_signal.csv", row.names = 1, header=T, sep=",")
p18Daphnia=plot(phylo.sig[,c(1,3)],
                xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Daphnia")

phylo.sig<-read.csv("TN_signal.csv", row.names = 1, header=T, sep=",")
p18TN=plot(phylo.sig[,c(1,3)],
           xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Total nitrogen")

phylo.sig<-read.csv("Depth_signal.csv", row.names = 1, header=T, sep=",")
p18Depth=plot(phylo.sig[,c(1,3)],
              xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Depth")

phylo.sig<-read.csv("Copepod_signal.csv", row.names = 1, header=T, sep=",")
p18Copepod=plot(phylo.sig[,c(1,3)],
                xlab="Phylogenetic Distance Class", ylab="Mantel Test Statistic") +lines(phylo.sig[,c(1,3)]) + points(phylo.sig[,c(1,3)], pch=21, bg="white", col="black", pty=4)+
  points(phylo.sig[,c(1,3)][phylo.sig[,5] < 0.05,], pch=21, bg="black", col="black", pty=4) + abline(h=0, lty=2, col="red") +
  title("Copepoda")

quartz.save("Phylo_signal18S.pdf", type="pdf")


# After I evaluated the existence of phylogenetic signal, QPE analyses can be performed
# First step of the QPE approach (it can be time-consuming process, so I run it on UPPMAX server)

# For 16S dataset
abund_table<-read.csv("OTU16S_rar.csv",row.names=1,check.names=FALSE)
row.names(abund_table) <- paste("'", row.names(abund_table),"'", sep = "")

library(ape)
OTU_tree <- read.tree("16Stree.tre")

metadata16S<-read.csv("ENV_16Ssubset.csv",header=T,row.names=1, sep=",")

meta_table=metadata16S
meta_table$Date
tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)


library(picante)
m=match.phylo.data(OTU_tree, abund_table)
OTU_tree=m$phy
abund_table=m$data

abund_table=t(abund_table)
row.names(abund_table)

tp814otu=t(subset(abund_table, row.names(abund_table) %in% tp814i))
tp818otu=t(subset(abund_table, row.names(abund_table) %in% tp818i))
tp822otu=t(subset(abund_table, row.names(abund_table) %in% tp822i))
tp826otu=t(subset(abund_table, row.names(abund_table) %in% tp826i))
tp830otu=t(subset(abund_table, row.names(abund_table) %in% tp830i))
tp903otu=t(subset(abund_table, row.names(abund_table) %in% tp903i))
tp907otu=t(subset(abund_table, row.names(abund_table) %in% tp907i))
tp911otu=t(subset(abund_table, row.names(abund_table) %in% tp911i))
tp915otu=t(subset(abund_table, row.names(abund_table) %in% tp915i))
tp919otu=t(subset(abund_table, row.names(abund_table) %in% tp919i))

m=match.phylo.data(OTU_tree, tp814otu)
abund_table814=data.frame(m$data)
abund_table814 = abund_table814[which(rowSums(abund_table814) != 0),]
m=match.phylo.data(OTU_tree, abund_table814)
OTU_tree814=m$phy

m=match.phylo.data(OTU_tree, tp818otu)
abund_table818=data.frame(m$data)
abund_table818 = abund_table818[which(rowSums(abund_table818) != 0),]
m=match.phylo.data(OTU_tree, abund_table818)
OTU_tree818=m$phy

m=match.phylo.data(OTU_tree, tp822otu)
abund_table822=data.frame(m$data)
abund_table822 = abund_table822[which(rowSums(abund_table822) != 0),]
m=match.phylo.data(OTU_tree, abund_table822)
OTU_tree822=m$phy

m=match.phylo.data(OTU_tree, tp826otu)
abund_table826=data.frame(m$data)
abund_table826 = abund_table826[which(rowSums(abund_table826) != 0),]
m=match.phylo.data(OTU_tree, abund_table826)
OTU_tree826=m$phy

m=match.phylo.data(OTU_tree, tp830otu)
abund_table830=data.frame(m$data)
abund_table830 = abund_table830[which(rowSums(abund_table830) != 0),]
m=match.phylo.data(OTU_tree, abund_table830)
OTU_tree830=m$phy

m=match.phylo.data(OTU_tree, tp903otu)
abund_table903=data.frame(m$data)
abund_table903 = abund_table903[which(rowSums(abund_table903) != 0),]
m=match.phylo.data(OTU_tree, abund_table903)
OTU_tree903=m$phy

m=match.phylo.data(OTU_tree, tp907otu)
abund_table907=data.frame(m$data)
abund_table907 = abund_table907[which(rowSums(abund_table907) != 0),]
m=match.phylo.data(OTU_tree, abund_table907)
OTU_tree907=m$phy

m=match.phylo.data(OTU_tree, tp911otu)
abund_table911=data.frame(m$data)
abund_table911 = abund_table911[which(rowSums(abund_table911) != 0),]
m=match.phylo.data(OTU_tree, abund_table911)
OTU_tree911=m$phy

m=match.phylo.data(OTU_tree, tp915otu)
abund_table915=data.frame(m$data)
abund_table915 = abund_table915[which(rowSums(abund_table915) != 0),]
m=match.phylo.data(OTU_tree, abund_table915)
OTU_tree915=m$phy

m=match.phylo.data(OTU_tree, tp919otu)
abund_table919=data.frame(m$data)
abund_table919 = abund_table919[which(rowSums(abund_table919) != 0),]
m=match.phylo.data(OTU_tree, abund_table919)
OTU_tree919=m$phy

#814
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table814),cophenetic(OTU_tree814),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table814),colnames(beta.mntd.weighted));
identical(colnames(abund_table814),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table814),ncol(abund_table814),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table814),taxaShuffle(cophenetic(OTU_tree814)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table814),ncol=ncol(abund_table814));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table814)-1)) {
  for (rows in (columns+1):ncol(abund_table814)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table814);
colnames(weighted.bNTI) = colnames(abund_table814);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI814.csv",quote=F)

#818
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table818),cophenetic(OTU_tree818),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table818),colnames(beta.mntd.weighted));
identical(colnames(abund_table818),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table818),ncol(abund_table818),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table818),taxaShuffle(cophenetic(OTU_tree818)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table818),ncol=ncol(abund_table818));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table818)-1)) {
  for (rows in (columns+1):ncol(abund_table818)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table818);
colnames(weighted.bNTI) = colnames(abund_table818);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI818.csv",quote=F)

#822
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table822),cophenetic(OTU_tree822),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table822),colnames(beta.mntd.weighted));
identical(colnames(abund_table822),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table822),ncol(abund_table822),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table822),taxaShuffle(cophenetic(OTU_tree822)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table822),ncol=ncol(abund_table822));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table822)-1)) {
  for (rows in (columns+1):ncol(abund_table822)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table822);
colnames(weighted.bNTI) = colnames(abund_table822);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI822.csv",quote=F)

#826
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table826),cophenetic(OTU_tree826),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table826),colnames(beta.mntd.weighted));
identical(colnames(abund_table826),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table826),ncol(abund_table826),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table826),taxaShuffle(cophenetic(OTU_tree826)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table826),ncol=ncol(abund_table826));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table826)-1)) {
  for (rows in (columns+1):ncol(abund_table826)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table826);
colnames(weighted.bNTI) = colnames(abund_table826);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI826.csv",quote=F)

#830
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table830),cophenetic(OTU_tree830),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table830),colnames(beta.mntd.weighted));
identical(colnames(abund_table830),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table830),ncol(abund_table830),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table830),taxaShuffle(cophenetic(OTU_tree830)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table830),ncol=ncol(abund_table830));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table830)-1)) {
  for (rows in (columns+1):ncol(abund_table830)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table830);
colnames(weighted.bNTI) = colnames(abund_table830);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI830.csv",quote=F)

#903
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table903),cophenetic(OTU_tree903),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table903),colnames(beta.mntd.weighted));
identical(colnames(abund_table903),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table903),ncol(abund_table903),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table903),taxaShuffle(cophenetic(OTU_tree903)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table903),ncol=ncol(abund_table903));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table903)-1)) {
  for (rows in (columns+1):ncol(abund_table903)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table903);
colnames(weighted.bNTI) = colnames(abund_table903);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI903.csv",quote=F)

#907
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table907),cophenetic(OTU_tree907),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table907),colnames(beta.mntd.weighted));
identical(colnames(abund_table907),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table907),ncol(abund_table907),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table907),taxaShuffle(cophenetic(OTU_tree907)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table907),ncol=ncol(abund_table907));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table907)-1)) {
  for (rows in (columns+1):ncol(abund_table907)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table907);
colnames(weighted.bNTI) = colnames(abund_table907);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI907.csv",quote=F)

#911
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table911),cophenetic(OTU_tree911),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table911),colnames(beta.mntd.weighted));
identical(colnames(abund_table911),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table911),ncol(abund_table911),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table911),taxaShuffle(cophenetic(OTU_tree911)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table911),ncol=ncol(abund_table911));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table911)-1)) {
  for (rows in (columns+1):ncol(abund_table911)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table911);
colnames(weighted.bNTI) = colnames(abund_table911);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI911.csv",quote=F)

#915
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table915),cophenetic(OTU_tree915),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table915),colnames(beta.mntd.weighted));
identical(colnames(abund_table915),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table915),ncol(abund_table915),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table915),taxaShuffle(cophenetic(OTU_tree915)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table915),ncol=ncol(abund_table915));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table915)-1)) {
  for (rows in (columns+1):ncol(abund_table915)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table915);
colnames(weighted.bNTI) = colnames(abund_table915);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI915.csv",quote=F)

#919
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table919),cophenetic(OTU_tree919),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table919),colnames(beta.mntd.weighted));
identical(colnames(abund_table919),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table919),ncol(abund_table919),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table919),taxaShuffle(cophenetic(OTU_tree919)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table919),ncol=ncol(abund_table919));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table919)-1)) {
  for (rows in (columns+1):ncol(abund_table919)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table919);
colnames(weighted.bNTI) = colnames(abund_table919);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI919.csv",quote=F)

# For 18S dataset
abund_table<-read.csv("OTU18S_rar.csv",row.names=1,check.names=FALSE)
row.names(abund_table) <- paste("'", row.names(abund_table),"'", sep = "")

library(ape)
OTU_tree <- read.tree("18Stree.tre")

metadata16S<-read.csv("ENV_18Ssubset.csv",header=T,row.names=1, sep=",")

meta_table=metadata16S
meta_table$Date
tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)


library(picante)
m=match.phylo.data(OTU_tree, abund_table)
OTU_tree=m$phy
abund_table=m$data

abund_table=t(abund_table)
row.names(abund_table)

tp814otu=t(subset(abund_table, row.names(abund_table) %in% tp814i))
tp818otu=t(subset(abund_table, row.names(abund_table) %in% tp818i))
tp822otu=t(subset(abund_table, row.names(abund_table) %in% tp822i))
tp826otu=t(subset(abund_table, row.names(abund_table) %in% tp826i))
tp830otu=t(subset(abund_table, row.names(abund_table) %in% tp830i))
tp903otu=t(subset(abund_table, row.names(abund_table) %in% tp903i))
tp907otu=t(subset(abund_table, row.names(abund_table) %in% tp907i))
tp911otu=t(subset(abund_table, row.names(abund_table) %in% tp911i))
tp915otu=t(subset(abund_table, row.names(abund_table) %in% tp915i))
tp919otu=t(subset(abund_table, row.names(abund_table) %in% tp919i))

m=match.phylo.data(OTU_tree, tp814otu)
abund_table814=data.frame(m$data)
abund_table814 = abund_table814[which(rowSums(abund_table814) != 0),]
m=match.phylo.data(OTU_tree, abund_table814)
OTU_tree814=m$phy

m=match.phylo.data(OTU_tree, tp818otu)
abund_table818=data.frame(m$data)
abund_table818 = abund_table818[which(rowSums(abund_table818) != 0),]
m=match.phylo.data(OTU_tree, abund_table818)
OTU_tree818=m$phy

m=match.phylo.data(OTU_tree, tp822otu)
abund_table822=data.frame(m$data)
abund_table822 = abund_table822[which(rowSums(abund_table822) != 0),]
m=match.phylo.data(OTU_tree, abund_table822)
OTU_tree822=m$phy

m=match.phylo.data(OTU_tree, tp826otu)
abund_table826=data.frame(m$data)
abund_table826 = abund_table826[which(rowSums(abund_table826) != 0),]
m=match.phylo.data(OTU_tree, abund_table826)
OTU_tree826=m$phy

m=match.phylo.data(OTU_tree, tp830otu)
abund_table830=data.frame(m$data)
abund_table830 = abund_table830[which(rowSums(abund_table830) != 0),]
m=match.phylo.data(OTU_tree, abund_table830)
OTU_tree830=m$phy

m=match.phylo.data(OTU_tree, tp903otu)
abund_table903=data.frame(m$data)
abund_table903 = abund_table903[which(rowSums(abund_table903) != 0),]
m=match.phylo.data(OTU_tree, abund_table903)
OTU_tree903=m$phy

m=match.phylo.data(OTU_tree, tp907otu)
abund_table907=data.frame(m$data)
abund_table907 = abund_table907[which(rowSums(abund_table907) != 0),]
m=match.phylo.data(OTU_tree, abund_table907)
OTU_tree907=m$phy

m=match.phylo.data(OTU_tree, tp911otu)
abund_table911=data.frame(m$data)
abund_table911 = abund_table911[which(rowSums(abund_table911) != 0),]
m=match.phylo.data(OTU_tree, abund_table911)
OTU_tree911=m$phy

m=match.phylo.data(OTU_tree, tp915otu)
abund_table915=data.frame(m$data)
abund_table915 = abund_table915[which(rowSums(abund_table915) != 0),]
m=match.phylo.data(OTU_tree, abund_table915)
OTU_tree915=m$phy

m=match.phylo.data(OTU_tree, tp919otu)
abund_table919=data.frame(m$data)
abund_table919 = abund_table919[which(rowSums(abund_table919) != 0),]
m=match.phylo.data(OTU_tree, abund_table919)
OTU_tree919=m$phy

#814
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table814),cophenetic(OTU_tree814),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table814),colnames(beta.mntd.weighted));
identical(colnames(abund_table814),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table814),ncol(abund_table814),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table814),taxaShuffle(cophenetic(OTU_tree814)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table814),ncol=ncol(abund_table814));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table814)-1)) {
  for (rows in (columns+1):ncol(abund_table814)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table814);
colnames(weighted.bNTI) = colnames(abund_table814);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI814.csv",quote=F)

#818
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table818),cophenetic(OTU_tree818),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table818),colnames(beta.mntd.weighted));
identical(colnames(abund_table818),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table818),ncol(abund_table818),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table818),taxaShuffle(cophenetic(OTU_tree818)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table818),ncol=ncol(abund_table818));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table818)-1)) {
  for (rows in (columns+1):ncol(abund_table818)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table818);
colnames(weighted.bNTI) = colnames(abund_table818);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI818.csv",quote=F)

#822
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table822),cophenetic(OTU_tree822),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table822),colnames(beta.mntd.weighted));
identical(colnames(abund_table822),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table822),ncol(abund_table822),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table822),taxaShuffle(cophenetic(OTU_tree822)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table822),ncol=ncol(abund_table822));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table822)-1)) {
  for (rows in (columns+1):ncol(abund_table822)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table822);
colnames(weighted.bNTI) = colnames(abund_table822);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI822.csv",quote=F)

#826
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table826),cophenetic(OTU_tree826),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table826),colnames(beta.mntd.weighted));
identical(colnames(abund_table826),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table826),ncol(abund_table826),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table826),taxaShuffle(cophenetic(OTU_tree826)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table826),ncol=ncol(abund_table826));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table826)-1)) {
  for (rows in (columns+1):ncol(abund_table826)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table826);
colnames(weighted.bNTI) = colnames(abund_table826);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI826.csv",quote=F)

#830
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table830),cophenetic(OTU_tree830),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table830),colnames(beta.mntd.weighted));
identical(colnames(abund_table830),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table830),ncol(abund_table830),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table830),taxaShuffle(cophenetic(OTU_tree830)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table830),ncol=ncol(abund_table830));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table830)-1)) {
  for (rows in (columns+1):ncol(abund_table830)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table830);
colnames(weighted.bNTI) = colnames(abund_table830);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI830.csv",quote=F)

#903
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table903),cophenetic(OTU_tree903),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table903),colnames(beta.mntd.weighted));
identical(colnames(abund_table903),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table903),ncol(abund_table903),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table903),taxaShuffle(cophenetic(OTU_tree903)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table903),ncol=ncol(abund_table903));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table903)-1)) {
  for (rows in (columns+1):ncol(abund_table903)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table903);
colnames(weighted.bNTI) = colnames(abund_table903);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI903.csv",quote=F)

#907
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table907),cophenetic(OTU_tree907),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table907),colnames(beta.mntd.weighted));
identical(colnames(abund_table907),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table907),ncol(abund_table907),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table907),taxaShuffle(cophenetic(OTU_tree907)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table907),ncol=ncol(abund_table907));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table907)-1)) {
  for (rows in (columns+1):ncol(abund_table907)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table907);
colnames(weighted.bNTI) = colnames(abund_table907);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI907.csv",quote=F)

#911
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table911),cophenetic(OTU_tree911),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table911),colnames(beta.mntd.weighted));
identical(colnames(abund_table911),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table911),ncol(abund_table911),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table911),taxaShuffle(cophenetic(OTU_tree911)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table911),ncol=ncol(abund_table911));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table911)-1)) {
  for (rows in (columns+1):ncol(abund_table911)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table911);
colnames(weighted.bNTI) = colnames(abund_table911);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI911.csv",quote=F)

#915
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table915),cophenetic(OTU_tree915),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table915),colnames(beta.mntd.weighted));
identical(colnames(abund_table915),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table915),ncol(abund_table915),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table915),taxaShuffle(cophenetic(OTU_tree915)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table915),ncol=ncol(abund_table915));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table915)-1)) {
  for (rows in (columns+1):ncol(abund_table915)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table915);
colnames(weighted.bNTI) = colnames(abund_table915);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI915.csv",quote=F)

#919
beta.mntd.weighted = as.matrix(comdistnt(t(abund_table919),cophenetic(OTU_tree919),abundance.weighted=T));
dim(beta.mntd.weighted);

identical(colnames(abund_table919),colnames(beta.mntd.weighted));
identical(colnames(abund_table919),rownames(beta.mntd.weighted));

beta.reps = 999; 
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table919),ncol(abund_table919),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table919),taxaShuffle(cophenetic(OTU_tree919)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table919),ncol=ncol(abund_table919));
dim(weighted.bNTI);

for (columns in 1:(ncol(abund_table919)-1)) {
  for (rows in (columns+1):ncol(abund_table919)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(abund_table919);
colnames(weighted.bNTI) = colnames(abund_table919);
weighted.bNTI;
write.csv(weighted.bNTI, "./weighted_bNTI919.csv",quote=F)

# Second step of the QPE approach (abundance-based Raup-Crick beta-diversity)
raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.
  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##make_null:
  
  ##looping over each pairwise community combination:
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_bray_curtis<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');
        
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); # null.spXsite;
        
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.spXsite,method='bray-curtis');
        
      }; # end reps loop
      
      ## empirically observed bray curtis
      obs.bray = distance(spXsite[c(null.one,null.two),],method='bray-curtis');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      
      print(c(null.one,null.two,date()));
      
    }; ## end null.two loop
    
  }; ## end null.one loop
  
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }	
  
  return(results)
  
}; ## end function

# For 16S dataset
comm<-read.csv("OTU16S_rar.csv",row.names=1,check.names=FALSE)
row.names(comm) <- paste("'", row.names(comm),"'", sep = "")

library(ape)
Tree <- read.tree("16Stree.tre")

library(picante)
mc=match.phylo.data(Tree, comm)

comm=mc$data
comm=data.frame(t(comm))

metadata16S<-read.csv("ENV_16Ssubset.csv",header=T,row.names=1, sep=";")
meta_table=metadata16S
meta_table$Date

tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]

tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)

tp814otu=subset(comm, row.names(comm) %in% tp814i)
tp818otu=subset(comm, row.names(comm) %in% tp818i)
tp822otu=subset(comm, row.names(comm) %in% tp822i)
tp826otu=subset(comm, row.names(comm) %in% tp826i)
tp830otu=subset(comm, row.names(comm) %in% tp830i)
tp903otu=subset(comm, row.names(comm) %in% tp903i)
tp907otu=subset(comm, row.names(comm) %in% tp907i)
tp911otu=subset(comm, row.names(comm) %in% tp911i)
tp915otu=subset(comm, row.names(comm) %in% tp915i)
tp919otu=subset(comm, row.names(comm) %in% tp919i)

abund_table814 = data.frame(tp814otu[,which(colSums(tp814otu) != 0)])
abund_table818 = data.frame(tp818otu[,which(colSums(tp818otu) != 0)])
abund_table822 = data.frame(tp822otu[,which(colSums(tp822otu) != 0)])
abund_table826 = data.frame(tp826otu[,which(colSums(tp826otu) != 0)])
abund_table830 = data.frame(tp830otu[,which(colSums(tp830otu) != 0)])
abund_table903 = data.frame(tp903otu[,which(colSums(tp903otu) != 0)])
abund_table907 = data.frame(tp907otu[,which(colSums(tp907otu) != 0)])
abund_table911 = data.frame(tp911otu[,which(colSums(tp911otu) != 0)])
abund_table915 = data.frame(tp915otu[,which(colSums(tp915otu) != 0)])
abund_table919 = data.frame(tp919otu[,which(colSums(tp919otu) != 0)])


results=raup_crick_abundance(abund_table814, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results814.csv")


results=raup_crick_abundance(abund_table818, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results818.csv")

results=raup_crick_abundance(abund_table822, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results822.csv")

results=raup_crick_abundance(abund_table826, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results826.csv")

results=raup_crick_abundance(abund_table830, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results830.csv")

results=raup_crick_abundance(abund_table903, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results903.csv")

results=raup_crick_abundance(abund_table907, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results907.csv")

results=raup_crick_abundance(abund_table911, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results911.csv")

results=raup_crick_abundance(abund_table915, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results915.csv")

results=raup_crick_abundance(abund_table919, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results919.csv")


# For 18S dataset
comm<-read.csv("OTU18S_rar.csv",row.names=1,check.names=FALSE)
row.names(comm) <- paste("'", row.names(comm),"'", sep = "")



library(ape)
Tree <- read.tree("18Stree.tre")


library(picante)
mc=match.phylo.data(Tree, comm)

comm=mc$data
comm=data.frame(t(comm))

metadata16S<-read.csv("ENV_18Ssubset.csv",header=T,row.names=1, sep=",")
meta_table=metadata16S
meta_table$Date

tp814=meta_table[which(meta_table$Date=="2015-08-14"), ]
tp818=meta_table[which(meta_table$Date=="2015-08-18"), ]
tp822=meta_table[which(meta_table$Date=="2015-08-22"), ]
tp826=meta_table[which(meta_table$Date=="2015-08-26"), ]
tp830=meta_table[which(meta_table$Date=="2015-08-30"), ]
tp903=meta_table[which(meta_table$Date=="2015-09-03"), ]
tp907=meta_table[which(meta_table$Date=="2015-09-07"), ]
tp911=meta_table[which(meta_table$Date=="2015-09-11"), ]
tp915=meta_table[which(meta_table$Date=="2015-09-15"), ]
tp919=meta_table[which(meta_table$Date=="2015-09-19"), ]


tp814i=row.names(tp814)
tp818i=row.names(tp818)
tp822i=row.names(tp822)
tp826i=row.names(tp826)
tp830i=row.names(tp830)
tp903i=row.names(tp903)
tp907i=row.names(tp907)
tp911i=row.names(tp911)
tp915i=row.names(tp915)
tp919i=row.names(tp919)


tp814otu=subset(comm, row.names(comm) %in% tp814i)
tp818otu=subset(comm, row.names(comm) %in% tp818i)
tp822otu=subset(comm, row.names(comm) %in% tp822i)
tp826otu=subset(comm, row.names(comm) %in% tp826i)
tp830otu=subset(comm, row.names(comm) %in% tp830i)
tp903otu=subset(comm, row.names(comm) %in% tp903i)
tp907otu=subset(comm, row.names(comm) %in% tp907i)
tp911otu=subset(comm, row.names(comm) %in% tp911i)
tp915otu=subset(comm, row.names(comm) %in% tp915i)
tp919otu=subset(comm, row.names(comm) %in% tp919i)

abund_table814 = data.frame(tp814otu[,which(colSums(tp814otu) != 0)])
abund_table818 = data.frame(tp818otu[,which(colSums(tp818otu) != 0)])
abund_table822 = data.frame(tp822otu[,which(colSums(tp822otu) != 0)])
abund_table826 = data.frame(tp826otu[,which(colSums(tp826otu) != 0)])
abund_table830 = data.frame(tp830otu[,which(colSums(tp830otu) != 0)])
abund_table903 = data.frame(tp903otu[,which(colSums(tp903otu) != 0)])
abund_table907 = data.frame(tp907otu[,which(colSums(tp907otu) != 0)])
abund_table911 = data.frame(tp911otu[,which(colSums(tp911otu) != 0)])
abund_table915 = data.frame(tp915otu[,which(colSums(tp915otu) != 0)])
abund_table919 = data.frame(tp919otu[,which(colSums(tp919otu) != 0)])


results=raup_crick(abund_table814, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results814.csv")


results=raup_crick(abund_table818, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results818.csv")

results=raup_crick(abund_table822, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results822.csv")

results=raup_crick(abund_table826, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results826.csv")

results=raup_crick(abund_table830, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results830.csv")

results=raup_crick(abund_table903, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results903.csv")

results=raup_crick(abund_table907, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results907.csv")

results=raup_crick(abund_table911, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results911.csv")

results=raup_crick(abund_table915, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results915.csv")

results=raup_crick(abund_table919, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
write.csv(as.matrix(results), "./results919.csv")

# End of the two-step procedure of QPE approach

# Summarized the outcomes and calculated manually the relative contribution of each assembly process

# QPE results
# Testing differences between periods by Kruskal-Wallis test
# I created an extra vector (in Excel) for the results to link QPE results to the two periods (dry vs. wet)
S16S=read.csv("period.csv", sep=";")

kruskal.test(S16S$Drift~S16S$period) # Here, change x ("Drift") if you want to test other process
hist(S16S$Drift)
boxplot(S16S$Drift~S16S$period)

# Plotting the results of QPE

library(reshape2)
library(ggplot2)
library(scales)

setwd("/Volumes/FREECOM HDD/Phd_project/Assemblymechanisms/16S/")
S16S=read.csv("results_16S_abund.csv", sep=";")
S16S=melt(S16S, id="date")

setwd("/Volumes/FREECOM HDD/Phd_project/Assemblymechanisms/18S/")
S18S=read.csv("results_18S_abund.csv", sep=";")
S18S=melt(S18S, id="date")

break.vec <- c(seq(from = as.Date("2015-08-14"), to = as.Date("2015-09-19"),
                   by = "4 day"))

#fill <- c("#334f3b","#65816d","#c8cbbf","#bd7423","#072f37")
fill <- c("#dfc27d","#a6611a","#80cdc1","#018571","black") 

conservation_status=c(Variable_selection="Variable selection", Homogeneous_selection="Homogeneous selection", 
                      Homogenizing_dispersal="Homogenizing dispersal",
                      Dispersal_limitation="Dispersal limitation or\nhistorical contingencies",Drift="Drift")


f_labs <- c(`dry` = "Dry \nperiod",
            `wet` = "Wet \nperiod")
p_title <- "Overall\nprocess estimates"
p_titleeuk <- "Overall\nprocess estimates"
positions <- c("dry","wet")

S16=ggplot() +
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_bar(data=S16S,aes(y=value, x=as.Date(date),fill=variable), stat = "identity")+
  facet_wrap(~variable, ncol=1, scales="free_y",labeller = labeller(variable=conservation_status))+
  scale_fill_manual(values = fill, labels=conservation_status)
S16=S16+labs(x="Study period", y="Proportion (%)")
S16=S16+theme_linedraw(base_size=12)+theme(legend.position = "bottom", 
                                           legend.direction = "horizontal", legend.title = element_blank(), 
                                           legend.text = element_text(size=12))+
  theme(axis.title = element_text(size=12),axis.text = element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))
S16=S16+theme(plot.title = element_text(size=12, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), vjust = 0),
        panel.grid = element_line(colour = "lightgray"),
        strip.text = element_text(size=12,margin = margin(0,0,0,0, "cm")),
        plot.margin = unit(c(0.2,0.5,0.5,0.5), "cm"),
        panel.spacing.y = unit(0.1, "lines"))+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")+
  guides(fill=guide_legend(nrow=2))+
  labs(subtitle = "     Dry period      Wet period", tag="C")+theme(plot.subtitle = element_text(face="bold"),
                                                                    plot.tag = element_text(face = "bold", size = 14))
S16


S18=ggplot() +
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_bar(data=S18S,aes(y=value, x=as.Date(date),fill=variable), stat = "identity")+
  facet_wrap(~variable, ncol=1, scales = "free_y", labeller = labeller(variable=conservation_status))+
  scale_fill_manual(values = fill, labels=conservation_status)
S18=S18+labs(x="Study period", y="Proportion (%)")+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")
S18=S18+theme_linedraw(base_size=12)+theme(legend.position = "bottom", legend.direction = "horizontal", 
                                           legend.title = element_blank(), legend.text = element_text(size=12))+
  theme(axis.title = element_text(size=12),axis.text = element_text(size=12),panel.border = element_rect(colour = "black", fill=NA, size=0.8))
S18=S18+
  theme(plot.title = element_text(size=12, face="bold"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), vjust = 0),
        axis.title.y = element_text(colour = "white"),
        panel.grid = element_line(colour = "lightgray"),
        strip.text.x = element_text(size=12,margin = margin(0,0,0,0, "cm")),
        plot.margin = unit(c(0.2,0.5,0.5,0.5), "cm"),
        panel.spacing.y = unit(0.1, "lines"))+
  guides(fill=guide_legend(nrow=2))+
  labs(subtitle = "     Dry period      Wet period", tag="D")+theme(plot.subtitle = element_text(face="bold"),
                                                                    plot.tag = element_text(face = "bold", size=14))

# Stacked barplot
setwd("/Volumes/FREECOM HDD/Phd_project/Assemblymechanisms/16S/")
pie16=read.csv("piechart16S.csv", sep=";", header = T)
setwd("/Volumes/FREECOM HDD/Phd_project/Assemblymechanisms/18S/")
pie18=read.csv("piechart18S.csv", sep=";", header = T)
pie16m=melt(pie16, id="period")
pie18m=melt(pie18, id="period")
p_title <- "A"
p_titleeuk <- "B"

p <- ggplot(pie16m, aes(y = value, x=period, fill = variable))
p=p + geom_bar(stat = "identity", color = "white", lwd=0.2) +
  scale_x_discrete(labels = as_labeller(f_labs), limits=positions) +
  scale_fill_manual(values = fill, labels=conservation_status)+
  guides(fill = guide_legend(override.aes = list(width = 1),
                             reverse = T,
                             title.position = "top.left",
                             label.position = "top",
                             nrow = 3)) +
  labs(x = NULL, y = "Proportion (%)",
       fill = NULL,
       tag = p_title) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        axis.ticks.length = unit(0, "cm"),
        axis.title = element_text(size=12),
        axis.title.y = element_text(margin = margin(t=0, r=0, b=0, l=0)),
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size=12),
        plot.title = element_text(size=14, face="bold"),
        plot.margin = unit(c(1,2,0,2), "cm"),
        plot.tag = element_text(face = "bold", size=14))


peuk <- ggplot(pie18m, aes(y = value, x=period, fill = variable))
peuk=peuk + geom_bar(stat = "identity", color = "white", lwd=0.2) +
  scale_x_discrete(labels = as_labeller(f_labs), limits=positions) +
  scale_fill_manual(values = fill, labels=conservation_status)+
  guides(fill = guide_legend(override.aes = list(width = 1),
                             reverse = T,
                             title.position = "top.left",
                             label.position = "top",
                             nrow = 3)) +
  labs(x = NULL, y = "Proportion (%)",
       fill = NULL,
       tag = p_titleeuk) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size=12),
        axis.title.y = element_text(margin = margin(t=0, r=0, b=0, l=0)),
        axis.ticks.length = unit(0, "cm"),
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size=12),
        plot.title = element_text(size=12),
        plot.margin = unit(c(1,2,0,2), "cm"),
        plot.tag = element_text(face = "bold", size=14))
peuk

#compile plots and save it
all_legend=get_legend(S16)

grid1=ggarrange(p,S16,ncol = 1, nrow = 2, heights = c(3,10))
grid11=annotate_figure(grid1, fig.lab="Bacteria", fig.lab.face = "bold", fig.lab.pos = "top.left")
grid2=ggarrange(peuk,S18,ncol = 1, nrow = 2, heights = c(3,10))
grid22=annotate_figure(grid2, fig.lab="Microeukaryotes", fig.lab.face = "bold", fig.lab.pos = "top.left")


ggarrange(all_legend, 
          ggarrange(p, peuk,S16, S18, ncol=2, nrow=2, heights = c(3,10)), nrow=2, heights=c(0.5,10))


quartz.save("Assembly_facets_MS21.pdf", type="pdf", width=8, height=15)

dev.size()
dev.off()


##### Complementary tests after QPE
### Distance-decay relationships (Abundance-based Raup-Crick beta-diversity vs. spatial distance)
library(vegan)
coord.dist=read.csv("./coordinates.csv") #GPS coordinates for each pool
spatial_dist=vegdist(coord.dist, "euc") 

RC=read.csv("results919.csv", header=T, row.names=1) # for the Mantel test, I used the abundance-based Raup-Crick beta-diversity; load the data separately for each time point
RC[!lower.tri(RC)] <- NA
RC[is.na(RC)] <- " "
RC=RC[,1:14]
plot(otu_dist~spatial_dist)
plot(as.dist(RC)~spatial_dist)
mantel(as.dist(RC),spatial_dist)
mantel(otu_dist,spatial_dist)

### Distance-decay relationships (Abundance-based Raup-Crick beta-diversity vs. environmental distance)
setwd("/Volumes/FREECOM HDD/Phd_project/Assemblymechanisms/18S/results")
RC=read.csv("results826.csv", header=T, row.names=1) # for the Mantel test, I used the abundance-based Raup-Crick beta-diversity; load the data separately for each time point
RC[!lower.tri(RC)] <- NA
RC[is.na(RC)] <- " "

mm=cbind(RC,tp826) # get the corresponding environmental data for the above specified time point (e.g. 08.26)

env.16S=mm[,c("Cond","Water_temp","TP","TN","Daphnia", "Depth","Copepod")] #subsetting for the structuring environmental variables that were selected previously by RDA
env=data.frame(scale(env.16S, center = T, scale = T))

env_dist=vegdist(env, "euc")

plot(as.dist(RC)~env_dist)
mantel(as.dist(RC),env_dist)

coord.dist=mm[,c("lat","long")]
spatial_dist=vegdist(coord.dist, "euc")
mantel(as.dist(RC),spatial_dist)

#partial Mantel tests
mantel.partial(as.dist(RC), env_dist, spatial_dist, method = "pearson", permutations = 999)
mantel.partial(as.dist(RC), spatial_dist,env_dist, method = "pearson", permutations = 999)

# I collected the outcomes of Mantel tests in an Excel doc and complied in a way that it is easy to plot
# Load the results of the Mantel tests
DDenv=read.csv("/Volumes/FREECOM HDD/Phd_project/ENVdecay_16S_RCbray.csv", header=T, sep=";")
DDspa=read.csv("/Volumes/FREECOM HDD/Phd_project/Spatialdecay_RCbray.csv", header=T, sep=";")

library(ggplot2)
env <- ggplot(DDenv, aes(x = as.Date(Date), y = rM, group = group) ) +
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_point(aes(colour = cut(p, c(-Inf, 0.05, 0.051, Inf))),size=3) +
  scale_color_manual(name = "",
                     values = c("(-Inf,0.05]" = "black",
                                "(0.051, Inf]" = "#747474"),
                     labels = c("significant", "non-significant"))+
  facet_wrap(group~., scales = "fixed", nrow = 2, 
             strip.position = "right")+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")+
  theme_linedraw(base_size = 10)+theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90),
                                       strip.text = element_text(size=10), axis.text = element_text(size=10))+
  labs(x="Sampling occasion")+ labs(subtitle = "      Dry period             Wet period", tag="A")+
  theme(plot.subtitle = element_text(face="bold"),
        plot.tag = element_text(face="bold"),
        legend.text = element_text(size=10), legend.position = "bottom")+
  ylab(expression(atop("Mantel correlation", paste(beta[RCbray]," vs. environmental distance"))))
env

spa <- ggplot(DDspa, aes(x = as.Date(Date), y = rM, group = group) ) +
  geom_vline(xintercept = as.numeric(as.Date("2015-09-01")), linetype="dashed", color = "black", size=0.5, alpha=0.8)+
  geom_point(aes(colour = cut(p, c(-Inf, 0.05, 0.051, Inf))),size=3) +
  scale_color_manual(name = "",
                     values = c("(-Inf,0.05]" = "black",
                                "(0.051, Inf]" = "#747474"),
                     labels = c("significant", "non-significant"))+
  facet_wrap(group~., scales = "fixed", nrow = 2, 
             strip.position = "right")+
  scale_x_date(breaks = break.vec, date_labels =  "%b %d")+
  theme_linedraw(base_size = 10)+theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90),
                                       strip.text = element_text(size=10), axis.text = element_text(size=10),
                                       legend.text = element_text(size=10),
                                       legend.position = "bottom")+
  labs(x="Sampling occasion")+ labs(subtitle = "      Dry period              Wet period", tag="B")+
  theme(plot.subtitle = element_text(face="bold"),
        plot.tag = element_text(face="bold"))+
  ylab(expression(atop("Mantel correlation", paste(beta[RCbray]," vs. spatial distance"))))
spa


library(ggpubr)
devtools::install_url("https://github.com/wilkelab/cowplot/archive/0.6.3.zip")
library(cowplot)
ggarrange(env,spa, common.legend = T, legend = "bottom")

quartz.save("Distancedecay_tagged.pdf", type="pdf")
quartz.save("Distancedecay_tagged.tiff", type="tiff")


# End of the script
