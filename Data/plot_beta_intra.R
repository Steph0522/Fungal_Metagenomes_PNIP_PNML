
intra_ro_wet_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=0-table_wet_ro.txt.csv") %>% mutate(qs="q0") %>% mutate(type="wet")
intra_ro_wet_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=1-table_wet_ro.txt.csv") %>% mutate(qs="q1")%>% mutate(type="wet")
intra_ro_wet_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=2-table_wet_ro.txt.csv") %>% mutate(qs="q2")%>% mutate(type="wet")

intra_ri_wet_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=0-table_wet_ri.txt.csv") %>% mutate(qs="q0")%>% mutate(type="wet")
intra_ri_wet_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=1-table_wet_ri.txt.csv") %>% mutate(qs="q1")%>% mutate(type="wet")
intra_ri_wet_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=2-table_wet_ri.txt.csv") %>% mutate(qs="q2")%>% mutate(type="wet")

intra_nr_wet_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=0-table_wet_nr.txt.csv") %>% mutate(qs="q0")%>% mutate(type="wet")
intra_nr_wet_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=1-table_wet_nr.txt.csv") %>% mutate(qs="q1")%>% mutate(type="wet")
intra_nr_wet_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=2-table_wet_nr.txt.csv") %>% mutate(qs="q2")%>% mutate(type="wet")


intra_ro_dry_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=0-table_dry_ro.txt.csv") %>% mutate(qs="q0")%>% mutate(type="dry")
intra_ro_dry_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=1-table_dry_ro.txt.csv") %>% mutate(qs="q1")%>% mutate(type="dry")
intra_ro_dry_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=2-table_dry_ro.txt.csv") %>% mutate(qs="q2")%>% mutate(type="dry")

intra_ri_dry_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=0-table_dry_ri.txt.csv") %>% mutate(qs="q0")%>% mutate(type="dry")
intra_ri_dry_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=1-table_dry_ri.txt.csv") %>% mutate(qs="q1")%>% mutate(type="dry")
intra_ri_dry_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=2-table_dry_ri.txt.csv") %>% mutate(qs="q2")%>% mutate(type="dry")

intra_nr_dry_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=0-table_dry_nr.txt.csv") %>% mutate(qs="q0")%>% mutate(type="dry")
intra_nr_dry_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=1-table_dry_nr.txt.csv") %>% mutate(qs="q1")%>% mutate(type="dry")
intra_nr_dry_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=2-table_dry_nr.txt.csv") %>% mutate(qs="q2")%>% mutate(type="dry")


intra_ro_exdry_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=0-table_exdry_ro.txt.csv") %>% mutate(qs="q0")%>% mutate(type="extreme-dry")
intra_ro_exdry_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=1-table_exdry_ro.txt.csv") %>% mutate(qs="q1")%>% mutate(type="extreme-dry")
intra_ro_exdry_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityro/Intra_Beta_Similarity-q=2-table_exdry_ro.txt.csv") %>% mutate(qs="q2")%>% mutate(type="extreme-dry")

intra_ri_exdry_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=0-table_exdry_ri.txt.csv") %>% mutate(qs="q0")%>% mutate(type="extreme-dry")
intra_ri_exdry_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=1-table_exdry_ri.txt.csv") %>% mutate(qs="q1")%>% mutate(type="extreme-dry")
intra_ri_exdry_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversityri/Intra_Beta_Similarity-q=2-table_exdry_ri.txt.csv") %>% mutate(qs="q2")%>% mutate(type="extreme-dry")

intra_nr_exdry_q0<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=0-table_exdry_nr.txt.csv") %>% mutate(qs="q0")%>% mutate(type="extreme-dry")
intra_nr_exdry_q1<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=1-table_exdry_nr.txt.csv") %>% mutate(qs="q1")%>% mutate(type="extreme-dry")
intra_nr_exdry_q2<- read.csv("/home/yendi/Downloads/men12923-sup-0002-supinfo/OSM-I-DiversityCode/Beta_diversitynr/Intra_Beta_Similarity-q=2-table_exdry_nr.txt.csv") %>% mutate(qs="q2")%>% mutate(type="extreme-dry")

library(ggpubr)

intra_ro<- rbind(intra_ro_wet_q0, intra_ro_dry_q0, intra_ro_exdry_q0,
                 intra_ro_wet_q1, intra_ro_dry_q1, intra_ro_exdry_q1,
                 intra_ro_wet_q2, intra_ro_dry_q2, intra_ro_exdry_q2)
i1<-intra_ro %>% ggpubr::ggboxplot(x = "type", y="TurnoverComp", fill = "type")+facet_wrap(
  ~qs, scales = "free")+stat_compare_means()+xlab("")+ylab("TurnOver-Roots")+scale_fill_manual(values = c("DarkGreen", "yellow", "red"))
intra_ri<- rbind(intra_ri_wet_q0, intra_ri_dry_q0, intra_ri_exdry_q0,
                 intra_ri_wet_q1, intra_ri_dry_q1, intra_ri_exdry_q1,
                 intra_ri_wet_q2, intra_ri_dry_q2, intra_ri_exdry_q2)
i2<-intra_ri %>% ggpubr::ggboxplot(x = "type", y="TurnoverComp", fill = "type")+facet_wrap(
  ~qs, scales = "free")+stat_compare_means()+xlab("")+ylab("TurnOver-Rhizosphere")+scale_fill_manual(values = c("DarkGreen", "yellow", "red"))
intra_nr<- rbind(intra_nr_wet_q0, intra_nr_dry_q0, intra_nr_exdry_q0,
                 intra_nr_wet_q1, intra_nr_dry_q1, intra_nr_exdry_q1,
                 intra_nr_wet_q2, intra_nr_dry_q2, intra_nr_exdry_q2)
i3<-intra_nr %>% ggpubr::ggboxplot(x = "type", y="TurnoverComp", fill = "type")+facet_wrap(
  ~qs, scales = "free")+stat_compare_means()+xlab("")+ylab("TurnOver-Bulk soil")+scale_fill_manual(values = c("DarkGreen", "yellow", "red"))

#intra_nr %>% ggpubr::ggbarplot(x = "type", y="TurnoverComp", fill = "type", stat="identity", add="mean_sd")+facet_wrap(~qs, scales = "free")
library(cowplot)
p<-plot_grid(i1+theme(legend.position = "none"),i2+theme(legend.position = "none"),i3+theme(legend.position = "none"), nrow = 3)
ggsave(plot=p, "/home/yendi/Desktop/intra-turnover.pdf", width = 16, height = 12)
