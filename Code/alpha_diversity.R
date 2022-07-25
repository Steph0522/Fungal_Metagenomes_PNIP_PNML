#Alpha diversity
library(tidyverse)
library(readxl)
library(qiime2R)
library(ggpubr)

metadata<-read_excel("Data/Metadatos.xlsx") %>% mutate_if(is.numeric, as.factor)
table_single_micop<- read.delim("Data/table_micop_single.txt") 
table_paired_micop<- read.delim("Data/table_micop_paired.txt")
table_qiime2<- data.frame(read_qza("Data/clustered_table_filter.qza")$data, 
                          check.names = F) %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "SampleID") %>% separate(
    SampleID, c(
      "id_metagenome", "R", "unmap", "Paired"), 
    sep = "_")%>% inner_join(metadata) %>% dplyr::select(
      -id_metagenome:-Paired, -id_sequence:-id_fisicoq) %>% column_to_rownames(
        var="SampleID") %>% t() %>% as.data.frame()

library(hillR)
q0<-hill_taxa(table_qiime2, q = 0, MARGIN = 2)
q1<-hill_taxa(table_qiime2, q = 1, MARGIN = 2)
q2<-hill_taxa(table_qiime2, q = 2, MARGIN = 2)

alphadiv_qiime2<- cbind(q0, q1, q2) %>%as.data.frame(
) %>% rownames_to_column(var = "SampleID") %>%  pivot_longer(
  .,cols = q0:q2, names_to =  "qs", values_to =  "value") %>% inner_join(
    metadata) %>% mutate(order=case_when(
      qs=="q0"~"q = 0",
      qs=="q1"~"q = 1",
      qs=="q2"~"q = 2"  ))

q0<-hill_taxa(table_single_micop, q = 0, MARGIN = 2)
q1<-hill_taxa(table_single_micop, q = 1, MARGIN = 2)
q2<-hill_taxa(table_single_micop, q = 2, MARGIN = 2)

alphadiv_single_micop<- cbind(q0, q1, q2) %>%as.data.frame(
) %>% rownames_to_column(var = "SampleID") %>%  pivot_longer(
  .,cols = q0:q2, names_to =  "qs", values_to =  "value") %>% inner_join(
    metadata) %>% mutate(order=case_when(
      qs=="q0"~"q = 0",
      qs=="q1"~"q = 1",
      qs=="q2"~"q = 2"  ))

q0<-hill_taxa(table_paired_micop, q = 0, MARGIN = 2)
q1<-hill_taxa(table_paired_micop, q = 1, MARGIN = 2)
q2<-hill_taxa(table_paired_micop, q = 2, MARGIN = 2)

alphadiv_paired_micop<- cbind(q0, q1, q2) %>%as.data.frame(
) %>% rownames_to_column(var = "SampleID") %>%  pivot_longer(
  .,cols = q0:q2, names_to =  "qs", values_to =  "value") %>% inner_join(
    metadata) %>% mutate(order=case_when(
      qs=="q0"~"q = 0",
      qs=="q1"~"q = 1",
      qs=="q2"~"q = 2"  ))

#qiime2
alphadiv_qiime2_q0<- alphadiv_qiime2 %>%filter(qs=="q0")
alphadiv_qiime2_q1<- alphadiv_qiime2 %>%filter(qs=="q1") 
alphadiv_qiime2_q2<- alphadiv_qiime2 %>%filter(qs=="q2") 

stat.test_qiime2_q0 <- compare_means(
  value ~ Poligono, data = alphadiv_qiime2_q0) %>%
  mutate(y.position = seq(460,600, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)
stat.test_qiime2_q1 <- compare_means(
  value ~ Poligono, data = alphadiv_qiime2_q1) %>%
  mutate(y.position = seq(260,400, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)
stat.test_qiime2_q2 <- compare_means(
  value ~ Poligono, data = alphadiv_qiime2_q2) %>%
  mutate(y.position = seq(460,600, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)

aq0<- ggboxplot(alphadiv_qiime2_q0, 
                x = "Poligono", y = "value", fill = "Poligono", facet.by = "order") +
  stat_pvalue_manual(stat.test_qiime2_q0, label = "p.signif", hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         strip.text.x = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")
aq1<- ggboxplot(alphadiv_qiime2_q1, 
                x = "Poligono", y = "value", fill = "Poligono", facet.by = "order") +
  stat_pvalue_manual(stat.test_qiime2_q1, label = "p.signif", hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         strip.text.x = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")

aq2<- ggboxplot(alphadiv_qiime2_q2 %>% mutate(methods="QIIME 2"), 
                x = "Poligono", y = "value", fill = "Poligono", facet.by = "order") +
  stat_pvalue_manual(stat.test_qiime2_q2, label = "p.signif", hide.ns = TRUE)+
  facet_grid(rows = vars(methods), cols = vars(order))+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")

qiime2al<-cowplot::plot_grid(aq0, aq1, aq2, ncol = 3, align = "h", rel_widths = c(1,1,1.1))

#single_micop
alphadiv_single_micop_q0<- alphadiv_single_micop %>%filter(qs=="q0")
alphadiv_single_micop_q1<- alphadiv_single_micop %>%filter(qs=="q1") 
alphadiv_single_micop_q2<- alphadiv_single_micop %>%filter(qs=="q2") 

stat.test_single_micop_q0 <- compare_means(
  value ~ Poligono, data = alphadiv_single_micop_q0) %>%
  mutate(y.position = seq(160,300, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)
stat.test_single_micop_q1 <- compare_means(
  value ~ Poligono, data = alphadiv_single_micop_q1) %>%
  mutate(y.position = seq(60,200, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)
stat.test_single_micop_q2 <- compare_means(
  value ~ Poligono, data = alphadiv_single_micop_q2) %>%
  mutate(y.position = seq(60,200, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)

as0<- ggboxplot(alphadiv_single_micop_q0, 
                x = "Poligono", y = "value", fill = "Poligono") +
  stat_pvalue_manual(stat.test_single_micop_q0, label = "p.signif", hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         strip.text.x = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")
as1<- ggboxplot(alphadiv_single_micop_q1, 
                x = "Poligono", y = "value", fill = "Poligono") +
  stat_pvalue_manual(stat.test_single_micop_q1, label = "p.signif", hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         strip.text.x = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")

as2<- ggboxplot(alphadiv_single_micop_q2 %>% mutate(methodS="SINGLE MICOP"), 
                x = "Poligono", y = "value", fill = "Poligono") +
  facet_grid( rows= vars(methodS))+
  stat_pvalue_manual(stat.test_single_micop_q2, label = "p.signif",
                     hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")

singleal<-cowplot::plot_grid(as0, as1, as2, ncol = 3, align = "h", rel_widths = c(1,1,1.1))

#paired_micop
alphadiv_paired_micop_q0<- alphadiv_paired_micop %>%filter(qs=="q0")
alphadiv_paired_micop_q1<- alphadiv_paired_micop %>%filter(qs=="q1") 
alphadiv_paired_micop_q2<- alphadiv_paired_micop %>%filter(qs=="q2") 

stat.test_paired_micop_q0 <- compare_means(
  value ~ Poligono, data = alphadiv_paired_micop_q0) %>%
  mutate(y.position = seq(160,300, by=10)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)
stat.test_paired_micop_q1 <- compare_means(
  value ~ Poligono, data = alphadiv_paired_micop_q1) %>%
  mutate(y.position = seq(0,28, by=2)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)
stat.test_paired_micop_q2 <- compare_means(
  value ~ Poligono, data = alphadiv_paired_micop_q2) %>%
  mutate(y.position = seq(0,28, by=2)) %>% mutate(
    "p.format2"=p.adj,"p.adj2"=p.format) %>% dplyr::select(
      -p.adj, -p.format) %>% rename(p.adj=p.adj2, p.format=p.format2)

ap0<- ggboxplot(alphadiv_paired_micop_q0, 
                x = "Poligono", y = "value", fill = "Poligono") +
  stat_pvalue_manual(stat.test_paired_micop_q0, label = "p.signif", hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_text(size = 12),
                         strip.text.x = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")
ap1<- ggboxplot(alphadiv_paired_micop_q1, 
                x = "Poligono", y = "value", fill = "Poligono") +
  stat_pvalue_manual(stat.test_paired_micop_q1, label = "p.signif", hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_text(size = 12),
                         strip.text.x = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")

ap2<- ggboxplot(alphadiv_paired_micop_q2 %>% mutate(methodS="PAIRED MICOP"), 
                x = "Poligono", y = "value", fill = "Poligono") +
  facet_grid(rows = vars(methodS))+
  stat_pvalue_manual(stat.test_paired_micop_q2, label = "p.signif",
                     hide.ns = TRUE)+
  theme_linedraw()+theme(legend.position = "none", 
                         axis.title = element_blank(),
                         axis.text.x = element_text(size = 12),
                         strip.text = element_text(size = 12, face = "bold.italic"),
                         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                         colour = "#E5E8E8"), 
                         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                         colour = "#E5E8E8"))+
  scale_fill_viridis_d(option ="turbo", name="Polygon")

pairedal<-cowplot::plot_grid(ap0, ap1, ap2, ncol = 3, align = "h", rel_widths = c(1,1,1.1))


cowplot::plot_grid(qiime2al, singleal, pairedal, nrow = 3)

all<- cowplot::plot_grid(aq0, aq1, aq2,
                         as0, as1, as2,
                         ap0, ap1, ap2, align = "hv", labels = c("a)", "b)", "c)",
                                                                 "d)", "e)", "f)",
                                                                 "g)", "h)", "i)"), 
                          scale = 0.9)
all
ggsave("Figures_final/Fig1.exploratory_alpha.pdf",width = 11, height = 7, dpi = 300, plot = all, device = "pdf")


