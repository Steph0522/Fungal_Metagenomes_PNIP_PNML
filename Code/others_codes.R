
second<- plot_grid(first, leg, nrow = 1, ncol = 2, rel_widths = c(1,0.1),
                   align = "hv", axis = "r" )

plot_distance<- lapply(bc.dist.tidy.filt2, distance.plot0)
plot_distance_env<- lapply(bc.dist.tidy.filt2, distance.plot1)

third<- plot_grid(
  plot_distance[[1]]+ ggtitle(stats_qiime2$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10),   
  plot_distance[[2]]+ ggtitle(stats_single$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10),  
  plot_distance[[3]]+ ggtitle(stats_paired$label)+theme(plot.title =  element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10), 
  plot_distance[[4]]+  ggtitle(stats_kraken$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -70)),#+theme(aspect.ratio =20/10), 
  nrow = 1)

blank_plot<-ggplot() + theme_void()+ggtitle("Environmental distance")
leg1<- get_title(blank_plot)
leg2<- plot_grid(NULL, NULL,leg1, NULL, ncol = 4)

fourth<- plot_grid(
  plot_distance_env[[1]]+  ggtitle(stats_qiime2_e$label)+theme(plot.title= element_text(size = 10,hjust = 0.1,vjust = -100), axis.title.x = element_blank()),
  plot_distance_env[[2]]+ ggtitle(stats_single_e$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),
  plot_distance_env[[3]]+ ggtitle(stats_paired_e$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title.x = element_blank()),  
  plot_distance_env[[4]]+  ggtitle(stats_kraken_e$label)+theme(plot.title = element_text(size = 10,hjust = 0.1,vjust = -100),axis.title = element_blank()),  
  nrow = 2,ncol = 2, labels = c("A) QIIME2", "B) SINGLE MICOP",
                                "C) PAIRED MICOP", "D) KRAKEN2"))

four<- plot_grid(fourth, leg2, nrow = 2, rel_heights = c(1,0.05))
four
#plot_grid(get_legend(b), NA, b,c, d,e, nrow = 3, rel_widths = c(5,5), labels = c("a", NA, "b","c", "d", "e"), label_size = 15)

all<-plot_grid(  fourth,third ,nrow = 2, axis = "lr",align = "hv",
                 rel_heights = c(1, 1))

all
#ggsave("Fig1.bray_distance_env.tiff",width = 8, height = 8.2, dpi = 300, plot = four, device = "tiff")



mantels<- data.frame(
  Method=c("QIIME2", "SINGLE MICOP", "PAIRED MICOP", "KRAKEN2"),
  r = c(mantel_qiime2[[3]], mantel_single[[3]], 
        mantel_paired[[3]], mantel_kraken[[3]]),
  p = c(mantel_qiime2[[4]], mantel_single[[4]], 
        mantel_paired[[4]], mantel_kraken[[4]])) %>% mutate_at(c("r"), funs(round(.,digits = 2)))


library(ggpubr)
mantel_q1<-mantels %>% ggtexttable(rows = NULL, theme = ttheme("blank"))%>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 2)

mantel_q2<-mantels %>% ggtexttable(rows = NULL, theme = ttheme("blank"))%>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)%>%
  tab_add_hline(at.row = 5, row.side = "bottom", linewidth = 2)

library(cowplot)
man<-plot_grid(mantel_q1,  mantel_q2, 
               labels = c("a) Bray curtis", "b)  Horn"), nrow = 1)


ver<-plot_grid(fourth,fourth, nrow = 2)
man
#ggsave("ver.png",width = 6, height = 10, dpi = 300, plot = ver, device = "png")
