#comparando

meta<-read_excel("Data/Metadatos.xlsx") %>% mutate_if( 
  is.numeric, as.factor)%>% mutate(
    Pol= case_when(
      Poligono==1~ 1,
      Poligono==2~ 2,
      Poligono==3~ 3,
      Poligono==4~ 4,
      Poligono==5~ 5,
      Poligono==6~ 6))

table_single_micop_beta<- read.delim("Data/table_micop_single.txt") %>% t() %>% as.data.frame() %>% mutate_if(is.numeric,  ~ifelse(.x==0, 0, 1)) 
table_paired_micop_beta<- read.delim("Data/table_micop_paired.txt") %>% t() %>% as.data.frame()%>% mutate_if(is.numeric,  ~ifelse(.x==0, 0, 1)) 
table_qiime2_beta<- data.frame(read_qza("Data/clustered_table_filter.qza")$data, check.names = F) %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "SampleID") %>% separate(
    SampleID, c(
      "id_metagenome", "R", "unmap", "Paired"), 
    sep = "_")%>% inner_join(meta) %>% dplyr::select(
      -id_metagenome:-Paired, -id_sequence:-id_fisicoq) %>% 
  column_to_rownames(
    var="SampleID") %>% mutate_if(is.numeric,  ~ifelse(.x==0, 0, 1)) %>% dplyr::select(-Sites,-Pol) %>% mutate_all(as.numeric)

library(hilldiv)
fung.dist.hilldiv0<- hilldiv::pair_dis(t(fungb), qvalue = 0)
fung.dist.hilldiv01<- hilldiv::pair_dis(t(fungb1), qvalue = 0)
fung.dist.hilldiv02<- hilldiv::pair_dis(t(fungb2), qvalue = 0)

#jaccard 

fung.dist.jac<-as.dist(fung.dist.hilldiv0$L1_UqN)
fung.dist.jac1<-as.dist(fung.dist.hilldiv01$L1_UqN)
fung.dist.jac2<-as.dist(fung.dist.hilldiv02$L1_UqN)

fung.dist.jtu<-as.dist(1-fung.dist.hilldiv0$L1_SqN)
fung.dist.jtu1<-as.dist(1-fung.dist.hilldiv01$L1_SqN)
fung.dist.jtu2<-as.dist(1-fung.dist.hilldiv02$L1_SqN)

fung.dist.jne<-as.dist(fung.dist.jac-fung.dist.jtu)
fung.dist.jne1<-as.dist(abs(fung.dist.jac1-fung.dist.jtu1))
fung.dist.jne2<-as.dist(abs(fung.dist.jac2-fung.dist.jtu2))

mat<- as.matrix(fung.dist.jac)
envi<- meta[match(colnames(mat), meta$SampleID),]


jd<-betadisper(fung.dist.jac,factor(envi$Pol))
jd1<-betadisper(fung.dist.jac1,factor(envi$Pol))
jd2<-betadisper(fung.dist.jac2,factor(envi$Pol))

td<-betadisper(fung.dist.jtu,factor(envi$Pol))
td1<-betadisper(fung.dist.jtu1,factor(envi$Pol))
td2<-betadisper(fung.dist.jtu2,factor(envi$Pol))

nd<-betadisper(fung.dist.jne,factor(envi$Pol))
nd1<-betadisper(fung.dist.jne1,factor(envi$Pol))
nd2<-betadisper(fung.dist.jne2,factor(envi$Pol))

adj<-adonis2(fung.dist.jac~envi$Pol)
adj1<-adonis2(fung.dist.jac1~envi$Pol)
adj2<-adonis2(fung.dist.jac2~envi$Pol)

adt<-adonis2(fung.dist.jtu~envi$Pol)
adt1<-adonis2(fung.dist.jtu1~envi$Pol)
adt2<-adonis2(fung.dist.jtu2~envi$Pol)

adn<-adonis2(fung.dist.jne~envi$Pol)
adn1<-adonis2(fung.dist.jne1~envi$Pol)
adn2<-adonis2(fung.dist.jne2~envi$Pol)



library(tidyverse)
library(ggordiplots)
function_plot_beta<- function(x){
  y <- gg_ordiplot(x, groups = envi$Poligono, hull = FALSE, spiders = TRUE, 
                   ellipse = FALSE, plot = FALSE, label = TRUE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group))+guides(
      color=guide_legend(title="Polygon"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Polygon")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  return(a)
}

function_plot_beta_micop<- function(x){
  y <- gg_ordiplot(x, groups = envi$Poligono, hull = FALSE, spiders = TRUE, 
                   ellipse = FALSE, plot = FALSE, label = TRUE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group))+guides(
      color=guide_legend(title="Polygon"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Polygon")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  return(a)
}

#jaccard
jaccard_qiime2 <- function_plot_beta(jd)
jac_qiime2 <- jaccard_qiime2+labs(title = paste("F = ",signif(adj$F[1], 3), ",",
                                                "p-value= ",signif(adj$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)
jaccard_single <- function_plot_beta_micop(jd1)
jac_single <- jaccard_single+labs(title = paste("F = ",signif(adj1$F[1], 3), ",",
                                                "p-value= ",signif(adj1$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)


jaccard_paired <- function_plot_beta_micop(jd2)
jac_paired <- jaccard_paired+labs(title = paste("F = ",signif(adj1$F[1], 3), ",",
                                                "p-value= ",signif(adj1$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)


#turnover
jacturn_qiime2 <- function_plot_beta(td)
turnqiime2 <- jacturn_qiime2+labs(title = paste("F = ",signif(adt$F[1], 3), ",",
                                                "p-value= ",signif(adt$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)
jacturn_single <- function_plot_beta_micop(td1)
turnsingle <- jacturn_single+labs(title = paste("F = ",signif(adt1$F[1], 3), ",",
                                                "p-value= ",signif(adt1$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)


jacturn_paired <- function_plot_beta_micop(td2)
turnpaired <- jacturn_paired+labs(title = paste("F = ",signif(adt1$F[1], 3), ",",
                                                "p-value= ",signif(adt1$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)

#nestedness
jacnest_qiime2 <- function_plot_beta(nd)
nestqiime2 <- jacnest_qiime2+labs(title = paste("F = ",signif(adn$F[1], 3), ",",
                                                "p-value= ",signif(adn$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)
jacnest_single <- function_plot_beta_micop(nd1)
nestsingle <- jacnest_single+labs(title = paste("F = ",signif(adn1$F[1], 3), ",",
                                                "p-value= ",signif(adn1$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)


jacnest_paired <- function_plot_beta_micop(nd2)
nestpaired <- jacnest_paired+labs(title = paste("F = ",signif(adn1$F[1], 3), ",",
                                                "p-value= ",signif(adn1$`Pr(>F)`[1], 5)))+ labs(x = "PCo 1", y = "PCo 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =6/10)


library(cowplot)
leg<- get_legend(jac_qiime2)
#a<-plot_grid(jac_qiime2+theme(legend.position = "none")+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+xlab(""), NULL,jac_single+theme(legend.position = "none")+ylab("")+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+xlab(""), NULL,jac_paired+theme(legend.position = "none")+ylab("")+xlab("")+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 3,nrow = 1 ,align = "hv", rel_widths = c(1,0,1,0,1))
#b<-plot_grid(turnqiime2+theme(legend.position = "none")+xlab(""), turnsingle+theme(legend.position = "none")+ylab("")+xlab(""), turnpaired+theme(legend.position = "none")+ylab("")+xlab(""), ncol = 3, align = "hv")
#c<-plot_grid(nestqiime2+theme(legend.position = "none"), nestsingle+theme(legend.position = "none")+ylab(""), nestpaired+theme(legend.position = "none")+ylab(""), ncol = 3, align = "hv")
d<- plot_grid(NULL, leg, NULL, ncol = 3)
#plot_grid(a,b,c, nrow = 3)

a<-plot_grid(jac_qiime2+theme(legend.position = "none")+xlab(""), jac_single+theme(legend.position = "none")+ylab("")+xlab(""), jac_paired+theme(legend.position = "none")+ylab("")+xlab(""), 
             turnqiime2+theme(legend.position = "none")+xlab(""), turnsingle+theme(legend.position = "none")+ylab("")+xlab(""), turnpaired+theme(legend.position = "none")+ylab("")+xlab(""),
             nestqiime2+theme(legend.position = "none"), nestsingle+theme(legend.position = "none")+ylab(""), nestpaired+theme(legend.position = "none")+ylab(""), ncol = 3, nrow = 3 ,align = "hv", greedy = FALSE)

e<-plot_grid(a, d, rel_widths = c(10,1))
e
ggsave("gao_beta_hill.pdf",width = 16, height = 10, dpi = 300, plot = e, device = "pdf")


