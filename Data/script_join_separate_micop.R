library(tidyverse)
library(ggpubr)
abund_micop_0<- read.delim("/home/yendi/Downloads/MiCoP/abundances_cut0.txt", check.names = F) %>% rename(cut0 = PERCENTAGE)
#abund_micop_15<- read.delim("/home/yendi/Downloads/MiCoP/abundances_cut15.txt", check.names = F) %>% rename(cut15 = PERCENTAGE)
#abund_micop_10<- read.delim("/home/yendi/Downloads/MiCoP/abundances_cut10.txt", check.names = F)%>% rename(cut10 = PERCENTAGE)
#abund_micop_20<- read.delim("/home/yendi/Downloads/MiCoP/abundances_cut20.txt", check.names = F)%>% rename(cut20 = PERCENTAGE)
#abund_micop_50<- read.delim("/home/yendi/Downloads/MiCoP/abundances_cut50.txt", check.names = F)%>% rename(cut50 = PERCENTAGE)
#abund_micop_100<- read.delim("/home/yendi/Downloads/MiCoP/abundances_cut100.txt", check.names = F)%>% rename(cut100= PERCENTAGE)
#abund_micop_r1_100<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r1_single_100.txt", check.names = F) %>% rename(r1cut100 = PERCENTAGE)
#abund_micop_r1_50<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r1_single_50.txt", check.names = F) %>% rename(r1rcut50 = PERCENTAGE)
#abund_micop_r1_15<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r1_single_15.txt", check.names = F) %>% rename(r1cut15 = PERCENTAGE)
#abund_micop_r2_50<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r2_single_50.txt", check.names = F) %>% rename(r2rcut50 = PERCENTAGE)
#abund_micop_r2_15<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r2_single_15.txt", check.names = F) %>% rename(r2cut15 = PERCENTAGE)
abund_micop_r2_0<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r2_single_0.txt", check.names = F) %>% rename(r2cut0 = PERCENTAGE)
abund_micop_r1_0<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r1_single_0.txt", check.names = F) %>% rename(r1cut0 = PERCENTAGE)
#abund_micop_r1_10<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r1_single_10.txt", check.names = F) %>% rename(r1cut10 = PERCENTAGE)
#abund_micop_r2_10<- read.delim("/home/yendi/Downloads/MiCoP/abundace_r2_single_10.txt", check.names = F) %>% rename(r2cut10 = PERCENTAGE)
abund_micop_both_0<- read.delim("/home/yendi/Downloads/MiCoP/abundances_both_0.txt", check.names = F) %>% rename(bothcut0 = PERCENTAGE)
abund_micop_paired_0<- read.delim("/home/yendi/Downloads/MiCoP/abundances_paired_script_0.txt", check.names = F) %>% rename(pairedcut0 = PERCENTAGE)


all_list<- list(abund_micop_0, abund_micop_10,abund_micop_15, abund_micop_20, abund_micop_50, abund_micop_100,
                abund_micop_r1_0, abund_micop_r1_10,abund_micop_r1_15, abund_micop_r1_50, abund_micop_r1_100,
                abund_micop_r2_0, abund_micop_r2_10,abund_micop_r2_15, abund_micop_r2_50)
all_list<- list(abund_micop_0, abund_micop_r1_0, abund_micop_r2_0, abund_micop_both_0, abund_micop_paired_0)                

all<- abund_micop_0 %>% full_join(abund_micop_10) %>% full_join(
  abund_micop_15)%>% full_join(abund_micop_20) %>% full_join(abund_micop_50)  %>% replace(is.na(.), 0)

all<- all_list %>% reduce(full_join) %>% replace(is.na(.), 0)

phylum<- all %>% filter(RANK =="phylum") %>% dplyr::select(-1:-3) %>% rename(taxa=TAXPATHSN) %>% separate(taxa, into = c("taxo", "phylum"))

phylum %>% pivot_longer(., cut0:pairedcut0, names_to = "names", values_to = "values") %>% ggbarplot(
  ., x="names", y="values", fill = "phylum")+ theme(legend.position="right")

classes<- all %>% filter(RANK =="class") %>% dplyr::select(-1:-3) %>% rename(taxa=TAXPATHSN) %>% separate(taxa, into = c("taxo", "phylum", "class"))
classes %>% pivot_longer(., cut0:pairedcut0, names_to = "names", values_to = "values") %>% ggbarplot(
  ., x="names", y="values", fill = "class")+ theme(legend.position="right")


orders<- all %>% filter(RANK =="order") %>% dplyr::select(-1:-3) %>% rename(taxa=TAXPATHSN) %>% separate(taxa, into = c("taxo", "phylum", "class", "order"))
  orders %>% pivot_longer(., cut0:pairedcut0, names_to = "names", values_to = "values")  %>% top_n(50) %>% ggbarplot(
  ., x="names", y="values", fill = "order")+ theme(legend.position="right")
  
familys<- all %>% filter(RANK =="family") %>% dplyr::select(-1:-3) %>% rename(taxa=TAXPATHSN) %>% separate(taxa, into = c("taxo", "phylum", "class", "order", "family"))
  familys %>% pivot_longer(., cut0:pairedcut0, names_to = "names", values_to = "values") %>% top_n(70) %>%  ggbarplot(
    ., x="names", y="values", fill = "family")+ theme(legend.position="right")

genus<- all %>% filter(RANK =="genus") %>% dplyr::select(-1:-3) %>% rename(taxa=TAXPATHSN) %>% separate(taxa, into = c("taxo", "phylum", "class", "order", "family", "genus"))
  genus %>% pivot_longer(., cut0:pairedcut0, names_to = "names", values_to = "values") %>% top_n(80) %>%  ggbarplot(
    ., x="names", y="values", fill = "genus")+ theme(legend.position="right")
  
  