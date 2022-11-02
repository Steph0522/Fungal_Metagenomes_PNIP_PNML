#OVERLAP
#Load libraries
library(broom)
library(cowplot)
library(tidyverse)
#source("../General/general_functions.R")


#Load data
#Load occupancy data frame
#Load vOTU overlap data
overlap.tidy <- pct.tidy

#Panel B
#Filter pairwise comparisons involving different time points
overlap.filt <- overlap.tidy %>% 
  filter(SampleID.x != SampleID.y) %>% 
  full_join(distance_dm) %>% rename(SpatialDistance = value) %>% 
  mutate(PercentShared = PercentShared * 100)

#Perform correlation analysis and regression
cor <- cor.test(overlap.filt$SpatialDistance, overlap.filt$PercentShared, method = "pearson", alternative = "two.sided") %>% tidy()
lm <- lm(PercentShared ~ SpatialDistance, data = overlap.filt) %>% tidy() %>% filter(term == "SpatialDistance")
overlap.stats <- data.frame(label = paste("r = ", signif(cor$estimate,3), 
                                          "\nslope = ", signif(lm$estimate, 3),
                                          "\nP = ", signif(cor$p.value, 3)))

#Plot
overlap.p <- overlap.filt %>% 
  ggplot(aes(SpatialDistance, PercentShared)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "gray25") +
  geom_text(data = overlap.stats, aes(x = 0, y = 30, label = label), hjust = 0) + 
  geom_smooth(color = "black", se = F, method = "lm") +
  # scale_x_continuous(breaks = seq(0, 18, by = 3)) +
  xlab("Distance between plots (m)") +
  ylab("Community overlap\n(% shared vOTUs)") +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "top")

overlap.p



#plot_grid(occupancy.p, overlap.p, nrow = 1, align = "v", axis = "lr", labels = c("a", "b"), label_size = 15)


