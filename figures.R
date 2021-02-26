library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggpmisc) # package to count measurements in quadrants
source('libs/func.R')
source('general.R')
# global Variables --------------------------------------------------------

TEs <- c('DNA', 'LINE', 'LTR', 'SINE')

Tools <- c('SalmonTE', 'Telescope', 'TEtranscripts', 'SQuIRE', 'TEtools')

# script relevant functions -----------------------------------------------

reOrderSetting <- function(df){
  
  df <- df %>% mutate(Setting = factor(Setting, levels = c('single', 'paired')))
  return(df)
}

reOrderSettingRev <- function(df){
  df <- df %>% mutate(Setting = factor(Setting, levels = c( 'paired','single')))
  return(df)
}

reOrderTools <- function(df){
  
  Tools <- c('salmonTE', 'Telescope', 'TEtranscripts', 'SQuIRE', 'TEtools')
  df <- df %>% mutate(Tool = factor(Tool, levels = Tools))
  return(df)
  
}


my_theme <- function(){
  
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white", color = "grey80"),
        strip.text = element_text(family = "CM Sans", size = 12),
        axis.text = element_text(family = "CM Sans", size = 10),
        axis.title = element_text(family = "CM Sans", size = 12),
        title = element_text(family = "CM Sans", size = 16),
        panel.spacing.x = unit(0.2, "lines"),
        panel.border =  element_rect(fill = "NA", color = "grey80"),
        legend.key         = element_rect(fill = "white"))
}


my_colors <- c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F", '#af7ac5', '#E09F3E', '#335C67')
names(my_colors) <- c('salmonTE', 'SQuIRE', 'TEtools', 'TEtranscripts', 'Telescope', 'single', 'paired')


# load general data -------------------------------------------------------
# The following data frames contain the raw counts of the tools. The count table
# are processed, which means that detected TEs with less than 5 reads in sum
# across all samples are removed. Furthermore the count tables are rearranged
# so that per detected TE per replicate, Tool and Setup a row exists where the
# recovered and the simulated count is stored and out of this a lot of 
# additional information were added. Out of this information I can calculated
# the binary measurement of TP, FP and FN.

df.single <- loadRdata(paste0(project, '/Data/single.combined.cnttbl.processed.Rdata'))
df.paired <- loadRdata(paste0(project, '/Data/paired.combined.cnttbl.processed.Rdata'))

df.all <- rbind(df.single, df.paired)


# F-Scores ----------------------------------------------------------------
#
# The F-Scores represent a trade-off between precision and recall. The higher
# the F-Score, the better the tool performs. The F-Score is calculated for
# different groups - all, small kimura, order, and order & small kimura.


## F-Scores - All ----------------------------------------------------------

metrics <- df.all %>% 
  count(Condition, Tool, sample, Simulation, Setting, name = 'count') %>%  
  spread(Condition, count) %>% 
  mutate(Sensitivity = TP/(TP+FN),
         Precision = TP/(TP+FP),
         F_score = 2*(Precision*Sensitivity)/(Precision+Sensitivity),
         FDR = FP/(TP+FP),
         Setting = as.factor(Setting),
         Tool = case_when(as.character(Tool) == 'salmonTE' ~ 'SalmonTE',
                          TRUE ~ Tool))

metrics.sum <- metrics %>% 
  group_by(Tool, Setting) %>%
  summarise(mean = mean(F_score), 
            sd = sd(F_score),
            se = sd/sqrt(n())) %>% 
  reOrderSetting()

metrics.sum$Setting <- ordered(metrics.sum$Setting, levels=c('paired', 'single'))
metrics.sum$Tool <- ordered(metrics.sum$Tool, levels = rev(Tools))

panel1.a <- ggplot(metrics.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge',
           aes(group=Setting), 
           color = 'black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2) +
  scale_y_continuous(limits = c(0,1.05), expand = c(0,0))+
  my_theme() +
  coord_flip() +
  scale_fill_manual(name = 'Setting', values = my_colors) +
  labs(y = 'mean(F-score)',
       title = 'all') +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size = 11, family = 'CM Sans', hjust = 1,  vjust = -1.5))

save(panel1.a, file = paste0(project, '/Figures/raw/figure_1a.raw.Rdata'))

png(paste0(project, '/Figures/figure_1a.png'), width = 8, height = 5, units = 'in', res = 300)
print(panel1.a)
dev.off()

write.csv(metrics.sum, file = paste0(project, '/Figures/tmp/figure_1a.csv'))



### Supplemental Figure F-Score by Kimura distance

metrics.kimura <- df.all %>% 
  count(Condition, Tool, sample, Simulation, Setting, Kim.grp, name = 'count') %>%  
  spread(Condition, count) %>% 
  mutate(Sensitivity = TP/(TP+FN),
         Precision = TP/(TP+FP),
         F_score = 2*(Precision*Sensitivity)/(Precision+Sensitivity),
         FDR = FP/(TP+FP),
         Setting = as.factor(Setting),
         Tool = case_when(as.character(Tool) == 'salmonTE' ~ 'SalmonTE',
                          TRUE ~ Tool))

metrics.sum.kimura <- metrics.kimura %>% 
  group_by(Tool, Setting, Kim.grp) %>%
  summarise(mean = mean(F_score), 
            sd = sd(F_score),
            se = sd/sqrt(n())) %>% 
  reOrderSetting()

metrics.kimura$Setting <- ordered(metrics.kimura$Setting, levels=c('single','paired'))
metrics.kimura$Tool <- ordered(metrics.kimura$Tool, levels = Tools)


pl.supplement.3 <- ggplot(metrics.kimura, aes(Kim.grp, F_score, color=Setting)) +
  geom_boxplot(width=0.5, position = "dodge") +
  facet_grid(Tool~Setting) +
  scale_color_manual(values = c('black', 'black')) +
  labs(y = 'F-score',
       x = 'Kimura distance intervals') +
  my_theme()+
  theme(panel.grid.major = element_line(colour = 'grey', linetype = 'dotted'),
        legend.position = 'none',
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(1, 'lines'),
        panel.spacing.y = unit(1, 'lines'),
        axis.text = element_text(family = "CM Sans", size = 16),
        axis.title = element_text(family = "CM Sans", size = 18),
        strip.text = element_text(family = "CM Sans", size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1)        )

png(paste0(project, "/Figures/figure_S2.png"), width = 12.5, height = 10, units = 'in', res = 600)
print(pl.supplement.3) # Make plot
dev.off()

ggsave(filename = paste0(project, "/Figures/figure_S2.pdf"), 
       print(pl.supplement.3),
       width = 12.5, height = 10, dpi = 600, units = "in", device=cairo_pdf)

save(pl.supplement.3, file = paste0(project, '/Figures/raw/figure_S2.raw.Rdata'))

write.csv(metrics.kimura, file = paste0(project, '/Figures/tmp/figure_S2.csv'))


#######
# Extract numbers
#

metrics.sum %>% group_by(Setting) %>% summarise(mean(mean))

#
####


## F-Scores - Order --------------------------------------------------------


metrics.order <- df.all %>% filter(order %in% TEs) %>% 
  count(Condition, Tool, sample, Simulation, Setting, order, name = 'count') %>%  
  spread(Condition, count) %>% 
  mutate(Sensitivity = TP/(TP+FN),
         Precision = TP/(TP+FP),
         F_score = 2*(Precision*Sensitivity)/(Precision+Sensitivity),
         FDR = FP/(TP+FP),
         Setting = as.factor(Setting))

metrics.order.sum <- metrics.order %>% 
  group_by(Tool, Setting, order) %>%
  summarise(mean = mean(F_score), 
            sd = sd(F_score),
            se = sd/sqrt(n())) %>% 
  reOrderSettingRev()

metrics.order.sum$Tool <- ordered(metrics.order.sum$Tool, levels = rev(Tools))

panel1.c <- ggplot(metrics.order.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge', 
           aes(group=Setting), 
           color = 'black') +
  facet_grid(.~order) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2) +
  scale_y_continuous(limits = c(0,1.1), expand = c(0,0))+
  my_theme() +
  scale_fill_manual(name = 'Setting', values = my_colors) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5))+
  labs(y = 'mean(F-score)',
       title = 'order') +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 11, family = 'CM Sans', hjust = 1,  vjust = -1.5),
        panel.spacing.x = unit(0.4, "lines"))

save(panel1.c, file = paste0(project, '/Figures/raw/figure_1c.raw.Rdata'))

png(paste0(project, '/Figures/figure_1c.png'), width = 8, height = 5, units = 'in', res = 300)
print(panel1.c)
dev.off()

write.csv(metrics.order.sum, file = paste0(project, '/Figures/tmp/figure_1c.csv'))

## F-Scores - small Kimura -------------------------------------------------


metrics.young <- df.all %>% filter(Kim.grp == '[0,5)') %>%  
  count(Condition, Tool, sample, Simulation, Setting, name = 'count') %>%  
  spread(Condition, count) %>% 
  mutate(Sensitivity = TP/(TP+FN),
         Precision = TP/(TP+FP),
         F_score = 2*(Precision*Sensitivity)/(Precision+Sensitivity),
         FDR = FP/(TP+FP),
         Setting = as.factor(Setting))

metrics.young.sum <- metrics.young %>% 
  group_by(Tool, Setting) %>%
  summarise(mean = mean(F_score), 
            sd = sd(F_score),
            se = sd/sqrt(n()))

metrics.young.sum$Tool <- ordered(metrics.young.sum$Tool, levels = rev(Tools))

#######
# Extract number
#

metrics.young.sum %>% group_by(Setting) %>% summarise(median(mean))

#
####


panel1.b <- ggplot(metrics.young.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge', 
           aes(group=Setting), 
           color = 'black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2) +
  scale_y_continuous(limits = c(0,1.05), expand = c(0,0))+
  my_theme() +
  coord_flip() +
  scale_fill_manual(name = 'Setup',
                    values = c('#335C67','#E09F3E'),
                    guide = guide_legend(reverse = TRUE))+
  labs(y = 'mean(F-score)',
       title = 'Kimura distance < 5') +
  theme(legend.justification=c(0.8,0.8), 
        legend.position=c(0.9,0.15),
        legend.title = element_text(size = 10, family = "CM Sans", face = 'bold'), 
        legend.text = element_text(size = 8, family = "CM Sans"),
        axis.title = element_blank(),
        plot.title = element_text(size = 11, family = 'CM Sans', hjust = 1, vjust = -1.5),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_blank(),
        legend.key.size = unit(0.6,"line")) 


save(panel1.b, file = paste0(project, '/Figures/raw/figure_1b.raw.Rdata'))

png(paste0(project, '/Figures/figure_1b.png'), width = 8, height = 5, units = 'in', res = 300)
print(panel1.b)
dev.off()

write.csv(metrics.young.sum, file = paste0(project, '/Figures/tmp/figure_1b.csv'))

## F-Score - small Kimura and order ---------------------------------------


#remove DNA transpons since there are no DNA transposons in the young group
metrics.young.order <- df.all %>% filter(Kim.grp == '[0,5)', order %in% TEs[2:4]) %>%  
  count(Condition, Tool, sample, Simulation, Setting, order, name = 'count') %>%  
  spread(Condition, count) %>% 
  mutate(Sensitivity = TP/(TP+FN),
         Precision = TP/(TP+FP),
         F_score = 2*(Precision*Sensitivity)/(Precision+Sensitivity),
         FDR = FP/(TP+FP),
         Setting = as.factor(Setting))

metrics.young.order.sum <- metrics.young.order %>% 
  group_by(Tool, Setting, order) %>%
  summarise(mean = mean(F_score), 
            sd = sd(F_score),
            se = sd/sqrt(n())) %>% 
  reOrderSettingRev()

metrics.young.order.sum$Tool <- ordered(metrics.young.order.sum$Tool, levels = rev(Tools))

metrics.young.order.sum %>% filter(Setting == 'single', order == 'LINE')

# How is the improvment of the F-Score in median across all tools for LINE?
#
metrics.young.order.sum %>% 
  filter(order == 'LINE') %>% 
  group_by(Setting) %>% 
  summarise(med = median(mean))

####

panel1.d <- ggplot(metrics.young.order.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge', 
           aes(group=Setting), 
           color = 'black') +
  facet_grid(.~order) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2) +
  scale_y_continuous(limits = c(0,1.05), expand = c(0,0))+
  my_theme() +
  scale_fill_manual(name = 'Setting', values = my_colors) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5))+
  labs(y = 'mean(F-score)',
       title = 'order & Kimura distance < 5') +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 11, family = 'CM Sans', hjust = 1, vjust = -1.5),
        panel.spacing.x = unit(0.4, "lines"))


# Store all relevant data
save(panel1.d, file = paste0(project, '/Figures/raw/figure_1d.raw.Rdata'))

png(paste0(project, '/Figures/figure_1d.png'), width = 8, height = 5, units = 'in', res = 300)
print(panel1.d)
dev.off()

write.csv(metrics.young.order.sum, file = paste0(project, '/Figures/tmp/figure_1d.csv'))


# Quantification of detected TEs ------------------------------------------


## Deviation of recovered counts -------------------------------------------


# How is the difference of the simulated counts compared to the recalled counts
# by the tools?
# I set for all instances with a mean of 0 the log is also set to 0.

total.deviation <- df.all %>% filter(Simulation == 'with DETEs')%>%
  group_by(Tool, TE, Setting) %>% 
  summarise(mean.recov = round(mean(recovered.counts),  digits = 0),
            mean.sim = round(mean(simulated.counts), digits = 0),
            mean.recov.log = case_when(mean.recov > 0 ~ log(mean.recov),
                                       TRUE ~ 0),
            mean.sim.log = case_when(mean.sim > 0 ~ log(mean.sim),
                                     TRUE ~ 0),
            .groups = 'drop') %>% 
  reOrderSetting() %>% 
  reOrderTools()

save(total.deviation, file = paste0(project, '/Data/total.deviation.set2.raw.Rdata'))

panel2.a <- ggplot(total.deviation, aes(mean.sim.log, mean.recov.log)) +
  geom_point(alpha = 0.025, color = 'grey') +
  geom_density2d(linetype=1, size= 0.25, bins=30, aes(color = Tool)) +
  geom_abline(lty=2, color='black', size=1) +
  geom_quadrant_lines(yintercept = log(5), xintercept = log(5)) +
  stat_quadrant_counts(xintercept = log(5), 
                       yintercept = log(5), 
                       family = "CM Sans", 
                       size = 2.5,
                       geom = "label_npc") +
  facet_grid(Tool~Setting)+
  my_theme() +
  scale_color_manual(name = 'Tool', values = my_colors) +
  labs(x = expression(paste("log(", baseMean[Simulated], ")")),
       y = expression(paste("log(", baseMean[Recovered], ")"))) +
  theme(legend.position = "none")

# add R^2 of the TP to the plot
# R^2 Calculation for true positives
# Do not need to filter for 'with DETEs' because it is already happened for the
# total.deviation data frame. First I calculate th R squared value by hand and
# add them afterwards to the plot.

r.squared <- total.deviation %>% 
  filter(mean.recov.log >= log(5), mean.sim.log >= log(5)) %>% 
  group_by(Tool, Setting) %>% 
  summarise(r.squared = summary(lm(mean.sim.log~mean.recov.log))$r.squared)

panel2.a <- 
  panel2.a + 
  geom_text(data = r.squared,
            aes(x=7.5, y=2.5, label = paste0('R^2 == ', round(r.squared,3))), 
            show.legend = F,
            parse = T,
            size = 4,
            family = "CM Sans")

save(panel2.a, file = paste0(project, '/Figures/raw/figure_2.raw.Rdata'))

png(paste0(project, "/Figures/figure_2.png"), width = 12.5, height = 15.5, units = 'in', res = 600)
print(panel2.a) # Make plot
dev.off()

write.csv(total.deviation, file = paste0(project, '/Figures/tmp/figure_2.csv'))
write.csv(r.squared, file = paste0(project, '/Figures/tmp/figure_2_Rsquared.csv'))

### Supplemental Figure for Set 1

total.deviation.set1 <- df.all %>% filter(Simulation == 'ETEs')%>%
  group_by(Tool, TE, Setting) %>% 
  summarise(mean.recov = round(mean(recovered.counts),  digits = 0),
            mean.sim = round(mean(simulated.counts), digits = 0),
            mean.recov.log = case_when(mean.recov > 0 ~ log(mean.recov),
                                       TRUE ~ 0),
            mean.sim.log = case_when(mean.sim > 0 ~ log(mean.sim),
                                     TRUE ~ 0),
            .groups = 'drop') %>% 
  reOrderSetting() %>% 
  reOrderTools()

save(total.deviation.set1, file = paste0(project, '/Data/total.deviation.set1.raw.Rdata'))

sup.fig.3 <- ggplot(total.deviation.set1, aes(mean.sim.log, mean.recov.log)) +
  geom_point(alpha = 0.025, color = 'grey') +
  geom_density2d(linetype=1, size= 0.25, bins=30, aes(color = Tool)) +
  geom_abline(lty=2, color='black', size=1) +
  geom_quadrant_lines(yintercept = log(5), xintercept = log(5)) +
  stat_quadrant_counts(xintercept = log(5), 
                       yintercept = log(5), 
                       family = "CM Sans", 
                       size = 2.5,
                       geom = "label_npc") +
  facet_grid(Tool~Setting)+
  my_theme() +
  scale_color_manual(name = 'Tool', values = my_colors) +
  labs(x = expression(paste("log(", baseMean[Simulated], ")")),
       y = expression(paste("log(", baseMean[Recovered], ")"))) +
  theme(legend.position = "none")



r.squared.set1 <- total.deviation.set1 %>% 
  filter(mean.recov.log >= log(5), mean.sim.log >= log(5)) %>% 
  group_by(Tool, Setting) %>% 
  summarise(r.squared = summary(lm(mean.sim.log~mean.recov.log))$r.squared)

sup.fig.3  <- 
  sup.fig.3 + 
  geom_text(data = r.squared.set1,
            aes(x=6.5, y=2.5, label = paste0('R^2 == ', round(r.squared,3))), 
            show.legend = F,
            parse = T,
            size = 4,
            family = "CM Sans")

save(sup.fig.3, file = paste0(project, '/Figures/raw/figure_S3.raw.Rdata'))

png(paste0(project, "/Figures/figure_S3.png"), width = 12.5, height = 15.5, units = 'in', res = 600)
print(sup.fig.4) # Make plot
dev.off()

write.csv(total.deviation.set1, file = paste0(project, '/Figures/tmp/figure_S3.csv'))
write.csv(r.squared.set1, file = paste0(project, '/Figures/tmp/figure_S3_Rsquared.csv'))

####
# Extract numbers 

r.squared %>% group_by(Setting) %>% summarise(med = median(r.squared))


## Proportion of different condition of all detected TEs -------------------


binary.df <- total.deviation %>% 
  mutate(binary.class = case_when(mean.recov.log >= log(5) & mean.sim.log >= log(5) ~ "TP",
                                  mean.recov.log <= log(5) & mean.sim.log >= log(5) ~ "FN",
                                  mean.recov.log >= log(5) & mean.sim.log <= log(5) ~ "FP",
                                  TRUE ~NA_character_)) %>% 
  filter(!is.na(binary.class)) %>% 
  dplyr::count(Tool, Setting, binary.class) %>% 
  group_by(Tool, Setting) %>% 
  mutate(percent = n/sum(n)) %>% 
  reOrderSetting() %>% 
  reOrderTools()

binary.col <-  c('#557174', '#3b6978','#84a9ac' )
names(binary.col) <- c('TP', 'FP', 'FN')       

panel2b <- ggplot(binary.df, aes(x = Setting, y = percent, fill = binary.class)) +
  geom_col() +
  facet_grid(.~Tool) +
  my_theme() +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.spacing.x=unit(0.0, "lines"),
        strip.background = element_blank()) +
  scale_fill_manual(name= 'Condition', values = binary.col) +
  labs(y = 'proportion of conditions') +
  theme(axis.title.x = element_blank(),
        title = element_text(family = "CM Sans", size = 12))

### Store raw plot
save(panel2b, file = paste0(project, '/Figures/panel2b.raw.Rdata'))

# Evaluate detection of DETEs ---------------------------------------------

## recall in contrast to fdr -----------------------------------------------

### Load Data
#
# I can take the ground truth of one setting because the ground truth is
# equal for both settings.
#

deseq.results.blank <- data.frame()

for(setting in c('single', 'paired')){
  
  df <- loadRdata(paste0(project, '/Data/', setting, '.combined.deseq2.res.processed.Rdata'))
  
  df$setting = setting
  deseq.results.blank <- rbind(deseq.results.blank, df)
  
  rm(df)
}

ground.truth <- loadRdata(paste0(project, '/Data/paired.trueDETEs.Rdata'))

deseq.results.blank <- deseq.results.blank %>%
  mutate(
    Condition = case_when(
      padj < 0.05 & TE %in% ground.truth ~ 'TP',
      padj < 0.05 & !(TE %in% ground.truth) ~ 'FP',
      padj > 0.05 & TE %in% ground.truth ~ 'FN',
      padj > 0.05 & !(TE %in% ground.truth) ~ 'TN',
      TRUE ~ NA_character_
    ),
    diff = case_when(TE %in% ground.truth ~ 1,
                     TRUE ~ 0)
  ) 


### for all -----------------------------------------------------------------

deseq.results <- deseq.results.blank %>% 
  dplyr::group_by(Tool, setting) %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    total = 1:n(),
    truepositive = base::cumsum(diff),
    falsepositive = total - truepositive,
    fdr = falsepositive / total,
    precision = truepositive / total,
    recall = truepositive / length(ground.truth))

deseq.results$setting <- ordered(deseq.results$setting, levels=c('single', 'paired'))



recall.fdr.all <- 
  ggplot(deseq.results, aes(x = recall, y = fdr, color = Tool)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.1,
             linetype = "dashed",
             color = "black") +
  scale_color_manual(name = 'Tool', values = my_colors) +
  coord_cartesian(ylim=c(0,1), expand = FALSE ) +
  scale_y_continuous(position = "right") +
  facet_grid(.~setting) +
  theme_bw() +
  theme(legend.justification=c(0,1), 
        legend.position=c(0.02,0.98),
        panel.grid.minor = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.background = element_blank(),
        strip.background=element_rect(colour="grey", fill="white"),
        strip.text = element_text(family = "CM Sans", size = 12),
        axis.text = element_text(family = "CM Sans", size = 10),
        axis.title = element_text(family = "CM Sans", size = 12),
        title = element_text(family = "CM Sans", size = 12),
        legend.key.size = unit(0.4,"line"),
        legend.title = element_text(size = 9, family = "CM Sans", face = 'bold'), 
        legend.text = element_text(size = 7, family = "CM Sans"))

### Store raw plot

save(recall.fdr.all, file = paste0(project, '/Figures/raw/figure_3a.raw.Rdata'))

png(paste0(project, "/Figures/figure_3a.png"), width = 8, height = 5, units = 'in', res = 600)
print(recall.fdr.all) # Make plot
dev.off()

write.csv(deseq.results, file = paste0(project, '/Figures/tmp/figure_3a.csv'))

### for TEs with a kimura smaller than 5 ------------------------------------

ground.truth.by.age <- data.frame(ground.truth) %>% 
  splitTEID('ground.truth') %>% 
  dplyr::mutate(Kimura = as.numeric(Kimura),
                Kim.grp = case_when(
                  Kimura < 5 ~ '[0,5)',
                  TRUE ~ '> 5'),
                Kim.grp = factor(
                  Kim.grp,
                  levels = c(
                    '[0,5)',
                    '> 5'
                  ))) %>% count(Kim.grp)


deseq.results.young <- deseq.results.blank %>% splitTEID('TE') %>% 
  dplyr::mutate(Kimura = as.numeric(Kimura),
                Kim.grp = case_when(
                  Kimura < 5 ~ '[0,5)',
                  TRUE ~ '> 5'),
                Kim.grp = factor(
                  Kim.grp,
                  levels = c(
                    '[0,5)',
                    '> 5'
                  ))) %>% 
  dplyr::group_by(Tool, setting, Kim.grp) %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    total = 1:n(),
    truepositive = base::cumsum(diff),
    falsepositive = total - truepositive,
    fdr = falsepositive / total,
    precision = truepositive / total,
    recall = case_when(Kim.grp == '[0,5)' ~ (truepositive / 434),
                       Kim.grp == '> 5' ~ (truepositive / 4028),
                       TRUE ~ NA_real_ )
  )


deseq.results.young$setting <- ordered(deseq.results.young$setting, levels=c('single', 'paired'))

recall.fdr.age <- 
  ggplot(deseq.results.young %>% filter(Kim.grp == '[0,5)'), aes(x = recall, y = fdr, color = Tool)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.1,
             linetype = "dashed",
             color = "grey80") +
  scale_color_manual(name = 'Tool', values = my_colors) +
  coord_cartesian(ylim=c(0,1), expand = FALSE ) +
  facet_grid(.~setting) +
  scale_y_continuous(position = "right") +
  theme_bw() +
  theme(legend.position='none',
        panel.grid.minor = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.background = element_blank(),
        strip.background=element_rect(colour="grey", fill="white"),
        strip.text = element_text(family = "CM Sans", size = 12),
        axis.text = element_text(family = "CM Sans", size = 10),
        axis.title = element_text(family = "CM Sans", size = 12),
        title = element_text(family = "CM Sans", size = 12))

######
# Extract some numbers

deseq.results.young %>% filter(Kim.grp == '[0,5)', between(fdr, 0.11, 0.09), setting == 'paired') %>% View()

deseq.results.young %>% filter(Kim.grp == '[0,5)', Tool == 'TEtranscripts', setting == 'paired') %>% View()

######
# Extract some numbers

deseq.results %>% filter(fdr == 0.1, setting == 'single') %>% View()

deseq.results.young %>% filter(Kim.grp == '[0,5)', Tool == 'TEtranscripts', setting == 'paired') %>% View()

### Store raw plot
save(recall.fdr.age, file = paste0(project, '/Figures/raw/figure_3b.Rdata'))

png(paste0(project, "/Figures/figure_3b.png"), width = 8, height = 5, units = 'in', res = 600)
print(recall.fdr.age) # Make plot
dev.off()

write.csv(deseq.results.young, file = paste0(project, '/Figures/tmp/figure_3b.csv'))

## log2Fold Change of FPs --------------------------------------------------

## Since false positives are responsible for the fdr rates in the plot before
## I want to check how is the expression and log2FC of such false positives?
## Therefore I filtered the data frame only for FP and take only that TEs
## with a median bigger than 0. The median ensures to take only instances that
## are detected in at least the half of the replicates. For the log2FC I 
## devide the mean expression of condition 1 by the mean expression of 
## condition 2. I added 1 to the mean before, to avoid the case of division by
## zero. The diagonal arrangement of dots are TEs that where only in one 
## condition detected. 

df.fp <- df.all %>% 
  filter(Condition == 'FP') %>% 
  dplyr::select(TE, sample, recovered.counts, Tool, Setting) %>% 
  spread(key = sample, value = recovered.counts) %>% # rearrange the count table
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% # set all NAs to 0
  rowwise() %>% 
  mutate(medianCount = median(c(sample_1,sample_2,sample_3,sample_4,sample_5,
                                sample_diff_1,sample_diff_2,sample_diff_3,sample_diff_4,sample_diff_5)),
         meanCount = mean(c(sample_1,sample_2,sample_3,sample_4,sample_5,
                                sample_diff_1,sample_diff_2,sample_diff_3,sample_diff_4,sample_diff_5))) %>% 
  filter(medianCount != 0) %>% 
  mutate(meanCondition1 = mean(c(sample_1,sample_2,sample_3,sample_4,sample_5) +1),
         meanCondition2 = mean(c(sample_diff_1,sample_diff_2,sample_diff_3,sample_diff_4,sample_diff_5)+1)) %>% 
  mutate(log2Fold = log(meanCondition2/meanCondition1))

df.fp$Tool <- factor(df.fp$Tool, levels=c('salmonTE', 'Telescope', 'TEtranscripts', 'SQuIRE', 'TEtools')) 
df.fp$Setting <- as.factor(df.fp$Setting)

df.fp$Setting <- factor(df.fp$Setting, levels=c('single', 'paired'))

log2FoldFP <- df.fp %>% filter(abs(log2Fold)>0.5) %>% 
  ggplot(aes(x=log2Fold, y=log(meanCount))) +
  geom_point(aes(color = log2Fold), alpha=0.1, size=0.5) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=1, bins=50, geom='polygon') +
  scale_fill_gradient(low = 'yellow', high = 'red') +
  scale_alpha(range = c(0.00, 1.0), guide = FALSE) +
  geom_density2d(colour='black', bins=30, size=0.25) +
  labs(x = expression(paste(log[2], "(fold change)")),
       y = expression(paste(log[2], '(mean count)'))) +
  guides(alpha=FALSE) +
  facet_grid(Tool~Setting) + 
  theme_bw() + 
  theme(legend.position = 'none')+
  theme(panel.grid.minor = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background=element_rect(colour="grey", fill="white"),
        strip.text = element_text(family = "CM Sans", size = 12),
        axis.text = element_text(family = "CM Sans", size = 10),
        axis.title = element_text(family = "CM Sans", size = 12),
        title = element_text(family = "CM Sans", size = 12))


### Store raw plot
save(log2FoldFP, file = paste0(project, '/Figures/raw/figure_3c.Rdata'))

png(paste0(project, "/Figures/figure_3c.png"), width = 12.5, height = 10, units = 'in', res = 600)
print(log2FoldFP) # Make plot
dev.off()

write.csv(df.fp %>% filter(abs(log2Fold)>0.5), 
          file = paste0(project, '/Figures/tmp/figure_3c.csv'))


## Supplement Age of Young elements


simdf <- loadRdata(paste0(project, '/Data/Simulation.paired.cnttbl.processed.Rdata'))


simdf <- simdf %>% splitTEID('TE')

simdf.young <- simdf %>% 
  filter(as.numeric(Kimura) < 5, order %in% c('SINE', 'LTR', 'LINE')) %>% 
  mutate(Kimura = as.numeric(Kimura),
         Kimura = round(Kimura, digits = 1))


pl.supplement.2 <- ggplot(simdf.young, aes(Kimura)) +
  geom_bar(width = 0.1, color = 'darkgrey', fill = 'grey')+
  facet_grid(.~order) +
  scale_y_continuous(limits = c(0,155), expand = c(0,0)) +
  theme_bw() +
  labs(x = 'Kimura distance',
       y = 'TE instance count') +
  theme(panel.grid.major = element_blank(),
        legend.position = 'none',
        panel.spacing.x = unit(1, 'lines'),
        panel.spacing.y = unit(1, 'lines'),
        strip.background=element_rect(colour="grey", fill="white"),
        axis.text = element_text(family = "CM Sans", size = 16),
        axis.title = element_text(family = "CM Sans", size = 18),
        strip.text = element_text(family = "CM Sans", size = 18))
  
save(pl.supplement.2, file = paste0(project, '/Figures/raw/figure_S1.Rdata'))

png(paste0(project, "/Figures/figure_S1.png"), width = 8, height = 7, units = 'in', res = 600)
print(pl.supplement.2) # Make plot
dev.off()



write.csv(simdf.young, 
          file = paste0(project, '/Figures/tmp/figure_S1.csv'))

  
