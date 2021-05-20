library(ggplot2)
library(tidyverse)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggpmisc) # package to count measurements in quadrants
library(ggtext) # to add r² in a nice form
library(glue) # glue text
#library(extrafont)
library(showtext) # to load the Source Sans Pro font
library(cowplot) #arrange plots
source('libs/func.R')
source('general.R')



# global Variables --------------------------------------------------------

## Load font from google and set it to default

#font_add_google(name = "Dancing Script", family = "Dancing Script")
font_add_google(name = "Source Sans Pro", family = "Source Sans Pro")
showtext_auto()

## TEs and tools of interest

TEs <- c('DNA', 'LINE', 'LTR', 'SINE')

Tools <- c('SalmonTE', 'Telescope', 'TEtranscripts', 'SQuIRE', 'TEtools')

project <- project.name

## Colors

#my_colors <- c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F", '#af7ac5', '#cbcbcb', '#909090')
my_colors <- c("#d95b35ff", "#0d0d49ff", "#64a183ff", "#e69c23ff", '#9f5cbaff', '#cbcbcb', '#909090')
names(my_colors) <- c('SalmonTE', 'SQuIRE', 'TEtools', 'TEtranscripts', 'Telescope', 'Single-end setup', 'Paired-end setup')

text.color <- c("#64a183ff", "#0d0d49ff",  "#e69c23ff", '#9f5cbaff', "#d95b35ff")#"#E07A5F")

# script relevant functions -----------------------------------------------

reOrderSetting <- function(df){
  
  df <- df %>% mutate(Setting = factor(Setting, levels = c('Single-end setup', 'Paired-end setup')))
  return(df)
}

reOrderSettingRev <- function(df){
  df <- df %>% mutate(Setting = factor(Setting, levels = c( 'Paired-end setup','Single-end setup')))
  return(df)
}

reOrderTools <- function(df){
  
  Tools <- c('SalmonTE', 'Telescope', 'TEtranscripts', 'SQuIRE', 'TEtools')
  df <- df %>% mutate(Tool = factor(Tool, levels = Tools))
  return(df)
  
}


my_theme <- function(){
  
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white", color = "grey80"),
        panel.spacing.x = unit(0.2, "lines"),
        panel.border =  element_rect(fill = "NA", color = "grey80"),
        legend.key         = element_rect(fill = "white"))
}

color_strips <- function(pl){
  
  cols <- c('#cbcbcb', '#909090', "#E07A5F", "#af7ac5", "#F2CC8F", "#3D405B", '#81B29A')
  txtcols <- c('#000000', '#000000', "#000000", "#000000", "#000000", "#ffffff", '#000000')
  g <- ggplot_gtable(ggplot_build(pl))
  
  strip_both <- which(grepl('strip-', g$layout$name))
  
  k <- 1
  
  for(i in seq_along(strip_both)){
    
    j <- which(grepl('rect', g$grobs[[strip_both[i]]]$grobs[[1]]$childrenOrder))
    l <- which(grepl('titleGrob', g$grobs[[strip_both[i]]]$grobs[[1]]$childrenOrder))
    
    g$grobs[[strip_both[i]]]$grobs[[1]]$children[[j]]$gp$fill <- cols[i]
    g$grobs[[strip_both[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- txtcols[i]
    k <- k+1
  }
  
  return(g)
}

shade_strips <- function(pl){
  
  cols <- c('#cbcbcb', '#909090')
  
  g <- ggplot_gtable(ggplot_build(pl))
  
  strip_both <- which(grepl('strip-', g$layout$name))
  
  k <- 1
  
  for(i in strip_both){
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- cols[k]
    k <- k+1
  }
  
  return(g)
}

fScore_theme <- function(pl, legend = FALSE){
  
  defW <- getOption("warn")
  options(warn=-1)
  pl <- pl + 
    scale_y_continuous(limits = c(0,1), 
                       expand = c(0,0),
                       breaks = c(0,0.25,0.5,0.75,1),
                       labels = c(0, "", 0.5, "",1)) +
    coord_flip() +
    my_theme() +
    scale_fill_manual(name = 'Setting', values = my_colors) +
    theme(legend.position = 'bottom',
          axis.title = element_blank(),
          plot.title = element_text(size = 12, family = "Source Sans Pro"),
          axis.text.y = element_text(colour = text.color, size = 10),
          panel.spacing.x = unit(0.5, "lines"),
          strip.background = element_rect(fill = 'white', colour = 'grey'),
          strip.text = element_text(family = "Source Sans Pro", size = 10),
          axis.text.x = element_text(family = "Source Sans Pro", size = 10),
          title = element_text(family = "Source Sans Pro", size = 12))
  
  
  if(!legend){
    pl <- pl + theme(legend.position = 'none')
  }
  
  options(warn = defW)
  return(pl)
}

dete_supplement_theme <- function(pl){
  
  pl <- pl +  
    geom_hline(yintercept = 0.1,
               linetype = "dashed",
               color = "black") +
    scale_color_manual(name = 'Tool', values = my_colors) +
    coord_cartesian(ylim=c(0,1), expand = FALSE ) +
    scale_x_continuous(breaks = c(0, 0.25,0.5, 0.75, 1), 
                       labels = c(0,"",0.5, "",1)) +
    scale_y_continuous(breaks = c(0, 0.25,0.5, 0.75, 1), 
                       labels = c(0,"",0.5, "",1)) +
    theme_bw() +
    theme(legend.position='none',
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'grey', linetype = 'dotted'),
          panel.spacing.x = unit(0.9, "lines"),
          panel.background = element_blank(),
          strip.background=element_rect(colour="grey", fill="white"),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 8),# face = 'bold'),
          axis.title = element_text(size = 10),# face = 'bold'),
          title = element_text(size = 10, color = '#414141'),
          legend.key.size = unit(0.4,"line"),
          legend.title = element_text(size = 9, face = 'bold'), 
          legend.text = element_text(size = 7 ))
  
  return(pl)

}

dete_main_theme <- function(pl){
  
  pl <- pl +  
    geom_line(size = 0.7) +
    geom_hline(yintercept = 0.1,
               linetype = "dashed",
               color = "black") +
    scale_color_manual(name = 'Tool', values = my_colors) +
    coord_cartesian(ylim=c(0,1), expand = FALSE ) +
    scale_x_continuous(breaks = c(0, 0.25,0.5, 0.75, 1), 
                       labels = c(0,"",0.5, "",1)) +
    scale_y_continuous(breaks = c(0, 0.25,0.5, 0.75, 1), 
                       labels = c(0,"",0.5, "",1)) +
    theme_bw() +
    theme(legend.position='none',
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'grey', linetype = 'dotted'),
          panel.spacing.x = unit(0.9, "lines"),
          panel.background = element_blank(),
          strip.background=element_rect(colour="grey", fill="white"),
          strip.text = element_text(size = 14),
          axis.text = element_text(size = 10, face = 'bold'),
          axis.title = element_text(size = 12, face = 'bold'),
          title = element_text(size = 12, color = '#414141'),
          legend.key.size = unit(0.4,"line"),
          legend.title = element_text(size = 9, face = 'bold'), 
          legend.text = element_text(size = 7 ))
  
  return(pl)
}


my_colors <- c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F", '#af7ac5', '#cbcbcb', '#909090')
names(my_colors) <- c('salmonTE', 'SQuIRE', 'TEtools', 'TEtranscripts', 'Telescope', 'single', 'paired')

text.color <- c("#81B29A", "#3D405B",  "#F2CC8F", '#af7ac5',"#E07A5F")
# load general data -------------------------------------------------------

# The following data frames contain the raw counts of the tools. The count tables
# are processed, which means that detected TEs with less than 5 reads in sum
# across all samples are removed. Furthermore the count tables are rearranged
# so that per detected TE per replicate, Tool and Setup a row exists where the
# recovered and the simulated counts are stored and beside these a lot of 
# additional information are added. These information is used to calculate
# binary measurement of TP, FP and FN, which is used to calculate the F-score.

df.single <- loadRdata(paste0(project, '/Data/single.combined.cnttbl.processed.Rdata'))
df.paired <- loadRdata(paste0(project, '/Data/paired.combined.cnttbl.processed.Rdata'))

df.all <- rbind(df.single, df.paired)

df.all <- df.all %>% 
  mutate(Tool = case_when(as.character(Tool) == 'salmonTE' ~ 'SalmonTE',
                                             TRUE ~ Tool),
         Setting = case_when(as.character(Setting) == 'single' ~ 'Single-end setup',
                             as.character(Setting) == 'paired' ~ 'Paired-end setup'))

# Main figures ------------------------------------------------------------

# Panel 1 -----------------------------------------------------------------

## F-Scores ----------------------------------------------------------------
#
# The F-Scores represent a trade-off between precision and recall. The higher
# the F-Score, the better the tool performs. The F-Score is calculated for
# different groups - all, small kimura, order, and order & small kimura.


### A) F-Scores - All ----------------------------------------------------------

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
  reOrderSetting() %>% 
  mutate(pseudoFacette = 'pseudo')

metrics.sum$Setting <- ordered(metrics.sum$Setting, levels=c('Paired-end setup', 'Single-end setup'))
metrics.sum$Tool <- ordered(metrics.sum$Tool, levels = rev(Tools))

panel1.a <- ggplot(metrics.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge',
           aes(group=Setting)) +

  facet_grid(.~pseudoFacette) +

  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2)  +
  labs(title = 'All')
  
panel1.a <- fScore_theme(panel1.a, legend = TRUE)

legend <- cowplot::get_legend(panel1.a)

panel1.a <- panel1.a +
  theme(legend.position = 'none',
        strip.text = element_text(color = 'white', family = "Source Sans Pro", size = 16, face = 'bold'),
        strip.background = element_rect(fill = 'white', colour = 'white'))

save(panel1.a, file = paste0(project, '/Figures/raw/figure_1a.raw.Rdata'))

png(paste0(project, '/Figures/figure_1a.png'), width = 550, height = 300, units = 'px', family = 'Source Sans Pro')
print(panel1.a)
dev.off()

write.csv(metrics.sum, file = paste0(project, '/Figures/tmp/figure_1a.csv'))


ggsave(filename = paste0(project, "/Figures/figure_1a.pdf"), 
       print(panel1.a),
       width = 3, height = 2, dpi = 600, units = "in", device=cairo_pdf)

# Extract numbers
#

metrics.sum %>% group_by(Setting) %>% summarise(mean(mean))

#
####


### C) F-Scores - Order --------------------------------------------------------


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
           aes(group=Setting)) +
  facet_grid(.~order) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2)+
  labs(title = 'By TE orders')

panel1.c <- fScore_theme(panel1.c)

save(panel1.c, file = paste0(project, '/Figures/raw/figure_1c.raw.Rdata'))

png(paste0(project, '/Figures/figure_1c.png'), width = 550, height = 300, units = 'px', family = 'Source Sans Pro')
print(panel1.c)
dev.off()

write.csv(metrics.order.sum, file = paste0(project, '/Figures/tmp/figure_1c.csv'))


ggsave(filename = paste0(project, "/Figures/figure_1c.pdf"), 
       print(panel1.c),
       width = 4.5, height = 2, dpi = 600, units = "in", device=cairo_pdf)


### B) F-Scores - small Kimura -------------------------------------------------


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
            se = sd/sqrt(n()))  %>% 
  mutate(pseudoFacette = 'pseudo')

metrics.young.sum$Tool <- ordered(metrics.young.sum$Tool, levels = rev(Tools))


# Extract number
#

metrics.young.sum %>% group_by(Setting) %>% summarise(median(mean))

#
####


panel1.b <- ggplot(metrics.young.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge', 
           aes(group=Setting)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2) +
  facet_grid(.~pseudoFacette) +
  labs(title = 'Kimura distance < 5')

panel1.b <- fScore_theme(panel1.b)

panel1.b <- panel1.b +
  theme(strip.text = element_text(color = 'white'),
        strip.background = element_rect(fill = 'white', colour = 'white'))

save(panel1.b, file = paste0(project, '/Figures/raw/figure_1b.raw.Rdata'))

png(paste0(project, '/Figures/figure_1b.png'), width = 550, height = 300, units = 'px', family = 'Source Sans Pro')
print(panel1.b)
dev.off()

write.csv(metrics.young.sum, file = paste0(project, '/Figures/tmp/figure_1b.csv'))


ggsave(filename = paste0(project, "/Figures/figure_1b.pdf"), 
       print(panel1.b),
       width = 3, height = 2, dpi = 600, units = "in", device=cairo_pdf)


### D) F-Score - small Kimura and order ---------------------------------------

metrics.young.order <- df.all %>% filter(Kim.grp == '[0,5)', order %in% TEs) %>%  
  count(Condition, Tool, sample, Simulation, Setting, order, name = 'count') %>%  
  spread(Condition, count) %>% 
  mutate_if(is.integer, ~replace(.,is.na(.), 0))  %>% 
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

# How is the improvement of the F-Score in median across all tools for LINE?

metrics.young.order.sum %>% 
  filter(order == 'LINE') %>% 
  group_by(Setting) %>% 
  summarise(med = median(mean))


panel1.d <- ggplot(metrics.young.order.sum, aes(x = Tool, mean, fill=Setting)) +
  geom_col(width=0.8, 
           position = 'dodge',
           aes(group=Setting)) +
  facet_grid(.~order) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(0.8),
                size = .6,
                color = 'black',
                width = .2) +
  labs(title = 'Kimura distance < 5 & by TE orders')


panel1.d <- fScore_theme(panel1.d)

# Store all relevant data
save(panel1.d, file = paste0(project, '/Figures/raw/figure_1d.raw.Rdata'))

png(paste0(project, '/Figures/figure_1d.png'),width = 550, height = 300, units = 'px', family = 'Source Sans Pro')
print(panel1.d)
dev.off()

write.csv(metrics.young.order.sum, file = paste0(project, '/Figures/tmp/figure_1d.csv'))

ggsave(filename = paste0(project, "/Figures/figure_1d.pdf"), 
       print(panel1.d),
       width = 4.5, height = 2, dpi = 600, units = "in", device=cairo_pdf)

# Arrange Panel 1 ---------------------------------------------------------

panel1 <- ggarrange(panel1.a, 
                    panel1.c, 
                    panel1.b, 
                    panel1.d,
                    font.label = (family = "Source Sans Pro"),
                    labels = c('A', 'C', 'B', 'D'))

panel1 <- annotate_figure(panel1,
                left = text_grob('Tool', rot= 90, size =16),
                bottom = text_grob('mean(F-score)', size= 16))

png(paste0(project, '/Figures/panel1.png'), width = 1200 , height = 650, 
    units = 'px',family = "Source Sans Pro")
print(panel1)
dev.off()


# Panel 2 -----------------------------------------------------------------

## Quantification of detected TEs ------------------------------------------

### Deviation of recovered counts -------------------------------------------

# How is the difference of the simulated counts compared to the counts recalled 
# by the tools?
# For all instances with a mean of 0 the log is also set to 0.

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


#total.deviation <- total.deviation %>% sample_n(10000)

# Calculate the correlation factor
# add R^2 of the TP to the plot
# R^2 Calculation for true positives
# Do not need to filter for 'with DETEs' because it is already happened for the
# total.deviation data frame. I calculate th R squared value and add the result
# to the respective plot.

tool_cor <- total.deviation %>% 
  filter(mean.recov.log >= log(5), mean.sim.log >= log(5)) %>% 
  group_by(Tool, Setting) %>% 
  summarise(r.squared = summary(lm(mean.sim.log~mean.recov.log))$r.squared) %>% 
  mutate(mean.sim.log = 9.5, 
         mean.recov.log = 2,
         label = glue("*r*<sup>2</sup>(TP) = {round(r.squared, 2)}"))


# generate plot

panel2 <- ggplot(total.deviation, aes(mean.sim.log, mean.recov.log)) +
  geom_point(alpha = 0.025, color = 'grey') +
  geom_density2d(linetype=1, size= 0.25, bins=30, aes(color = Tool)) +
  geom_abline(lty=2, color='black', size=1) +
  geom_quadrant_lines(yintercept = log(5), xintercept = log(5)) +
  stat_quadrant_counts(xintercept = log(5), 
                       yintercept = log(5), 
                       size = 5.6, #5.6
                       geom = "label_npc",
                       quadrants = c(1,2,4)) +
  facet_grid(Tool~Setting)+
  my_theme() +
  scale_color_manual(name = 'Tool', values = my_colors) +
  labs(x = expression(paste("log(", baseMean[Simulated], ")")),
       y = expression(paste("log(", baseMean[Recovered], ")"))) +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) + 
  scale_y_continuous(breaks = c(0,2,4,6,8),
                     labels = c(0,2,4,6,8)) +
  scale_x_continuous(breaks = c(0,2,4,6,8),
                     labels = c(0,2,4,6,8)) +
  aes(color = Tool) +
  geom_richtext(data = tool_cor,
                aes(label = label, fill = after_scale(alpha(color, .2))),
                text.color = 'black',
                size = 5, #5.6,
                family = 'Source Sans Pro',
                hjust = 1, vjust = 0) 

panel2 <- color_strips(panel2)

save(panel2, file = paste0(project, '/Figures/raw/figure_2.raw.Rdata'))

png(paste0(project, "/Figures/figure_2.png"), width = 718, height = 1046, units = 'px', family = "Source Sans Pro")
grid.draw(panel2) # Make plot
dev.off()

write.csv(total.deviation, file = paste0(project, '/Figures/tmp/figure_2.csv'))
write.csv(tool_cor, file = paste0(project, '/Figures/tmp/figure_2_Rsquared.csv'))


####
# Extract numbers 

tool_cor %>% group_by(Setting) %>% summarise(med = median(r.squared))

# Panel 3 -----------------------------------------------------------------

## Evaluate detection of DETEs - TPR vs FDR--------------------------------

### All Instances -----------------------------------------------------------


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
                     TRUE ~ 0),
    Tool = case_when(Tool == 'salmonTE' ~ 'SalmonTE',
                     TRUE ~ Tool),
    Setting = case_when(as.character(setting) == 'single' ~ 'Single-end setup',
                        as.character(setting) == 'paired' ~ 'Paired-end setup')
  ) 

deseq.results <- deseq.results.blank %>% 
  dplyr::group_by(Tool, Setting) %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    total = 1:n(),
    truepositive = base::cumsum(diff),
    falsepositive = total - truepositive,
    FDR = falsepositive / total,
    precision = truepositive / total,
    TPR = truepositive / length(ground.truth))

deseq.results$Setting <- ordered(deseq.results$Setting, 
                                 levels=c('Single-end setup', 'Paired-end setup'))

panel3.a <- 
  ggplot(deseq.results, aes(x = TPR, y = FDR, color = Tool)) +
  geom_line(size = 0.5) +
  facet_grid(.~Setting)


if(project != 'polyester_mm_100bp_fG'){

  panel3.a <- dete_supplement_theme(panel3.a)

}else{
  
  panel3.a <- dete_main_theme(panel3.a)
}



### Store raw plot

save(panel3.a, file = paste0(project, '/Figures/raw/figure_3a.raw.Rdata'))



panel3.a <- shade_strips(panel3.a)

png(paste0(project, "/Figures/figure_3a.png"), 
    width = 4.1, height = 2.4, units = 'in', res = 150)
grid.draw(panel3.a) # Make plot
dev.off()

write.csv(deseq.results, file = paste0(project, '/Figures/tmp/figure_3a.csv'))

### for TEs with a kimura smaller than 5 ------------------------------------

n.young <- as.data.frame(ground.truth) %>% 
  splitTEID('ground.truth') %>% 
  mutate(Kimura = as.numeric(Kimura)) %>% 
  filter(Kimura < 5) %>% nrow()

n.old <- as.data.frame(ground.truth) %>% 
  splitTEID('ground.truth') %>% 
  mutate(Kimura = as.numeric(Kimura)) %>% 
  filter(Kimura >= 5) %>% nrow()

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
  dplyr::group_by(Tool, Setting, Kim.grp) %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    total = 1:n(),
    truepositive = base::cumsum(diff),
    falsepositive = total - truepositive,
    FDR = falsepositive / total,
    precision = truepositive / total,
    TPR = case_when(Kim.grp == '[0,5)' ~ (truepositive / n.young),
                       Kim.grp == '> 5' ~ (truepositive / n.old),
                       TRUE ~ NA_real_ )
  )


deseq.results.young$Setting <- ordered(deseq.results.young$Setting, levels=c('Single-end setup', 'Paired-end setup'))

panel3.b <- 
  ggplot(deseq.results.young %>% filter(Kim.grp == '[0,5)'), aes(x = TPR, y = FDR, color = Tool)) +
  geom_line(size = 0.5) +
  facet_grid(.~Setting)

if(project != 'polyester_mm_100bp_fG'){
  
  panel3.b <- dete_supplement_theme(panel3.b)
  
}else{
  
  panel3.b <- dete_main_theme(panel3.b)
  
}



### Store raw plot

save(panel3.b, file = paste0(project, '/Figures/raw/figure_3b.Rdata'))

panel3.b <- shade_strips(panel3.b)


png(paste0(project, "/Figures/figure_3b.png"), 
    width = 4.1, height = 2.4, units = 'in', res = 150)
grid.draw(panel3.b) # Make plot
dev.off()



write.csv(deseq.results.young, file = paste0(project, '/Figures/tmp/figure_3b.csv'))

## log2Fold Change of FPs --------------------------------------------------

## Since especially false positives are responsible for the fdr rates in the 
## plot before I want to check how is the expression and log2FC of such false 
## positives?
## For this I filtered the data frame to get only FP with a median bigger than 
## 0. The median ensures to take only instances that are detected in at least 
## the half of the replicates. For the log2FC I divide the mean expression of 
## condition 1 by the mean expression of condition 2. 
## I added 1 to the mean before, to avoid the case of division by zero. 
## The diagonal arrangement of dots are TEs that were only detected in one
## condition. 

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

df.fp$Tool <- factor(df.fp$Tool, levels=c('SalmonTE', 'Telescope', 'TEtranscripts', 'SQuIRE', 'TEtools')) 
df.fp$Setting <- as.factor(df.fp$Setting)

df.fp$Setting <- factor(df.fp$Setting, levels=c('Single-end setup', 'Paired-end setup'))

panel3.c <- df.fp %>% filter(abs(log2Fold)>0.5) %>% 
  ggplot(aes(x=log2Fold, y=log(meanCount))) +
  geom_point(aes(color = log2Fold), alpha=0.1, size=0.5) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=1, bins=50, geom='polygon') +
  scale_fill_gradient(low = 'yellow', high = 'red') +
  scale_alpha(range = c(0.00, 1.0), guide = FALSE) +
  geom_density2d(colour='black', bins=30, size=0.25) +
  scale_x_continuous(breaks = c(-6 , -3, 0, 3, 6), 
                     labels = c(-6,-3,0,3,6)) + 
  labs(x = "log fold change",
       y = "log mean read count") +
       #title = 'False positives (|log fold change| > 0.5)') +
  guides(alpha=FALSE) +
  facet_grid(Tool~Setting) + 
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.minor = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background=element_rect(colour="grey", fill="white"),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = 'bold'),
        title = element_text(size = 12, color = '#414141'))




### Store raw plot


save(panel3.c, file = paste0(project, '/Figures/raw/figure_3c.Rdata'))

panel3.c <- color_strips(panel3.c)

png(paste0(project, "/Figures/figure_3c.png"), 
    width = 8, height = 8, units = 'in', res = 600)
grid.draw(panel3.c) # Make plot
dev.off()

write.csv(df.fp %>% filter(abs(log2Fold)>0.5), 
          file = paste0(project, '/Figures/tmp/figure_3c.csv'))




## Arrange Panel 3 ---------------------------------------------------------


panel3 <- ggdraw()+
  draw_plot(panel3.a, x=0.05, y =.80, width = 0.40, height = 0.18) +
  draw_plot(panel3.b, x=0.55, y =.80, width = 0.40, height = 0.18) +
  draw_plot(panel3.c, x=0, y =0, width = 1, height = 0.78) +
  draw_plot_label(label = c('A', 'B', 'C'), size = 16, 
                  family = "Source Sans Pro",
                  x = c(0,0.5,0), y = c(1,1,0.8))

png(paste0(project, "/Figures/figure_3.png"), width = 1158, height = 1479, units = 'px', family=  "Source Sans Pro")
panel3
dev.off()
  

# Supplemental figures -----------------------------------------------------



## F-Score of Kimura intervals ----

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

metrics.kimura$Setting <- ordered(metrics.kimura$Setting, levels=c('Single-end setup', 'Paired-end setup'))
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
        axis.text = element_text(family = "Source Sans Pro", size = 10),
        axis.title = element_text(family = "Source Sans Pro", size = 12),
        strip.text = element_text(family = "Source Sans Pro", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)        )

pl.supplement.3 <- color_strips(pl.supplement.3)

png(paste0(project, "/Figures/figure_S2.png"), width = 12.5, height = 10, units = 'in', res = 600)
grid.draw(pl.supplement.3) # Make plot
dev.off()

ggsave(filename = paste0(project, "/Figures/figure_S2.pdf"), 
       grid.draw(pl.supplement.3),
       width = 10, height = 7, dpi = 600, units = "in", device=cairo_pdf)

save(pl.supplement.3, file = paste0(project, '/Figures/raw/figure_S2.raw.Rdata'))

write.csv(metrics.kimura, file = paste0(project, '/Figures/tmp/figure_S2.csv'))


## Deviation figure for Set 1 ----------

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

tool_cor.set1 <- total.deviation.set1 %>% 
  filter(mean.recov.log >= log(5), mean.sim.log >= log(5)) %>% 
  group_by(Tool, Setting) %>% 
  summarise(r.squared = summary(lm(mean.sim.log~mean.recov.log))$r.squared) %>% 
  mutate(mean.sim.log = 9.5, 
         mean.recov.log = 2,
         label = glue("*r*<sup>2</sup>(TP) = {round(r.squared, 2)}"))

save(total.deviation.set1, file = paste0(project, '/Data/total.deviation.set1.raw.Rdata'))


sup.fig.3 <- ggplot(total.deviation.set1, aes(mean.sim.log, mean.recov.log)) +
  geom_point(alpha = 0.025, color = 'grey') +
  geom_density2d(linetype=1, size= 0.25, bins=30, aes(color = Tool)) +
  geom_abline(lty=2, color='black', size=1) +
  geom_quadrant_lines(yintercept = log(5), xintercept = log(5)) +
  stat_quadrant_counts(xintercept = log(5), 
                       yintercept = log(5), 
                       size = 6.5, #5.6
                       geom = "label_npc",
                       quadrants = c(1,2,4)) +
  facet_grid(Tool~Setting)+
  my_theme() +
  scale_y_continuous(breaks = c(0,2,4,6,8),
                     labels = c(0,2,4,6,8)) +
  scale_x_continuous(breaks = c(0,2,4,6,8),
                     labels = c(0,2,4,6,8)) +
  scale_color_manual(name = 'Tool', values = my_colors) +
  labs(x = expression(paste("log(", baseMean[Simulated], ")")),
       y = expression(paste("log(", baseMean[Recovered], ")"))) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 22)) + 
  aes(color = Tool) +
  geom_richtext(data = tool_cor.set1,
                aes(label = label, fill = after_scale(alpha(color, .2))),
                text.color = 'black',
                size = 5, #5.6,
                family = 'Source Sans Pro',
                hjust = 1, vjust = 0) 

sup.fig.3  <- color_strips(sup.fig.3 )

save(sup.fig.3 , file = paste0(project, '/Figures/raw/figure_S3.raw.Rdata'))

png(paste0(project, "/Figures/figure_S3.png"), width = 718, height = 1046, units = 'px', family = "Source Sans Pro")
grid.draw(sup.fig.3) # Make plot
dev.off()

write.csv(total.deviation.set1, file = paste0(project, '/Figures/tmp/figure_S3.csv'))
write.csv(tool_cor.set1, file = paste0(project, '/Figures/tmp/figure_S3_Rsquared.csv'))




## Kimura distance shift ----



simdf <- loadRdata(paste0(project, '/Data/Simulation.paired.cnttbl.processed.Rdata'))


simdf <- simdf %>% splitTEID('TE')

if(project != "polyester_nfu_100bp_fG"){
simdf.young <- simdf %>% 
  filter(as.numeric(Kimura) < 5, order %in% c('SINE', 'LTR', 'LINE')) %>% 
  mutate(Kimura = as.numeric(Kimura),
         Kimura = round(Kimura, digits = 1))
}else{
  simdf.young <- simdf %>% 
    filter(as.numeric(Kimura) < 5, order %in% c('SINE', 'LTR', 'LINE', 'DNA')) %>% 
    mutate(Kimura = as.numeric(Kimura),
           Kimura = round(Kimura, digits = 1))
}


pl.supplement.2 <- ggplot(simdf.young, aes(round(Kimura, 1))) +
  geom_bar(width = 0.1, color = 'darkgrey', fill = 'grey')+
  facet_grid(.~order) +
  #scale_y_continuous(limits = c(0,155), expand = c(0,0)) +
  theme_bw() +
  labs(x = 'Kimura distance',
       y = 'TE instance count') +
  theme(panel.grid.major = element_blank(),
        legend.position = 'none',
        panel.spacing.x = unit(1, 'lines'),
        panel.spacing.y = unit(1, 'lines'),
        strip.background=element_rect(colour="grey", fill="white"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 18))

save(pl.supplement.2, file = paste0(project, '/Figures/raw/figure_S1.Rdata'))

png(paste0(project, "/Figures/figure_S1.png"), width = 780, height = 300, units = 'px', family = 'Source Sans Pro')
print(pl.supplement.2) # Make plot
dev.off()

write.csv(simdf.young, 
          file = paste0(project, '/Figures/tmp/figure_S1.csv'))

# Additional Figures ------------------------------------------------------


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

write.csv(total.deviation.set1, file = paste0(project, '/Figures/tmp/figure_S3.csv'))
write.csv(r.squared.set1, file = paste0(project, '/Figures/tmp/figure_S3_Rsquared.csv'))

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


sink("sessionInfo.log")
sessionInfo()
sink() 
