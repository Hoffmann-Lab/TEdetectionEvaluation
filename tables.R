library(tidyverse)
library(janitor)
source('libs/func.R')
source('general.R')

# General Stuff ----

project <- project.name
## Simulation ----
### How much Instances per order were simulated separated by their Kimura group
### assignment. When I load the processed simulation data than instance that
### are not considered as expressed (therefore as not simulated) are already
### removed.


df.order.kim.count <- loadRdata(paste0(project, '/Data/Simulation.paired.cnttbl.processed.Rdata'))


df.order.kim.count <- df.order.kim.count %>% 
  splitTEID('TE') %>% 
  mutate(Kimura = as.numeric(Kimura),
      Kim.grp = case_when(
      Kimura < 5 ~ '[0,5)',
      Kimura >= 5 & Kimura < 10 ~ '[5,10)',
      Kimura >= 10 & Kimura < 15 ~ '[10,15)',
      Kimura >= 15 & Kimura < 20 ~ '[15,20)',
      Kimura >= 20 & Kimura < 25 ~ '[20,25)',
      Kimura >= 25 & Kimura < 30 ~ '[25,30)',
      Kimura >= 30 & Kimura < 35 ~ '[30,35)',
      Kimura >= 35 & Kimura < 40 ~ '[35,40)',
      Kimura >= 40 & Kimura < 45 ~ '[40,45)',
      TRUE ~ '> 45'
    ),
    Kim.grp = factor(
      Kim.grp,
      levels = c(
        '[0,5)',
        '[5,10)',
        '[10,15)',
        '[15,20)',
        '[20,25)',
        '[25,30)',
        '[30,35)',
        '[35,40)',
        '[40,45)',
        '> 45')))

df.order.kim.count <- df.order.kim.count %>% 
  filter(order %in% c('DNA', 'LINE', 'LTR', 'SINE')) %>% 
  count(order, Kim.grp) %>%  spread(Kim.grp, n) %>%  
  mutate(sum = rowSums(.[-1], na.rm = TRUE)) %>%
  adorn_totals("row")

write.csv(df.order.kim.count, file = paste0(project, '/tables/order.kim.count.csv'))


## Mapping Rate ----
# How much reads were mapped per tool and setup? Therefore, I just sum up
# the read counts per column in the count table and calculate the propoetion
# in relation to the simulated ones.

files <- collectFiles(paste0(project, '/Data'), 'mapped.reads.Rdata')

df <- do.call('rbind', lapply(files, function(x) loadRdata(x)))

sim <- df %>% 
  dplyr::filter(Tool == 'Simulation') %>% 
  dplyr::select(Sample, Setting, mappedReads) %>% 
  dplyr::rename(simulatedReads = mappedReads)

df.fig <- df %>% 
  dplyr::filter(Tool != 'Simulation') %>% 
  merge(sim, by=c('Sample', 'Setting')) %>% 
  mutate(per_recovered = (mappedReads/simulatedReads)*100) %>% 
  group_by(Tool, Setting) %>% 
  summarise(sum.alignments = sum(mappedReads),
            sum.simulated = sum(simulatedReads),
            percentage = (sum.alignments/sum.simulated)*100)

write.csv(df.fig, file=paste0(project, '/tables/mapping.rate.csv'))

## Deviation stuff - Figure 2 ----------------------------------------------
### Questions: 
### - How much true positives are underestimated?
### - How much true positives are overestimated?
### - What is the median percent of over and underestimated instances?

### General Condition counts ----
# Determine the mean read counts across the replicates. Determine the Conditions
# of TP, FN and FP. Count the conditions and calculate the percent and store
# them in a .csv file.

df.single <- loadRdata(paste0(project, '/Data/single.combined.cnttbl.processed.Rdata'))
df.paired <- loadRdata(paste0(project, '/Data/paired.combined.cnttbl.processed.Rdata'))

df.all <- rbind(df.single, df.paired)

df.condition <- df.all %>%
  group_by(Tool, TE, Setting, Simulation) %>% 
  summarise(mean.recov = round(mean(recovered.counts),  digits = 0),
            mean.sim = round(mean(simulated.counts), digits = 0),
            mean.recov.log = case_when(mean.recov > 0 ~ log(mean.recov),
                                       TRUE ~ 0),
            mean.sim.log = case_when(mean.sim > 0 ~ log(mean.sim),
                                     TRUE ~ 0),
            .groups = 'drop') %>% 
  mutate(Condition = case_when(mean.recov.log >= log(5) & mean.sim.log >= log(5) ~ "TP",
                               mean.recov.log < log(5) & mean.sim.log >= log(5) ~ "FN",
                               mean.recov.log >= log(5) & mean.sim.log < log(5) ~ "FP",
                               TRUE ~ 'TN')) %>% 
  count(Condition, Tool, Simulation, Setting, name = 'count') %>% 
  spread(Condition, count) %>% 
  mutate(Sum = TP+FP+FN,
         TP.per = (TP/Sum)*100,
         FP.per = (FP/Sum)*100,
         FN.per = (FN/Sum)*100)

df.condition.final <- df.condition %>% 
  mutate(TP = paste(TP, paste0('(', round(TP.per, digits = 1), ')')),
                        FP = paste(FP,  paste0('(', round(FP.per, digits = 1), ')')),
                        FN = paste(FN,  paste0('(', round(FN.per, digits = 1), ')'))) %>% 
  dplyr::select(Tool, Simulation, Setting, TP, FP, FN, Sum)

# Sort the Table
df.condition.final <- df.condition.final[with(df.condition.final, order(Setting, Simulation)),]

write.csv(df.condition.final, paste0(project, '/tables/condition.count.csv'))

### True positive stuff ----

df.TP.behavior <- df.all %>% 
  group_by(Tool, TE, Setting, Simulation) %>% 
  summarise(mean.recov = round(mean(recovered.counts),  digits = 0),
            mean.sim = round(mean(simulated.counts), digits = 0),
            mean.recov.log = case_when(mean.recov > 0 ~ log(mean.recov),
                                       TRUE ~ 0),
            mean.sim.log = case_when(mean.sim > 0 ~ log(mean.sim),
                                     TRUE ~ 0),
            .groups = 'drop') %>% 
  mutate(Condition = case_when(mean.recov.log >= log(5) & mean.sim.log >= log(5) ~ "TP",
                               mean.recov.log < log(5) & mean.sim.log >= log(5) ~ "FN",
                               mean.recov.log >= log(5) & mean.sim.log < log(5) ~ "FP",
                               TRUE ~ 'TN')) %>%
  filter(Condition == 'TP') %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'correctly')) %>% 
  count(estimation, Tool, Simulation, Setting, name = 'count') %>% 
  spread(estimation, count) %>%
  mutate(Sum = overestimated+underestimated+correctly,
         overestimated.per = (overestimated/Sum)*100,
         underestimated.per = (underestimated/Sum)*100,
         correctly.per = (correctly/Sum)*100)

df.TP.behavior.final <- df.TP.behavior %>% 
  mutate(correctly = paste(correctly, paste0('(', round(correctly.per, digits = 1), ')')),
         overestimated = paste(overestimated,  paste0('(', round(overestimated.per, digits = 1), ')')),
         underestimated = paste(underestimated,  paste0('(', round(underestimated.per, digits = 1), ')'))) %>% 
  dplyr::select(Tool, Simulation, Setting, correctly, overestimated, underestimated, Sum)

# Sort the Table
df.TP.behavior.final <- df.TP.behavior.final[with(df.TP.behavior.final, order(Setting, Simulation)),]

write.csv(df.TP.behavior.final, paste0(project, '/tables/TP.behavior.csv'))



df <- loadRdata(paste0(project, '/Data/total.deviation.set1.raw.Rdata'))

df.tp <- df %>% filter(mean.recov.log >= log(5), mean.sim.log >= log(5))

df.tp <- df.tp %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'strike'))


# How much strikes, over- and unterestimated TEs per tool and setting?

df.tp.counts <- df.tp %>% 
  count(name="count", Tool, Setting, estimation) %>% 
  group_by(Tool, Setting) %>% 
  mutate(percent = (count/(sum(count))*100)) 

write.csv(df.tp.counts, file = paste0(project, '/tables/tp.counts.estimation.csv'))

# Median of underestimation

df.tp.counts %>% 
  ungroup() %>% 
  filter(estimation == 'underestimated') %>%
  group_by(Setting) %>% 
  summarise(median = median(percent))

# Median of overestimation

df.tp.counts %>% 
  ungroup() %>% 
  filter(estimation == 'overestimated') %>% 
  group_by(Setting) %>% 
  summarise(median = median(percent))

# Percent of overestimated TP in specific ranges

# range between 1.6 and 3
df %>% 
  filter(mean.recov.log >= log(5), between(mean.sim.log, 1.6, 3)) %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'strike')) %>% 
  #filter(estimation == 'overestimated') %>%
  group_by(Tool, Setting) %>% 
  #summarise(deviation = mean(diff)-100) %>% 
  count(name="count", Tool, Setting, estimation) %>% 
  #group_by(Tool, Setting) %>% 
  mutate(percent = (count/(sum(count))*100)) %>% View()

# range between 1.6 and 5
df %>% 
  filter(mean.recov.log >= log(5), between(mean.sim.log, 1.6, 5)) %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'strike')) %>% 
  group_by(Tool, Setting) %>% 
  count(name="count", Tool, Setting, estimation) %>% 
  mutate(percent = (count/(sum(count))*100)) %>% View()

# range between 1.6 and 7
df %>% 
  filter(mean.recov.log >= log(5), between(mean.sim.log, 1.6, 7)) %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'strike')) %>% 
  group_by(Tool, Setting) %>% 
  count(name="count", Tool, Setting, estimation) %>% 
  mutate(percent = (count/(sum(count))*100)) %>% View()

# range between 2.5 and 7
df %>% 
  filter(mean.recov.log >= log(5), between(mean.sim.log, 2.6, 7)) %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'strike')) %>% 
  group_by(Tool, Setting) %>% 
  count(name="count", Tool, Setting, estimation) %>% 
  mutate(percent = (count/(sum(count))*100)) %>% View()

# full range
df %>% 
  filter(mean.recov.log >= log(5), mean.sim.log >= log(5)) %>% 
  mutate(diff = (mean.recov/mean.sim)*100,
         estimation = case_when(diff > 100 ~ 'overestimated',
                                diff < 100 ~ 'underestimated',
                                TRUE ~ 'strike')) %>% 
  group_by(Tool, Setting) %>% 
  count(name="count", Tool, Setting, estimation) %>% 
  mutate(percent = (count/(sum(count))*100)) %>% View()


# How is the mean diviation?

df.tp %>% 
  filter(estimation == 'overestimated') %>%
  mutate(per.deviation = 100 - diff) %>% 
  group_by(Tool, Setting) %>% 
  summarise(median = median(per.deviation)) 


# Portions of Condition

df.count.condition <- df %>% mutate(condition = case_when(
  
  mean.recov.log >= log(5) & mean.sim.log >= log(5) ~ 'TP',
  mean.recov.log >= log(5) & mean.sim.log < log(5) ~ 'FP',
  mean.recov.log < log(5) & mean.sim.log > log(5) ~ 'FN',
  TRUE ~ 'TN'
  
)) %>% count(condition, Tool, Setting)


df.count.condition %>% 
  filter(condition == 'FP', Setting == 'single') %>% 
  summarise(median(n))

df.count.condition %>% 
  filter(condition == 'FN', Setting == 'single') %>% 
  summarise(median(n))


sink("sessionInfo.log")
sessionInfo()
sink()
