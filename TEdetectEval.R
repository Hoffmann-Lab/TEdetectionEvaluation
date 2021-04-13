library(tidyverse)
library(parallel)

# load settings
source("general.R")

invisible(sapply(list.files('libs','*.R', full.names = TRUE), source))

conflictHandling()

#================================= ToDO =======================================
#
# - add the possibility to avoid the translation process for SQuIRE
# - Implement a if else condition that checks if the combined tool data already
#   exists. Additionally, you can add a condition if you want to calculate the
#   data again. This would avoid the afford to remove the data all time when
#   you want to recalculate the things.
#==============================================================================


#=============================== Questions ====================================
#
# Why there are TE ids named as NA within the SQuIRE count table? Is this a 
# result because of the translation process?
#
# Why I do filter for simulated.expressed and recovered.expressed? This should
# be done by the data processing before, shouldn't be?
#
# Mon Oct  5 17:07:37 2020 ------------------------------
# 
# Currently, I remove each instance with an NA in their TE column, which affects
# especially SQuIRE. Does this have an impact to the results? The reason for
# the NAs is the translation process. When their is no intersection between an
# instance of my annotation and and instances of squires annotation there will
# be an NA generated in the TE column.
#
# I think this is a misunderstanding. I keep the squire annotation when the
# ID is not intersect with an instance of my annotation, therefore I have
# most likely to types of IDs in the squire data frame. One separated by 
# colons (Squire anno) and the other only separated by pipe symbols.
#
#==============================================================================

if(project.name == 'dev'){
  print('ATTENTION THE DEV MODUS IS STARTED!')
}

#============================ Expression detection ============================


for(setting in 1:length(data)){
  
  data.df <- data[[setting]]
  
  if(!file.exists(paste0(result.dirs[['data']],
                         levels(as.factor(data.df$setting)),
                         '.combined', 
                         file.extension))) {
    
      print('Load the count tables of the tools.')
    
      tools <- data.df %>% 
        rownames_to_column(var = 'Tool') %>% 
        pull(Tool)
    
      data.combined <- do.call("rbind",
                              mclapply(tools,
                                        function(x)
                                          countTableHandler(data.df[x,],
                                                            read.threshold),
                                        mc.cores = determineCores(nrow(data.df))))
    
      print('Done.')
    
      print('Rearrange the count table and generate table for expression analysis')
      
      if(project.name == 'dev'){
        print('Dev Modus') 
        print('20.000 items are sampled...')
        data.combined <- sample_n(data.combined, 20000)
      }
      
      data.combined.rotated <- generateExpAnaTab(data.combined)
      data.combined.rotated$Setting <- levels(as.factor(data.df$setting))
    
      print('Done.')
    
      print('Store the final count table.')
    
      save(data.combined.rotated,
          file = paste0(result.dirs[['data']],
                        levels(as.factor(data.df$setting)),
                        '.combined',
                        file.extension))
    
      rm(data.combined.rotated, data.combined)
    
  } else{
    
      print('Table already exists.')
    
  }
  
    print('Done.')
  
}


#========================== DESeq 2 ============================================
# 

deseq.res.all <- data.frame()

for(setting in 1:length(data)){
  
  data.df <- data[[setting]]
  
  if(!file.exists(paste0(result.dirs[['data']],
                         levels(as.factor(data.df$setting)),
                         ".combined.deseq2.res.Rdata"))) {

  print('Run DESeq2...')
    
  tools <- data.df %>% 
      rownames_to_column(var = 'Tool') %>% 
      pull(Tool)

   deseq.results.combined  <- do.call("rbind",
                                     mclapply(tools,
                                              function(x)
                                                applyDESeq(data.df[x, ],
                                                           x,
                                                           sample.assignment,
                                                           project.name),
                                              mc.cores = determineCores(nrow(data.df))))

    
    
    save(deseq.results.combined, file = paste0(result.dirs[["data"]],
                                               levels(as.factor(data.df$setting)),
                                               '.combined.deseq2.res.Rdata'))
    
    deseq.results.combined$setting <- levels(as.factor(data.df$setting))
    deseq.res.all <- rbind(deseq.res.all, deseq.results.combined)
    
    
    
  }else{

    print('Load combined DESeq2 Table.')

    deseq.results.combined <- loadRdata(paste0(result.dirs[["data"]],
                                               levels(as.factor(data.df$setting)),
                                               '.combined.deseq2.res.Rdata'))
    
    deseq.results.combined$setting <- levels(as.factor(data.df$setting))
    deseq.res.all <- rbind(deseq.res.all, deseq.results.combined)

  }
}

print('Done.')

 
print('Prepare data frame for recall vs fdr plot')

# Instances that were detected as differentially expressed out of the simulation
# table were considered as truly differentially expressed for the following
# analysis.

for(setting in 1:length(data)){
  
  data.df <- data[[setting]]
  set <- levels(as.factor(data.df$setting))
                
  ground.truth.set <- determine.ground.truth(deseq.res.all %>%
                                             dplyr::filter(Tool == 'Simulation',
                                                           setting == set), set)

 
  save(ground.truth.set, file = paste0(result.dirs[['data']],set, '.trueDETEs.Rdata'))

}

# I add the Simulation data frame to the list of dfs. But only these instances
# that are contained in the ground truth. I think by adding this set to the
# list I can avoid a second step.

deseq.res.all <- data.frame()

for(setting in 1:length(data)){
  
  data.df <- data[[setting]]
  set <- levels(as.factor(data.df$setting))
  
  ground.truth.set <- loadRdata(file = paste0(result.dirs[['data']],
                                              set, 
                                              '.trueDETEs.Rdata'))
  
  deseq.results.combined <- loadRdata(file = paste0(result.dirs[["data"]],
                                                    set,
                                                    '.combined.deseq2.res.Rdata'))
  list.of.df <- split(deseq.results.combined %>%
                        dplyr::select(c('TE', 'Tool')),
                        deseq.results.combined$Tool)

  list.of.df[['Simulation']] <- list.of.df[['Simulation']] %>%
    dplyr::filter(TE %in% ground.truth.set)

  deseq.results.updated <- fill.up.dfs(list.of.df, deseq.results.combined)

  rm(deseq.results.combined)

  deseq.results.updated <- deseq.results.updated %>%
    mutate(padj = as.numeric(padj),
           padj = case_when(is.na(padj) ~ 1,TRUE ~ padj)) %>%
    filter(Tool != 'Simulation')

  save(deseq.results.updated, file = paste0(result.dirs[['data']],
                                            levels(as.factor(data.df$setting)),
                                            '.combined.deseq2.res.processed.Rdata'))

  
  rm(deseq.results.updated)

}



sink("sessionInfo.log")
sessionInfo()
sink()


