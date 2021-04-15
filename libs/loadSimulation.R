# This functions load the count table of the simulated instance and extract
# the pure counts and removes information like if an instance is differentially
# expressed for example.

loadSimCounts <- function(simFile){
  
  simulated_read_counts <- read.table(simFile, header = TRUE, sep = ',')
  x <-  nrow(simulated_read_counts)
  
  print(paste(as.character(x), "instances were initially simulated"))
  
  simulated_read_counts <- simulated_read_counts %>% 
    dplyr::rename(instance = X)
    #rownames_to_column(var = 'instance') %>% 
    #tidyr::separate(instance, c('instance'), sep = ":")
  
  simulated_read_counts <- removeDoubles(simulated_read_counts)
  
  samples <- extractSampleNames(simulated_read_counts)
  
  return(simulated_read_counts)
  
}

simulationHandler <- function(tool.information){
  
  
  simCounts <- loadSimCounts(paste0(tool.information$path, tool.information$file))
  
  save(simCounts,
       file = paste0(tool.information$output,
                     tool.information$setting,
                     ".simCounts.full.Rdata"))
  
  simCountsRaw <- extractSimCounts(simCounts)
  
  names(simCountsRaw)[names(simCountsRaw) == "instance"] <- "TE" 
  return(simCountsRaw)
  
}

removeDoubles <- function(simulation.df){
  
  doubles <- simulation.df %>% 
    filter(duplicated(instance)) %>% 
    pull(instance)
  
  print(paste(as.character(length(doubles)*2), "double entries were removed"))
  simulation.df <- simulation.df %>% 
    filter(!(instance %in% doubles))
  
  return(simulation.df)
}

extractSimCounts <- function(sim.df, id = "instance"){
  
  samples <- extractSampleNames(sim.df)
  
  selected.columns <- append(c("instance"), samples) 
  sim.df <- sim.df[selected.columns]
  
  return(sim.df)
}