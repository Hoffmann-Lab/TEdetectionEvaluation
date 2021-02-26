# This functions load the count table of the simulated instance and extract
# the pure counts and removes information like if an instance is differentially
# expressed for example.

loadSimCounts <- function(simFile){
  
  ######
  # This function loads the information table from the simulation. Since the
  # names of the instances are prolonged by bedtools the id's do not fit to
  # the origin. This function removes the extended tail, which is separated 
  # by colons. Since the align file can contain instances multiple times
  # it is possible that instances were simulated multiple times. To avoid
  # trouble I remove these double entries.
  #
  # Mon Oct  5 10:56:43 2020 ------------------------------
  #
  # Removed the step to filter out instances with less than 5 read counts across
  # all samples. It is handled now late in my benchmarking process, so that I
  # can be sure that each table is handled equally and can change it at
  # only one position in my code.
  #
  # Thu Oct 29 15:18:32 2020 ------------------------------
  #
  # Removed the part where 1 read were substracted from each read count. I fixed
  # that issue in my readiator. Now, the counts of the read tables are exactly
  # the number of reads that are contained in the fastq files.
  #
  #
  #####
  
  simulated_read_counts <- read.table(simFile, header = TRUE)
  x <-  nrow(simulated_read_counts)
  
  print(paste(as.character(x), "instances were initially simulated"))
  
  simulated_read_counts <- simulated_read_counts %>% 
    rownames_to_column(var = 'instance') %>% 
    tidyr::separate(instance, c('instance'), sep = ":")
  
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