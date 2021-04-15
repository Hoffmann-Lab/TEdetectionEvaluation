
loadTable <- function(sample){
  
  require(stringi)
  
  tab <- read.csv(paste0(sample, '/telescope-telescope_report.tsv'), 
                  sep = '\t', 
                  skip = 1,
                  header = TRUE,
                  stringsAsFactors = TRUE)
  
  tab <- tab %>% dplyr::select(transcript, final_count) %>%
    dplyr::rename(TE = transcript, count = final_count)
  
  sample <-  tail(strsplit(sample, "/")[[1]], n = 1)
  
  if(stri_sub(sample, -1) == '_'){
    
    tab$sample = stri_sub(sample, 1, -2)
  
  }else if(identical(base::strsplit(sample, "[.]")[[1]][2], 'fastq')){
    
    tab$sample = strsplit(sample, "[.]")[[1]][1]
  
  }else{
    
    tab$sample = sample
  }
  
  return(tab)
  
}

readTelescope <- function(result.path){
 
  samples <- list.files(result.path, include.dirs = TRUE, full.names = TRUE)
  count.table <- do.call('rbind', lapply(samples, function(x) loadTable(x)))
  
  count.table <- count.table %>% tidyr::spread(sample, count) %>% 
    mutate(across(everything(), ~replace_na(.x, 0)))
  
  count.table <- count.table[!count.table$TE == '__no_feature',]
  
  return(count.table)
  
}


telescopeHandler <- function(tool.information){
  
  telescope.cnttbl <- readTelescope(tool.information$path)
  
  return(telescope.cnttbl)
  
}
