#============================ general stuff ====================================

createResultDirs <- function(path){
 # creates result directory for data and figures that were produced during 
 # the evaluation process
  
  evalData <- paste0(path, "/Data/")
  dir.create(evalData, recursive = TRUE)
  
  evalFigures <-paste0(path, "/Figures/")
  dir.create(evalFigures, recursive = TRUE)
  dir.create(paste0(evalFigures, 'tmp/'), recursive = TRUE)
  dir.create(paste0(evalFigures, 'raw/'), recursive = TRUE)
  dir <- c(evalData, evalFigures)
  names(dir) <- c("data", "figures")
  
  return(dir)
} 

conflictHandling <- function(){
  
  require(conflicted)
  
  conflict_prefer('select', 'dplyr')
  conflict_prefer('count', 'dplyr')
  conflict_prefer('filter', 'dplyr')
  conflict_prefer('rename', 'dplyr')
  
}

collectFiles <- function(path = NULL, file.type = NULL){
  
  ######### 
  # This method collects certain types of files and stores them in
  # a list that is returned.
  #######
  
  if(is.null(file.type) | is.null(path)){
    message("Path or Filetype is missing")
    return(NA)
  }else{
    pattern = paste0("*", file.type, "$")
    list.files(path = path, 
               pattern = pattern,
               full.names = T)
  }
  
}

determineCores <- function(number.of.tasks){
  # puffer to keep some cores free when the script is running
  puffer.cores = 4
  
  if(number.of.tasks < (detectCores()-puffer.cores)){
    
    number.of.cores <- number.of.tasks
    
  }else{
    
    number.of.cores <- detectCores() - puffer.cores
  }
  
  return(number.of.cores)
}

loadRdata <- function (filename) 
{
  
  load(filename)
  get(ls()[ls() != "filename"])

}
#============================= table stuff =====================================

rotateCountTable <- function(df, k = "sample", h = "simulated.counts"){
  
  require(lazyeval) 
  require(tidyr)
  
  gather_(df, names(df)[2:ncol(df)], key = k, value = h)
  
  
}

filterSumCounts <- function(df, threshold=5){
  # this function removes all copies where the rowsum (counts across samples) is
  # not above the threshold
  
  df[rowSums(df[2:ncol(df)]) > threshold,]
  
}

splitTEID <- function(df, column){
  
  require(tidyverse)
  if(column %in% colnames(df)){
    df %>% separate(column,c("chromosome",
                             "strand",
                             "start",
                             "end",
                             "order",
                             "super_family",
                             "family",
                             "id",
                             "Kimura"), 
                    sep = "[|]",
                    remove = F) 
  }else{
    print("Column doesn't exist!")
  }
 
}

tabRS <- function(df = NULL){
  require(kableExtra)
  kable(df) %>% 
    kable_styling(bootstrap_options = c("hover", "condensed", "responsive"))
}

extractSampleNames <- function(simulation.df){
  
  colnames(simulation.df[,grepl('^sample_', colnames(simulation.df))])
}

countMappedReads <- function(df, Tool = NA_character_, Setting = NA_character_){
  
  mappedRead.df = data.frame(Sample = names(df[2:ncol(df)]),
                             Tool = Tool,
                             Setting = Setting,
                             mappedReads = colSums(df[2:ncol(df)])
                             ) %>% remove_rownames()
  
  return(mappedRead.df)
}






























