#================================= ToDO =======================================
#
# - read a file where the settings are stored
# - sort the columns of the count tables in the same order --> is it necessary??
# - add the possibility to avoid the translation process for SQuIRE
# - Implement a if else condition that checks if the combined tool data already
#   exists. Additionally, you can add a condition if you want to calculate the
#   data again. This would avoid the afford to remove the data all time when
#   you want to recalculate the things.
# - implement a method in bash or python to create the translation process
#   and store the path in a certain file
# - consider the file names, they are currently hard coded but that doesn't work
#   for other data sets
# - write a method for the task to add the mass of information to the rotated
#   count table. And ask if the file already exists, if not than don't do that
#   step.
#
#
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

#==============================================================================


#=============================== Functions ====================================





tetranscriptsHandler <- function(tool.information){
  
  if(!file.exists(paste0(as.character(tool.information[,3]) ,
                         "TEtranscripts", as.character(tool.information[,4])))){
    
    tetranscripts.cnttbl <- readTEtranscripts(paste0(as.character(tool.information[,2]),
                                                     "/mm_simulated.cntTable")) %>% 
      filterSumCounts(threshold = read.threshold)
    
    names(tetranscripts.cnttbl)[names(tetranscripts.cnttbl) == "TE.id"] <- "TE"
    save(tetranscripts.cnttbl, file = paste0(as.character(tool.information[,3]), 
                                             "TEtranscripts",
                                             as.character(tool.information[,4])))
  }else{
    
    tetranscripts.cnttbl <- loadRdata(paste0(as.character(tool.information[,3]), 
                                             "TEtranscripts",
                                             as.character(tool.information[,4])))
  }
  
   
  tetranscripts.cnttbl$Tool <- "TEtranscripts"
  return(tetranscripts.cnttbl)
}

squireHandler <- function(tool.information){
  
  squire.cnttbl <- readSquire(tool.information$path, tool.information$file)
  
  if (squire.translation) {
    translation.map <- loadRdata(tool.information$additional.files)
    translation.map$alternative.annotation <- as.character(translation.map$alternative.annotation) 
    squire.cnttbl <- translate.squire.ids(squire.cnttbl,
                                          translation.map)
    
    rm(translation.map)
    
  }
  
  squire.cnttbl <-
    squire.cnttbl %>% filter(!(is.na(TE) | duplicated(TE)))
  
  return(squire.cnttbl)
  
}

loadFastqOrder <- function(file){
  
  sample.lst <- tryCatch(
    {
      scan(file, what="", sep = "\n", quiet = TRUE)
    },
    error=function(cond){
      message("File orderOfFastqs.txt does not seem to exist.")
      message("Original error message:")
      return(NA)
    }
  )
  return(sample.lst)
}

tetoolsHandler <- function(tool.information){
  
  TETools.sample.lst <- loadFastqOrder(paste0(tool.information$path, 'orderOfFastqs.txt'))
  
  tetools.cnttbl <- readTETools(paste0(tool.information$path,
                                       tool.information$file),
                                TETools.sample.lst)
  
  return(tetools.cnttbl)
  
}

countTableHandler <- function(tool.information, sumCount.threshold){
  
  tool = row.names(tool.information)[1]
  
  print(paste('Load the count table of', tool))
  
  if(!file.exists(paste0(
    tool.information$output,
    tool,
    tool.information$file.extension
  ))) {
    switch (
      tool,
      "Simulation" = count.tbl <-
        {
          simulationHandler(tool.information)
        },
      "salmonTE" = count.tbl <-
        {
          salmonTEHandler(tool.information)
        },
      "TEtranscripts" = count.tbl <-
        {
          tetranscriptsHandler(tool.information)
        },
      "SQuIRE" = count.tbl <- {
        squireHandler(tool.information)
      },
      "TEtools" = count.tbl <-
        {
          tetoolsHandler(tool.information)
        },
      stop("No function for that tool exists.")
    )
    
    count.tbl <-
      count.tbl %>% filterSumCounts(threshold = sumCount.threshold)
    
    save(count.tbl, file = paste0(
      tool.information$output,
      tool,
      tool.information$file.extension
    ))   
    
    count.tbl <- rotateCountTable(count.tbl, h = 'Counts')
    count.tbl$Tool = row.names(tool.information)[1]
    
    
  }else{
    
    count.tbl <- loadRdata(paste0(
      tool.information$output,
      tool,
      tool.information$file.extension
    ))
    
    count.tbl <- rotateCountTable(count.tbl, h = 'Counts')
    count.tbl$Tool = row.names(tool.information)[1]
  }
  
  
  return(count.tbl)
  
}

mergeWithSimulation <- function(tool, countTab, simTab){
  
  #df.a <- as.data.frame(df.a[[tool]])
  
  df <- merge(countTab, simTab, by=c("TE", "sample"), all = TRUE) %>%
    dplyr::rename(recovered.counts = Counts.x,
           simulated.counts = Counts.y,
           Tool = Tool.x) %>%
    select(-c('Tool.y')) %>%
    mutate(Tool = tool)
  
  return(df)
}

