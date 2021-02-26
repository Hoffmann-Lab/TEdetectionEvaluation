# This function loads the certain count tables of the tools that are
# stored in the that frame that this function get
# After Loading the respective count table each table is handled in the same
# way, which means to remove elements that have less than 5 reads in sum across
# all replicates.


countTableHandler <- function(tool.information, sumCount.threshold){
  
  tool = row.names(tool.information)[1]
  
  print(paste('Load the count table of', tool))
  
  if(!file.exists(paste0(
    tool.information$output,
    tool,
    ".",
    tool.information$setting,
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
      "Telescope" = count.tbl <- 
        {
          telescopeHandler(tool.information)
        },
      stop("No function for that tool exists.")
    )
    
    mappedReads <- countMappedReads(count.tbl, tool, tool.information$setting)
   
    save(mappedReads, file = paste0(
      tool.information$output,
      tool,
      ".",
      tool.information$setting,
      "mapped.reads.Rdata"
    )) 
    
    count.tbl <-
      count.tbl %>% filterSumCounts(threshold = sumCount.threshold)
    
    save(count.tbl, file = paste0(
      tool.information$output,
      tool,
      ".",
      tool.information$setting,
      tool.information$file.extension
    ))   
    
    count.tbl <- rotateCountTable(count.tbl, h = 'Counts')
    count.tbl$Tool = tool
    count.tbl$Setting = tool.information$setting
    
    
  }else{
    
    count.tbl <- loadRdata(paste0(
      tool.information$output,
      tool,
      ".",
      tool.information$setting,
      tool.information$file.extension
    ))
    
    count.tbl <- rotateCountTable(count.tbl, h = 'Counts')
    count.tbl$Tool = tool
    count.tbl$Setting = tool.information$setting
  }
  
  
  return(count.tbl)
  
}


