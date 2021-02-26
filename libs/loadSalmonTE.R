# This function loads the count table from SalmonTE

readSalmonTE <- function(path, file = NULL ){
  
  if(is.null(file)){
    salmonTE.cnttbl <- tryCatch(
      {
        read.csv(paste0(path, "/EXPR.csv"),
                 check.names = FALSE)
      },
      error=function(cond){
        message("The standard EXPR.csv file does not seem to exist.")
        message("Original error message:")
        return(NA)
      })
  }else{
    salmonTE.cnttbl <- tryCatch(
      {
        read.csv(paste0(path, file),
                 check.names = FALSE)
      },
      error=function(cond){
        message("Given count table of salmonTE does not seem to exist.")
        message("Original error message:")
        return(NA)
      })
  }
  return(salmonTE.cnttbl)
}

salmonTEHandler <- function(tool.information){
  
  
  salmonTE.cnttbl <- readSalmonTE(tool.information$path, tool.information$file)
  
  return(salmonTE.cnttbl)
  
}