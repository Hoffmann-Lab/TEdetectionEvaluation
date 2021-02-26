# This function loads the count tables from TEtools

readTETools <- function(count.tbl, sample.lst = NULL){
  
  if(is.null(sample.lst)){
    
    print("Sample list does not exist, therefore samples  are named with
          increasing numbers.")
    
    tetools.count.tbl <- tryCatch(
      {
        read.csv(count.tbl, sep=" ", check.names = FALSE)
      },
      error=function(cond){
        message('The count table of TEtools does not seem to exist.')
        message("Original error message:")
        return(NA)
      }
    )
      
    names(tetools.count.tbl) <- c("TE", 
                                    sprintf("sample_%02d", 
                                            seq(1, ncol(tetools.count.tbl)-2)),
                                  "Sum")
    
  }else{
    
    if(length(read.csv(count.tbl, sep = " ", nrows = 1))-2 != length(sample.lst)){
      
      stop("Your sample list has not the correct length.")
    }
    
    tetools.count.tbl <- tryCatch(
      {
        read.csv(count.tbl, sep=" ", check.names = FALSE)
      },
      error=function(cond){
        message('The count table of TEtools does not seem to exist.')
        message("Original error message:")
        return(NA)
      }
      
    )
    
    names(tetools.count.tbl) <- c("TE", sample.lst, "Sum")
    
    }
    
  tetools.count.tbl$Sum <- NULL
  
  return(tetools.count.tbl)  
} 

loadFastqOrder <- function(file){
  
  sample.lst <- tryCatch(
    {
      scan(file, what="", sep = "\n", quiet = TRUE)
    },
    error=function(cond){
      message("File orderOfFastqs.txt does not seem to exist.")
      message("Original error message:")
      return(NULL)
    }
  )
  return(sample.lst)
}

tetoolsHandler <- function(tool.information){
  
  TETools.sample.lst <- loadFastqOrder(tool.information$additional.files)
  
  tetools.cnttbl <- readTETools(paste0(tool.information$path,
                                       tool.information$file),
                                TETools.sample.lst)
  
  return(tetools.cnttbl)
  
}