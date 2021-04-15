

#================================= To-Do ======================================
#================================ Functions ==================================

runDESeq <- function(count.table, 
                     col.data, 
                     condition,
                     results=TRUE){
  require(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = count.table,
                                colData = col.data,
                                design = ~condition)
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  print(summary(res))
  
  if(results){
    return(res)
  }else{
    return(dds)
  } 
}



applyDESeq <- function(tool.information, 
                       tool, 
                       sample.assignment, 
                       project.name = NULL){
 
  
  print(paste0(tool.information$output,
               tool,
               ".",
               tool.information$setting,
               tool.information$file.extension))

  count.table <- loadRdata(paste0(tool.information$output,
                                  tool,
                                  ".",
                                  tool.information$setting,
                                  tool.information$file.extension))


  if(count.table %>% filter(duplicated(TE)) %>% nrow() > 0){
    stop("There exist double entries in your count table.")
  }else if(count.table %>% filter(is.na(TE)) %>% nrow() > 0){
    stop("There are NA within your TE column.")
  }
  
  if(project.name == 'dev'){
    print('Dev Modus') 
    print('20.000 items are sampled...')
    count.table <- sample_n(count.table, 20000)
  }

  count.table <- count.table %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    filter(!is.na(TE))

  count.table <- column_to_rownames(remove_rownames(count.table),
                                            var = "TE")


  col.data <- sample.assignment

  col.data <- col.data[order(match(col.data[,1], colnames(count.table))),]

  col.data <- column_to_rownames(remove_rownames(col.data),
                                 var = 'Samples')

  deseq.results <- runDESeq(round(count.table),
                             col.data,
                             "condition")

  deseq.res <- as.data.frame(deseq.results)

  deseq.res$Tool <- tool
  
  deseq.res <- rownames_to_column(deseq.res, var = 'TE')

  save(deseq.res, file = paste0(tool.information$output,
                                tool,
                                ".",
                                tool.information$setting,
                                ".deseq2.res.Rdata"))

  return(deseq.res)
}

merge.list <- function(x, ...) {
  f <- function(x, y) merge(x, y, ...)
  Reduce(f, x)
}

join.list <- function(x, ...){
  
  f <- function(x,y) full_join(x, y, ...)
  Reduce(f, x)

}

createRow <- function(tmp.row, dummy.row){
  
  new.row <- dummy.row
  
  new.row[1:ncol(new.row)] <- NA
  new.row$TE <- tmp.row$TE
  new.row$Tool <- tmp.row$Tool
  new.row$padj <- 1
  return(new.row)
  
}

missingInstancesToDF <- function(df){
  
  missing.instances <- df[[names(df)]] %>% 
    mutate(TE = as.character(TE)) %>% 
    pull(TE)
  
  print(length(missing.instances)) 
  
  df.add <- data.frame(TE = missing.instances,
                       baseMean = NA,
                       log2FoldChange = NA,
                       lfcSE = NA,
                       stat = NA,
                       pvalue = NA,
                       padj = 1,
                       Tool = names(df)
                      )
  
  return(df.add)

}

fill.up.dfs <- function(list.of.dfs, df.to.extend){
  
  # This function merges a list of data frames and stores the output in a tmp 
  # data frame. This tmp df is used to figure out which instance wasn't
  # detected by a certain tool and is then used to create a new row of the 
  # undetected instance for the certain tool, which is added the the original 
  # df. The new line is filled up with NAs except for the TE id, Tool and the
  # pvalue is set to 1 (considered as not detected).
  # I used a double for loop to interate across the cells of the tmp df. 
  #
  #
  
  print('Merge lists.')
  
  tmp <- join.list(list.of.dfs, by = 'TE')
  names(tmp)[2:ncol(tmp)] <- names(list.of.df)
  tmp.NAs <- apply(tmp, 1, function(x){any(is.na(x))})
  
  tmp.filtered <- tmp[tmp.NAs,] %>% 
    gather("Tool", "Present", 2:ncol(tmp)) %>% 
    filter(is.na(Present), Tool != 'Simulation') #Remove simulation because to fill up this one is not necessary. 
   
  tmp.filtered <- split(tmp.filtered, tmp.filtered$Tool)
  
  print('Done. Fill up the exisiting table')
  
  
  df.to.add <- do.call('rbind', mclapply(1:length(tmp.filtered),
                                       function(x) missingInstancesToDF(tmp.filtered[x]),
                                       mc.cores = 1))

  df.to.extend <- rbind(df.to.extend %>% filter(Tool != 'Simulation'), df.to.add)

  return(df.to.extend)

}

determine.ground.truth <- function(simulation.deseq.results, set = NULL){
  
  detected.detes <- simulation.deseq.results %>% 
    dplyr::filter(padj < 0.05) %>% 
    pull(TE)
  
  if(length(detected.detes) > 5000){
    
    return(detected.detes)
    
  }else{
    simulated.detes <- loadRdata(paste0(result.dirs[['data']], set, ".simCounts.full.Rdata")) %>% 
      dplyr::filter(diff == 'TRUE') %>% 
      pull(instance)
    
    return(intersect(simulated.detes, detected.detes))
  }  
}


