mergeWithSimulation <- function(tool, countTab, simTab){
  
  df <- merge(countTab, simTab, by=c("TE", "sample"), all = TRUE) %>%
    dplyr::rename(recovered.counts = Counts.x,
                  simulated.counts = Counts.y,
                  Tool = Tool.x) %>%
    dplyr::select(-c('Tool.y')) %>%
    mutate(Tool = tool)
  
  return(df)
}




# Here is the table generated that can be used to generate plots with respect
# to the analysis how good the tools detect the expression of elements.

generateExpAnaTab <- function(count.table){
  
  require(parallel)
  
  count.tables.lst <- split(count.table, count.table$Tool)
  simulation.df <- count.tables.lst[['Simulation']]
  simulation.df$Setting <- NULL
  count.tables.lst[['Simulation']] <- NULL
  
  #print(names(count.tables.lst)) 
  
  cl <- makeCluster(determineCores(nrow(data.df)))
  
  data.combined.rotated <-
    do.call("rbind",
            mclapply(1:length(count.tables.lst),
                     function(x)
                       mergeWithSimulation(tool = names(count.tables.lst[x]),
                                           countTab =  count.tables.lst[[x]],
                                           simTab = simulation.df)))
  
  stopCluster(cl)
  
  
  data.combined.rotated <- data.combined.rotated %>%
    tidyr::separate(
      TE,
      c(
        "chromosome",
        "strand",
        "start",
        "end",
        "order",
        "super_family",
        "family",
        "id",
        "Kimura"
      ),
      sep = "[|]",
      remove = FALSE
    ) %>%
    mutate(
      Simulated = case_when(!is.na(simulated.counts) ~ TRUE,
                            TRUE ~ FALSE),
      sim.expressed = case_when(
        !is.na(simulated.counts) & simulated.counts >= 5 ~ TRUE,
        TRUE ~ FALSE
      ),
      recovered.expressed = case_when(recovered.counts >= 5 ~ TRUE,
                                      TRUE ~ FALSE),
      Condition = case_when(
        sim.expressed & recovered.expressed ~ 'TP',!(sim.expressed) &
          recovered.expressed ~ 'FP',
        sim.expressed &
          !recovered.expressed ~ 'FN',
        TRUE ~ NA_character_
      ),
      Simulation = case_when(grepl('^sample_diff', sample) ~ 'with DETEs',
                             TRUE ~ 'ETEs'),
      Kimura = as.numeric(Kimura),
      Kim.grp = case_when(
        Kimura < 5 ~ '[0,5)',
        Kimura >= 5 & Kimura < 10 ~ '[5,10)',
        Kimura >= 10 & Kimura < 15 ~ '[10,15)',
        Kimura >= 15 & Kimura < 20 ~ '[15,20)',
        Kimura >= 20 & Kimura < 25 ~ '[20,25)',
        Kimura >= 25 & Kimura < 30 ~ '[25,30)',
        Kimura >= 30 & Kimura < 35 ~ '[30,35)',
        Kimura >= 35 & Kimura < 40 ~ '[35,40)',
        Kimura >= 40 & Kimura < 45 ~ '[40,45)',
        TRUE ~ '> 45'
      ),
      Kim.grp = factor(
        Kim.grp,
        levels = c(
          '[0,5)',
          '[5,10)',
          '[10,15)',
          '[15,20)',
          '[20,25)',
          '[25,30)',
          '[30,35)',
          '[35,40)',
          '[40,45)',
          '> 45'
        )
      )#,
      #FC = recovered.counts / simulated.counts
    ) %>%
    filter(sim.expressed | recovered.expressed)
  
  return(data.combined.rotated)
}
