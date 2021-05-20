
#=============================== Settings =====================================

project.name = 'polyester_mm_100bp_fG_3reps'
read.threshold = 10
file.extension = ".cnttbl.processed.Rdata"
squire.translation = TRUE
organism = 'mm'

# Sample Assignment for DESeq2
# 
sample.assignment <-
  data.frame(
    Samples = c(
      'sample_1',
      'sample_2',
      'sample_3',
      'sample_diff_1',
      'sample_diff_2',
      'sample_diff_3'
    ),
    condition = c(
      'control',
      'control',
      'control',
      'diff',
      'diff',
      'diff'
    )
  )

# sample.assignment <-
#   data.frame(
#     Samples = c(
#       'sample_1',
#       'sample_2',
#       'sample_3',
#       'sample_4',
#       'sample_5',
#       'sample_diff_1',
#       'sample_diff_2',
#       'sample_diff_3',
#       'sample_diff_4',
#       'sample_diff_5'
#     ),
#     condition = c(
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff'
#     )
#   )
# 
# sample.assignment <-
#   data.frame(
#     Samples = c(
#       'sample_1',
#       'sample_2',
#       'sample_3',
#       'sample_4',
#       'sample_5',
#       'sample_6',
#       'sample_7',
#       'sample_8',
#       'sample_9',
#       'sample_10',
#       'sample_diff_1',
#       'sample_diff_2',
#       'sample_diff_3',
#       'sample_diff_4',
#       'sample_diff_5',
#       'sample_diff_6',
#       'sample_diff_7',
#       'sample_diff_8',
#       'sample_diff_9',
#       'sample_diff_10'
#     ),
#     condition = c(
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'control',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff',
#       'diff'
#     )
#     )
#=============================== generate data df ==============================

source('libs/func.R')

data <- read.csv('dataInfo.csv', sep = "\t", stringsAsFactors = F)

data <- split(data, data$setting)

result.dirs <- createResultDirs(project.name)

for(i in range(1,length(data))){
 
  data[[i]] <- data[[i]] %>% remove_rownames() %>% column_to_rownames(var = 'Tool')
  data[[i]]$file.extension = file.extension
  data[[i]]$output = result.dirs[['data']]
  
}


