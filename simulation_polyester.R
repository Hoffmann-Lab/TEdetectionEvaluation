library(polyester)
library(Biostrings)
library(tidyverse)
library(argparser)

loadRdata <- function (filename) 
{
  
  load(filename)
  get(ls()[ls() != "filename"])
  
}

random_fold_change <- function(number.of.instances, percent.of.dete = 5, seed = 3){
  
  number.of.detes <- floor(percent.of.dete*number.of.instances/100)
  
  detes.up <- floor(number.of.detes/2)
  detes.down <- number.of.detes - detes.up
  
  set.seed(seed)
  log2FC.up <- sample(2:10, detes.up, replace = TRUE) 
  
  set.seed(seed+1)
  
  log2FC.down <- sample(2:10, detes.down, replace = TRUE)
  
  
  fold_changes <- matrix(c(log2FC.down, 
                           rep(1,number.of.instances - detes.down),
                           rep(1,detes.down),
                           log2FC.up, 
                           rep(1,number.of.instances - number.of.detes)), 
                         nrow = number.of.instances)
  
  # shuffle the rows of the matrix
  
  set.seed(seed+2)
  
  rand <- sample(nrow(fold_changes))
  
  fold_changes <- fold_changes[rand, ]
  
  return(fold_changes)
  
}


parser <- argparser::arg_parser("Simulation with polyester")
parser <- argparser::add_argument(parser, "--dete", help="percentage of DETEs", default=5, type = "numeric")
parser <- argparser::add_argument(parser, "--fa", help="fasta file with instances that shall be simulated", type = "character")
parser <- argparser::add_argument(parser, "--replicates", help="number of replicates per group", type = "numeric", default = 5)
parser <- argparser::add_argument(parser, "--setup", help="setup of simulation singel/paired", type = "character", default = "single")
parser <- argparser::add_argument(parser, "--length", help="read length", type = "numeric", default = 100)
parser <- argparser::add_argument(parser, "--output", help='Name of the result directory', type = 'character', default = 'simulated_data_set')
argv <- parse_args(parser)



repeat.fasta <- readDNAStringSet(argv$fa)

fold_changes <- random_fold_change(length(repeat.fasta), percent.of.dete = argv$dete)

if(argv$setup == 'single'){
  setup.paired = FALSE
}else{
  print('paired end data will be simulated')
  setup.paired = TRUE
}

readspertx = round(20 * width(repeat.fasta) / 100)

#simulation without gc_bias because an error occurs
simulate_experiment(argv$fa,
                    num_reps = c(argv$replicates, argv$replicates),
                    fold_changes = fold_changes,
                    reads_per_transcript = readspertx,
                    #meanmodel = TRUE,
                    outdir = paste0(argv$output, '_tmp'),
                    paired = setup.paired,
                    readlen = argv$length,
                    error_model = 'illumina5',
                    seed=234)

# re-read the count table of the simulation before and add a gc bias and simulate
# again

count_matrix <- loadRdata(paste0(argv$output, '_tmp/sim_counts_matrix.rda'))

# Generate a GC bias 
set.seed(5)
gcmodel <- sample(1:7, 2*argv$replicates, replace = TRUE)


count_matrix.biased <- add_gc_bias(count_matrix, gcmodel, repeat.fasta)

# Substitute NAs by 0 and add 1 too all entries, because it doesn't work with NA or 0.
count_matrix.biased[is.na(count_matrix.biased)] <- 0
count_matrix.biased <- count_matrix.biased + 1

simulate_experiment_countmat(argv$fa,
                             readmat = count_matrix.biased,
                             outdir = argv$output,
                             paired = setup.paired,
                             readlen = argv$length,
                             error_model = 'illumina5',
                             seed = 234)

# Evaluation specific


count_matrix.biased <- as.data.frame(count_matrix.biased)

system(sprintf("cp %s_tmp/sim_tx_info.txt %s/", argv$output, argv$output))

fold_change_info <- read.csv(paste0(argv$output, '_tmp/sim_tx_info.txt'), sep = '\t')

fold_change_info <- fold_change_info %>% 
  dplyr::rename(TE = transcriptid) %>% 
  mutate(diff = case_when(DEstatus.1 | DEstatus.2 ~ TRUE,
                          TRUE ~ FALSE)) %>% select(TE, diff)
  
count_matrix.biased <- merge(fold_change_info, count_matrix.biased, 
                             by.x = 'TE', by.y = "row.names") %>%  
  column_to_rownames(var = "TE")

colnames(count_matrix.biased) <- c('diff',
                                   'sample_1',
                                   'sample_2',
                                   'sample_3',
                                   'sample_4',
                                   'sample_5',
                                   'sample_diff_1',
                                   'sample_diff_2',
                                   'sample_diff_3',
                                   'sample_diff_4',
                                   'sample_diff_5')

write.csv(count_matrix.biased, 
          file = paste0(argv$output, '/countTable.csv'),
          quote = F)
