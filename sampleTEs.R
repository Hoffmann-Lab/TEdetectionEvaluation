library(polyester)
library(Biostrings)
library(data.table)
library(tidyverse)
library(argparser)

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

getChromosomes <- function(reference.fasta){
  
  rnd <- floor(runif(1, min=0, max=123943))
  system(sprintf("grep '^>' %s | sed 's/>//' > chrom%s.tmp", reference.fasta, rnd))
  
  chrom <- read.csv(paste0("chrom", rnd, ".tmp"), header = F)
  
  system(sprintf("rm chrom%s.tmp", rnd))
  
  return(as.vector(chrom$V1))

}

pullRandomTes <- function(bed.file.name, 
                          reference.fasta,
                          organism,
                          repeats.of.interest = c('LINE', 'SINE', 'LTR', 'DNA'),
                          min.length = 100,
                          number.of.instances = 1000,
                          res.dir = './',
                          seed = 42){
  
  
  
  bed.file <- fread(bed.file.name, sep = '\t', header = F)
  total.TEs  <-  nrow(bed.file)
  chrom <- getChromosomes(reference.fasta)
  
  # remove all TEs which are not annotated on a chromosome that are in the
  # reference fastq and not contained in repeats of interest
  bed.file <- bed.file %>% 
    filter(V1 %in% chrom) %>% 
    filter(grepl(paste(repeats.of.interest, collapse = "|"), V4))
  
  # Calculate and filter by length and Kimura distance
  bed.file <- bed.file %>%
    splitTEID('V4') %>%
    mutate(length = as.numeric(end) - as.numeric(start)) %>%
    filter(length > min.length, Kimura != "NA") %>%
    select(V1, V2, V3, V4, V5, V6)

  filtered.TEs <- nrow(bed.file)
  
  set.seed(seed)
  sampled.TEs <- sample_n(bed.file, number.of.instances)
  
  write.table(sampled.TEs, file = paste0(res.dir, '/', organism,'.sampled.TEs.bed'), 
              sep = '\t', col.names = F, row.names = F, quote = F)
  
  
  log.df <- data.frame(setting = c('bed-file', 
                                   'fasta-file', 
                                   'species',
                                   'min.length', 
                                   'num.of.drawn.instances', 
                                   'res.dir', 
                                   'seed', 
                                   'total.TEs', 
                                   'total.filtered.TEs'),
                       values = c(bed.file.name, 
                                  reference.fasta, 
                                  organism, 
                                  min.length, 
                                  number.of.instances, 
                                  res.dir, 
                                  seed, 
                                  total.TEs, 
                                  filtered.TEs))
  
  write.table(log.df, 
              file = paste0(res.dir, '/.TEsampling.log'), 
              sep = '\t', 
              col.names = F, 
              row.names = F, 
              quote = F)
}





parser <- argparser::arg_parser("Sample a random set of TEs")
parser <- argparser::add_argument(parser, "--min", help="minimal length of TE", default=100, type = "numeric")
parser <- argparser::add_argument(parser, "--species", help="species where the TEs are from", type = "character", default = 'None')
parser <- argparser::add_argument(parser, "--output", help='Name of the result directory', type = 'character', default = './')
parser <- argparser::add_argument(parser, "--count", help = 'Count of TEs that shall be extracted', type = 'numeric', default = 1000)
parser <- argparser::add_argument(parser, "--seed", help="seed; default = 42", type = "numeric", default = 42)
parser <- argparser::add_argument(parser, "--bed", help="bed file that contains the annotated TEs", type = "character")
parser <- argparser::add_argument(parser, "--fa", help="reference genome where the TEs are frome", type = "character")
argv <- parse_args(parser)



if(argv$output != './'){
  
  dir.create(argv$output)
  
}

bed.filtered <- pullRandomTes(bed.file = argv$bed,
                              reference.fasta = argv$fa,
                              organism = argv$species,
                              min.length = argv$min,
                              number.of.instances = argv$count,
                              res.dir = argv$output,
                              seed = argv$seed)

