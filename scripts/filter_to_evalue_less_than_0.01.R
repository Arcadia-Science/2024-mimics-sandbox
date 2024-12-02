library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--input"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output CSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

EVALUE <- 0.01

foldseek_results <- read_tsv(args$input, show_col_types = FALSE) 

foldseek_results_filtered <- foldseek_results %>%
  filter(evalue < EVALUE)

write_csv(foldseek_results_filtered, args$output)
