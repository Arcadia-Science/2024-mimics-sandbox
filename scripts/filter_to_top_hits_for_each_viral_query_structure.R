library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--input"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output CSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

foldseek_results <- read_tsv(args$input, show_col_types = FALSE) 

# Select the top hit for each query protein. Start with selecting the lowest
# evalue. If there is a tie, use the maximum TM Score. If there is still a tie,
# arbitrarily select one of the hits.

foldseek_results_filtered <- foldseek_results %>%
  rowwise() %>%
  mutate(min_tmscore = min(alntmscore, qtmscore, ttmscore),
         max_tmscore = max(alntmscore, qtmscore, ttmscore),
         avg_tmscore = sum(alntmscore, qtmscore, ttmscore) / 3 ) %>%
  ungroup() %>%
  group_by(query) %>%
  slice_min(evalue) %>%
  slice_max(max_tmscore) %>%
  slice_head(n = 1) %>%
  ungroup()

write_csv(foldseek_results_filtered, args$output)
