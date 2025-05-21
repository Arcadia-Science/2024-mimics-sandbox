library(tidyverse)
library(janitor)
setwd("~/github/2024-parasite-human-mimics/")

files <- Sys.glob("outputs/structural_comparison/foldseek_with_uniprot_metadata/*with_uniprot_metadata.tsv")
# col_types <- read_tsv(files[1], guess_max = 10000) %>%
#   map_chr(~class(.)) %>%
#   map(cols)
# 
# foldseek_results <- files %>%
#   map_dfr(read_tsv, col_types = col_types)

col_types_inferred <- read_tsv(files[1], guess_max = 10000) %>%
  map(~cols(.default = col_guess()))

# Use inferred column types to read all files
foldseek_results <- files %>%
  map_dfr(~read_tsv(.x, col_types = col_types_inferred[[1]]))

foldseek_results_top_matches <- foldseek_results %>%
  group_by(query) %>%
  slice_max(qtmscore) %>%
  filter(qtmscore > 0.7 | alntmscore > 0.7 | ttmscore > 0.7) 

interleukin
