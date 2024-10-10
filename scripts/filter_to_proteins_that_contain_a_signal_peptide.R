library(tidyverse)
library(janitor)
library(optparse)

option_list <- list(
  make_option(c("--input"), type="character",
              help="Path to uniprot metadata TSV file."),
  make_option(c("--output_tsv"), type="character",
              help="Path to output metadata TSV file."),
  make_option(c("--output_txt"), type="character",
              help="Path to output protein identifier TXT file.")
)

args <- parse_args(OptionParser(option_list=option_list))


uniprot_metadata_with_signal_peptide <- read_tsv(args$input) %>%
  clean_names() %>%
  filter(!is.na(signal_peptide)) %>%
  filter(str_detect(string = signal_peptide, pattern = "SIGNAL"))

write_tsv(uniprot_metadata_with_signal_peptide, file = args$output_tsv)
write_tsv(uniprot_metadata_with_signal_peptide %>% select(protid),
          file = args$output_txt, col_names = FALSE)
