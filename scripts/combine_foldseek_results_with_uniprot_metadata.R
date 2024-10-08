library(tidyverse)
library(janitor)
library(optparse)

option_list <- list(
  make_option(c("--input_foldseek_results"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--input_human_metadata"), type="character",
              help="Path to human uniprot metadata TSV file."),
  make_option(c("--input_query_metadata"), type="character",
              help="Path to query uniprot metadata TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

foldseek_results <- read_tsv(args$input_foldseek_results) %>%
  mutate(query = str_remove(string = query, pattern = "\\.pdb"),
         target = str_remove(string = target, pattern = "-F[0-9]-model_v4\\.pdb\\.gz"),
         target = str_remove(string = target, pattern = "-F[0-9][0-9]-model_v4.pdb.gz"),
         target = str_remove(string = target, pattern = "-F[0-9][0-9][0-9]-model_v4.pdb.gz"),
         target = str_remove(string = target, pattern = "AF-"))

human_uniprot_metadata <- read_tsv(args$input_human_metadata) %>%
  clean_names() %>%
  rename_with(.cols = everything(), function(x){paste0("human_", x)})

query_uniprot_metadata <- read_tsv(args$input_query_metadata) %>%
  clean_names() %>%
  rename_with(.cols = everything(), function(x){paste0("query_", x)})

foldseek_results <- foldseek_results %>%
  left_join(human_uniprot_metadata, by = c("target" = "human_protid")) %>%
  left_join(query_uniprot_metadata, by = c("query" = "query_protid")) %>%
  dplyr::relocate(query_organism) %>%
  dplyr::relocate(all_of(c("human_gene_names_primary", "human_function_cc",
                           "human_tissue_specificity", "human_subcellular_location_cc")),
                  .before = lddt) %>%
  arrange(desc(alntmscore))

# foldseek_results_top_matches <- foldseek_results %>%
#   group_by(query) %>%
#   slice_max(qtmscore) %>%
#   filter(qtmscore > 0.7 | alntmscore > 0.7 | ttmscore > 0.7)

write_tsv(foldseek_results, args$output)