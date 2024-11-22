library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--host"), type="character",
              help="Name of host organism."),
  make_option(c("--input"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output CSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

foldseek_results <- read_tsv(args$input, show_col_types = FALSE) 

foldseek_results_filtered <- foldseek_results %>%
  # filter on quality and strength of foldseek match
  filter(tcov > 0.25) %>%
  filter(qcov > 0.25) %>%
  filter(alnlen > 100) %>%
  filter(evalue < 0.0001) %>%
  filter(host_pdb_plddt > 50) %>%
  filter(query_pdb_plddt > 50) %>%
  filter(lddt > 0.5) %>%
  rowwise() %>%
  mutate(min_tmscore = min(alntmscore, qtmscore, ttmscore),
         max_tmscore = max(alntmscore, qtmscore, ttmscore),
         avg_tmscore = sum(alntmscore, qtmscore, ttmscore) / 3 ) %>%
  ungroup() %>%
  filter(max_tmscore > 0.4) %>%
  dplyr::relocate(all_of(c("avg_tmscore", "min_tmscore")),
                  .after = ttmscore) %>%
  # remove viral matches that hit DNA/RNA replicative machinery
  filter(!grepl(pattern = "polymerase", x = host_function_cc))

if(args$host == "human"){
  foldseek_results_filtered <- foldseek_results_filtered %>%
    # filter to host proteins that are extracellular 
    filter(any_extracellular_location == "extracellular") %>%
    # bring forward the immune and mouse information
    dplyr::relocate(immune_involvement, mouse_immune_phenotype, 
                    .after = host_subcellular_location_cc)
} else {
  foldseek_results_filtered <- foldseek_results_filtered %>%
    # filter to host proteins that have signal peptides OR proteins that are
    # annotated as extracellular or membrane targets
    mutate(host_subcellular_location_cc = tolower(host_subcellular_location_cc)) %>%
    filter(
      grepl(pattern = "SIGNAL", x = host_signal_peptide) |
        grepl(pattern = "cell membrane|cell junction|cell projection|subcellular location: membrane|secreted",
              x = host_subcellular_location_cc)
    )
}

foldseek_results_filtered <- foldseek_results_filtered %>%
  # change order of columns so it's easier to see the host organism
  dplyr::relocate(host_organism, .before = target) %>%
  # select only the top hit for each query
  arrange(query, evalue) %>%
  group_by(query, query_species) %>%
  slice_min(evalue) %>%
  # sometimes, there are two hits with the same evalue
  # select the one with the highest max tmscore
  slice_max(max_tmscore) 

write_csv(foldseek_results_filtered, args$output)
