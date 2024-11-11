library(tidyverse)
library(janitor)
library(optparse)

option_list <- list(
  make_option(c("--input_foldseek_results"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--input_host_metadata"), type="character",
              help="Path to host uniprot metadata TSV file."),
  make_option(c("--input_host_lddt"), type="character",
              help="Path to host structure quality measurement TSV file."),
  make_option(c("--input_query_metadata"), type="character",
              help="Path to query metadata TSV file."),
  make_option(c("--input_query_lddt"), type="character",
              help="Path to query structure quality measurement TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

#foldseek_results <- read_tsv("outputs/viral/rat/nomburg_virus_matches.tsv") %>%
foldseek_results <- read_tsv(args$input_foldseek_results, show_col_types = FALSE) %>%
  mutate(query = str_remove(string = query, pattern = "\\.pdb"),
         target = str_remove(string = target, pattern = "\\.pdb"))


#host_pdb_lddt <- read_tsv("outputs/viral/rat/host_proteome_pdb_structure_quality.tsv") %>%
host_pdb_lddt <- read_tsv(args$input_host_lddt, show_col_types = FALSE) %>%
  select(protid, pdb_plddt = pdb_confidence)

#host_uniprot_metadata <- read_tsv("outputs/viral/rat/host_proteome_protein_features.tsv") %>%
host_uniprot_metadata <- read_tsv(args$input_host_metadata, show_col_types = FALSE) %>%
  clean_names() %>%
  left_join(host_pdb_lddt, by = c("protid")) %>%
  rename_with(.cols = everything(), function(x){paste0("host_", x)})

#query_pdb_lddt <- read_tsv("outputs/viral/rat/viral_structure_quality.tsv") %>%
query_pdb_lddt <- read_tsv(args$input_query_lddt, show_col_types = FALSE) %>%
  select(protid, pdb_plddt = pdb_confidence)

#query_metadata <- read_tsv("outputs/viral/rat/viral_structure_metadata.tsv") %>%
query_metadata <- read_tsv(args$input_query_metadata, show_col_types = FALSE) %>%
  clean_names() %>%
  left_join(query_pdb_lddt, by = c("nomburg_protein_name" = "protid")) %>%
  rename_with(.cols = everything(), function(x){paste0("query_", x)}) %>%
  # Some NCBI accessions had multiple uniprot accessions.
  # Below we filter to keep the one with the most complete annotations (or
  # arbitrarily select one if they're equally annotated) so that our results
  # aren't duplicated upon joining to the metadata.
  mutate(num_NAs = rowSums(is.na(across(.cols = -query_nomburg_protein_name)))) %>%
  group_by(query_nomburg_protein_name) %>%
  filter(num_NAs == min(num_NAs)) %>%
  # If there are ties, keep the first one
  slice(1) %>%
  ungroup() %>%
  select(-num_NAs, query_structure_filepaths)
  

foldseek_results <- foldseek_results %>%
  left_join(host_uniprot_metadata, by = c("target" = "host_protid")) %>%
  left_join(query_metadata, by = c("query" = "query_nomburg_protein_name")) %>%
  dplyr::relocate(query_species) %>%
  dplyr::relocate(all_of(c("host_gene_names_primary", "host_function_cc",
                           "host_tissue_specificity", "host_subcellular_location_cc")),
                  .before = lddt) %>%
  arrange(desc(alntmscore))

write_tsv(foldseek_results, args$output)
