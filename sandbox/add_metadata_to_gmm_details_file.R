library(tidyverse)
library(janitor)
library(optparse)

option_list <- list(
  make_option(c("--host"), type="character",
              help="Name of host organism."),
  make_option(c("--input_results"), type="character",
              help="Path to details CSV file from GMM"),
  make_option(c("--input_human_metadata"), type="character",
              help="Path to human metadata CSV file."),
  make_option(c("--input_host_metadata"), type="character",
              help="Path to host uniprot metadata TSV file."),
  make_option(c("--input_host_lddt"), type="character",
              help="Path to host structure quality measurement TSV file."),
  make_option(c("--input_query_metadata"), type="character",
              help="Path to query metadata TSV file."),
  make_option(c("--input_query_uniprot_metadata"), type="character",
              help="Path to query metadata TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

# args$host <- "human"
# args$input_results <- "~/Downloads/gmmviro3d_full032525_evaluedetail.csv"
# args$input_human_metadata <- "outputs/metadata/human_metadata_combined.csv"
# args$input_host_metadata <- "outputs/viral/human/host_proteome_protein_features.tsv"
# args$input_host_lddt <- "outputs/viral/human/host_proteome_pdb_structure_quality.tsv"
# args$input_query_metadata <- "inputs/viral/merged_viral_metadata_human.tsv"
# args$input_query_uniprot_metadata <- "outputs/viral/human/viral_protein_features.tsv"
# args$output <- "~/Downloads/gmmviro3d_full032525_evalue_all_details.tsv"


results <- read_csv(args$input_results, show_col_types = FALSE)

if(nrow(results) > 0){

  host_pdb_lddt <- read_tsv(args$input_host_lddt, show_col_types = FALSE) %>%
    select(protid, pdb_plddt = pdb_confidence) %>%
    distinct()
  
  query_uniprot_metadata <- read_tsv(args$input_query_uniprot_metadata, show_col_types = FALSE) %>%
    clean_names() %>%
    select(-protid) %>%
    distinct() %>%
    group_by(entry) %>%
    slice_head(n = 1)
  
  query_metadata <- read_tsv(args$input_query_metadata, show_col_types = FALSE) %>%
    clean_names() %>%
    rename(uniprot_id = uni_prot_id) %>%
    left_join(query_uniprot_metadata, by = c("uniprot_id" = "entry"), relationship = "many-to-one") %>%
    # names might not match, so edit to original file names
    mutate(protid = str_remove(string = structure_file, pattern = "\\.pdb$")) %>%
    rename_with(.cols = everything(), function(x){paste0("query_", x)}) %>%
    distinct()
  
  # We have more metadata for human (UniProt, Protein Atlas, OpenTargets) than we
  # do for other species (UniProt). Treat this differently but try to end up with
  # the same column names when the same information is represented (UniProt)
  if(args$host == "human"){
    host_metadata <- read_csv(args$input_human_metadata, show_col_types = FALSE) %>% 
      rename(protid = uniprot) %>%
      left_join(host_pdb_lddt, by = c("protid")) %>%
      # We want the metadata columns that are the same between species to have the
      # same column names.
      rename(host_protid = protid, host_pdb_plddt = pdb_plddt) %>%
      rename_with(.cols = starts_with("uniprot"), function(x){gsub("uniprot_", "host_", x)}) %>%
      group_by(host_protid) %>%
      slice_head(n = 1) %>%
      distinct()
  } else {
    host_metadata <- read_tsv(args$input_host_metadata, show_col_types = FALSE) %>%
      clean_names() %>%
      left_join(host_pdb_lddt, by = c("protid")) %>%
      rename_with(.cols = everything(), function(x){paste0("host_", x)}) %>%
      distinct()
  }
  
  results <- results %>%
    left_join(host_metadata, by = c("target" = "host_protid")) %>%
    left_join(query_metadata, by = c("query" = "query_protid")) %>%
    dplyr::relocate(query_virus_name) %>%
    dplyr::relocate(all_of(c("host_gene_names_primary",
                             "host_protein_names")),
                    .after = protein_id) %>%
    arrange(evalue) %>%
    distinct()
  
  write_tsv(results, args$output)

} else {
  file.create(args$output_full)
  file.create(args$output)
}

