library(tidyverse)
library(janitor)
setwd("~/github/2024-parasite-human-mimics/")


# read in results ---------------------------------------------------------
# results <- read_csv("~/Downloads/per_cluster_gmm_results030725.csv")
# details <- read_csv("~/Downloads/best_cluster_rows_030725_detailed.csv")
results <- read_csv("~/Downloads/gmmviro3d_full032525_evalue.csv")
details <- read_tsv("~/Downloads/gmmviro3d_full032525_evalue_all_details.tsv")

# summarize features ------------------------------------------------------

summarize_feature <- function(data, feature) {
  data %>%
    mutate(has_feature = !is.na(.[[feature]])) %>%
    group_by(original_cluster_id, best_cluster) %>%
    summarise(total = n(),
              with_feature = sum(has_feature, na.rm = TRUE),
              ratio = paste0(with_feature, "/", total),
              fraction = with_feature / total,
              .groups = 'drop') %>%
    rename(!!paste0(feature, "_ratio") := ratio,
           !!paste0(feature, "_fraction") := fraction) %>%
    select(-total, -with_feature) # Optionally remove the intermediate columns if not needed
}

host_transmembrane <- summarize_feature(details, "host_transmembrane")
host_glycosylation <- summarize_feature(details, "host_glycosylation")
host_disulfide_bond <- summarize_feature(details, "host_disulfide_bond")
host_lipidation <- summarize_feature(details, "host_lipidation")
host_signal_peptide <- summarize_feature(details, "host_signal_peptide")
host_peptide <- summarize_feature(details, "host_peptide")
host_post_translational_modification <- summarize_feature(details, "host_post_translational_modification")
host_propeptide <- summarize_feature(details, "host_propeptide")
host_transit_peptide <- summarize_feature(details, "host_transit_peptide")


query_transmembrane <- summarize_feature(details, "query_transmembrane")
query_glycosylation <- summarize_feature(details, "query_glycosylation")
query_disulfide_bond <- summarize_feature(details, "query_disulfide_bond")
query_lipidation <- summarize_feature(details, "query_lipidation")
query_signal_peptide <- summarize_feature(details, "query_signal_peptide")
query_peptide <- summarize_feature(details, "query_peptide")
query_post_translational_modification <- summarize_feature(details, "query_post_translational_modification")
query_propeptide <- summarize_feature(details, "query_propeptide")
query_transit_peptide <- summarize_feature(details, "query_transit_peptide")


# summarize character features --------------------------------------------

# Define a function to summarize character features based on a "no " prefix
summarize_character_feature <- function(data, feature) {
  data %>%
    mutate(is_negative = str_starts(.[[feature]], "no "),
           is_positive = !is.na(.[[feature]]) & !is_negative) %>%
    group_by(original_cluster_id, best_cluster) %>%
    summarise(total = n(),
              positive_count = sum(is_positive, na.rm = TRUE),
              positive_ratio = paste0(positive_count, "/", total),
              positive_fraction = positive_count / total,
              .groups = 'drop') %>%
    rename(!!paste0(feature, "_ratio") := positive_ratio,
           !!paste0(feature, "_fraction") := positive_fraction) %>%
    select(-total, -positive_count)
}

host_immune_involvement <- summarize_character_feature(details, "immune_involvement")
host_mouse_immune_phenotype <- summarize_character_feature(details, "mouse_immune_phenotype")

# averages ----------------------------------------------------------------

calculate_average <- function(data, feature) {
  data %>%
    group_by(original_cluster_id, best_cluster) %>%
    summarise(!!paste0(feature, "_average") := mean(get(feature), na.rm = TRUE),
              .groups = 'drop')
}

average_tlen <- calculate_average(details, "host_length")
average_qlen <- calculate_average(details, "query_length")
host_average_plddt <- calculate_average(details, "host_pdb_plddt")
host_average_annotation <- calculate_average(details, "host_annotation")

query_plddt_info <- details %>%
  select(original_cluster_id, best_cluster, query, query_chosen_method, query_colab_fold_p_lddt, query_esm_fold_p_lddt) %>%
  mutate(query_plddt = ifelse(query_chosen_method == "ColabFold", query_colab_fold_p_lddt, query_esm_fold_p_lddt)) 

query_average_plddt <- calculate_average(query_plddt_info, "query_plddt")
query_average_annotation <- calculate_average(details, "query_annotation")

# location ----------------------------------------------------------------

analyze_active_sites <- function(data) {
  data %>%
    mutate(active_sites = str_extract_all(host_active_site, "\\d+"),
           active_sites = map(active_sites, ~ifelse(length(.x) == 0, NA, as.numeric(.x))),
           has_active_site = map_lgl(active_sites, ~any(!is.na(.x)))) %>%
    rowwise() %>%
    mutate(in_window = sum(map_lgl(active_sites, ~ .x >= tstart & .x <= tend), na.rm = TRUE),
           total_sites = max(length(unlist(active_sites)), 1)) %>%
    group_by(original_cluster_id, best_cluster) %>%
    summarise(total_proteins = n(),
              proteins_with_active_sites = sum(has_active_site, na.rm = TRUE),
              total = proteins_with_active_sites,
              with_feature = sum(in_window),
              ratio = paste0(with_feature, "/", total),
              fraction = with_feature / total,
              .groups = 'drop') %>%
    mutate(ratio = ifelse(is.na(ratio), "0/0", ratio),
           fraction = replace_na(fraction, 0)) %>%
    rename(host_active_site_in_viral_aln_ratio = ratio,
           host_active_site_in_viral_aln_fraction = fraction,
           host_active_site_presence = proteins_with_active_sites) %>%
    select(-total, -with_feature)  %>%
    mutate(host_active_site_in_viral_aln_ratio = ifelse(host_active_site_presence == 0, "no host active site", host_active_site_in_viral_aln_ratio),
           host_active_site_in_viral_aln_fraction = ifelse(host_active_site_presence == 0, NA, host_active_site_in_viral_aln_fraction))
}


active_site_analysis <- analyze_active_sites(details)



# concatenation -----------------------------------------------------------

concatenate_feature <- function(data, feature) {
  data %>%
    group_by(original_cluster_id, best_cluster) %>%
    summarise(summary = paste(unique(na.omit(get(feature))), collapse = "; "),
              .groups = 'drop') %>%
    rename_with(~paste0(feature, "_summary"), .cols = "summary")
}

host_ec_number <- concatenate_feature(details, "host_ec_number")
host_interacts_with <- concatenate_feature(details, "host_interacts_with")
host_single_cell_expression_cluster <- concatenate_feature(details, "proteinatlas_single_cell_expression_cluster")
host_rna_single_cell_type_specific_n_tpm <- concatenate_feature(details, "proteinatlas_rna_single_cell_type_specific_n_tpm")
host_rna_tissue_specific_n_tpm <- concatenate_feature(details, "proteinatlas_rna_tissue_specific_n_tpm")

#query_interacts_with <- concatenate_feature(details, "query_interacts_with")

query_family <- concatenate_feature(details, "query_family")
query_interacts_with <- concatenate_feature(details, "query_interacts_with")

# pull out bits -----------------------------------------------------------

detect_extracellular_from_topo_domain <- function(data, feature) {
  data %>%
    group_by(original_cluster_id, best_cluster) %>%
    summarise(
      total = n(),
      contains_extracellular = sum(str_detect(get(feature), "Extracellular"), na.rm = TRUE),
      ratio = paste0(contains_extracellular, "/", total),
      fraction = contains_extracellular / total,
      .groups = 'drop'
    ) %>%
    rename(
      !!paste0(feature, "_extracellular_ratio") := ratio,
      !!paste0(feature, "_extracellular_fraction") := fraction
    ) %>%
    select(-total, -contains_extracellular) 
}

host_any_topo_domain_extracellular_count <- detect_extracellular_from_topo_domain(details, "host_topological_domain")



# combine summary data ----------------------------------------------------


combined_data <- host_transmembrane %>%
  left_join(query_transmembrane, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_glycosylation, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_glycosylation, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_disulfide_bond, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_disulfide_bond, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_lipidation, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_lipidation, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_signal_peptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_signal_peptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_peptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_peptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_post_translational_modification, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_post_translational_modification, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_propeptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_propeptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_transit_peptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_transit_peptide, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_immune_involvement, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_mouse_immune_phenotype, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(average_tlen, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(average_qlen, by = c("original_cluster_id", "best_cluster")) %>%
  mutate(length_ratio = query_length_average / host_length_average) %>%
  select(-query_length_average, -host_length_average) %>%
  left_join(host_average_plddt, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_average_plddt, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_average_annotation, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_average_annotation, by = c("original_cluster_id", "best_cluster")) %>%
  #left_join(active_site_analysis, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_ec_number, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_interacts_with, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_interacts_with, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_single_cell_expression_cluster, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_rna_single_cell_type_specific_n_tpm, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_rna_tissue_specific_n_tpm, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(query_family, by = c("original_cluster_id", "best_cluster")) %>%
  left_join(host_any_topo_domain_extracellular_count, by = c("original_cluster_id", "best_cluster"))




combined_data <- combined_data %>%
  mutate(transmembrane_fractional_change = query_transmembrane_fraction - host_transmembrane_fraction,
         glycosylation_fractional_change = query_glycosylation_fraction - host_glycosylation_fraction,
         disulfide_bond_fractional_change = query_disulfide_bond_fraction - host_disulfide_bond_fraction,
         lipidation_fractional_change = query_lipidation_fraction - host_lipidation_fraction,
         signal_peptide_fractional_change = query_signal_peptide_fraction - host_signal_peptide_fraction,
         peptide_fractional_change = query_peptide_fraction - host_peptide_fraction,
         post_translational_modification_fractional_change = query_post_translational_modification_fraction - host_post_translational_modification_fraction,
         propeptide_fractional_change = query_propeptide_fraction - host_propeptide_fraction,
         transit_peptide_fractional_change = query_transit_peptide_fraction - host_transit_peptide_fraction) %>%
  select(original_cluster_id, best_cluster, length_ratio, 
         transmembrane_fractional_change, glycosylation_fractional_change,
         disulfide_bond_fractional_change, lipidation_fractional_change,
         signal_peptide_fractional_change, peptide_fractional_change, 
         post_translational_modification_fractional_change,
         propeptide_fractional_change, transit_peptide_fractional_change,
         immune_involvement_fraction,
         mouse_immune_phenotype_fraction, host_pdb_plddt_average, query_plddt_average, 
         host_annotation_average, query_annotation_average, 
         #host_active_site_in_viral_aln_fraction,
         host_ec_number_summary, 
         host_interacts_with_summary, query_interacts_with_summary,
         proteinatlas_single_cell_expression_cluster_summary, 
         proteinatlas_rna_single_cell_type_specific_n_tpm_summary, 
         proteinatlas_rna_tissue_specific_n_tpm_summary, 
         query_family_summary, 
         host_topological_domain_extracellular_fraction)


count_shared_interactions <- function(data) {
  data %>%
    mutate(
      host_interacts_with_list = str_split(host_interacts_with_summary, ";\\s*"),
      query_interacts_with_list = str_split(query_interacts_with_summary, ";\\s*"),
      num_interacts_with_shared_between_query_and_host = map2_int(host_interacts_with_list, query_interacts_with_list, 
                                                                  ~length(intersect(
                                                                    .x[.x != ""],
                                                                    .y[.y != ""]
                                                                  )))
    ) %>%
    select(-host_interacts_with_list, -query_interacts_with_list)
}

combined_data <- count_shared_interactions(combined_data) %>%
  relocate(num_interacts_with_shared_between_query_and_host, .after = "query_interacts_with_summary")


# switch from uniprot ids to gene symbols ---------------------------------


replace_uniprot_with_gene_names <- function(summary_column, mapping_df) {
  str_split(summary_column, ";\\s*") %>%  # Split each summary by semicolons
    map_chr(function(ids_list) {
      # Extract potential UniProt IDs (including those within brackets)
      ids_cleaned <- str_extract_all(ids_list, "[A-Z0-9-]+(?=\\]|$)") %>% 
        unlist() %>% 
        na.omit() %>% 
        unique() %>%
        map_chr(~ str_remove(.x, "-.*"))
      
      # Map IDs to gene names using the mapping dataframe
      ids_mapped <- map_chr(ids_cleaned, function(id) {
        gene_name <- mapping_df$uniprot_gene_names_primary[mapping_df$uniprot == id]
        if (length(gene_name) == 0) id else gene_name  # Return the original id if no mapping is found
      })
      
      # Reconstruct the string with semicolons
      paste(ids_mapped, collapse = "; ")
    })
}

human_metadata <- read_csv("outputs/metadata/human_metadata_combined.csv") %>%
  select(uniprot, uniprot_gene_names_primary)

virus_metadata <- read_tsv("outputs/viral/human/viral_protein_features.tsv") %>%
  clean_names() %>%
  select(uniprot = entry, uniprot_gene_names_primary = gene_names_primary)

uniprot_metadata <- bind_rows(human_metadata, virus_metadata) %>%
  distinct()


combined_data$host_interacts_with_summary <- replace_uniprot_with_gene_names(combined_data$host_interacts_with_summary,
                                                                             mapping_df = uniprot_metadata)
combined_data$query_interacts_with_summary <- replace_uniprot_with_gene_names(combined_data$query_interacts_with_summary,
                                                                              mapping_df = uniprot_metadata)



# cluster summary ---------------------------------------------------------



final <- left_join(results, combined_data, by = c("original_cluster_id", "best_cluster"))
write_csv(final, "~/Downloads/gmmviro3d_full032525_evalue_with_summary.csv")
