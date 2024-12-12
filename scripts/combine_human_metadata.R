library(tidyverse)
library(janitor)
library(optparse)


# set CLI arguments -------------------------------------------------------

option_list <- list(
  make_option(c("--input_uniprot"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--input_protein_atlas"), type="character",
              help="Path to human uniprot metadata TSV file."),
  make_option(c("--input_mouse_opentargets"), type="character",
              help="Path to query uniprot metadata CSV file."),
  make_option(c("--output"), type="character",
              help="Path to output CSV file."),
  make_option(c("--output_filtered"), type="character",
              help="Path to output CSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

# Define sets of key words for protein location ---------------------------
# Also create a regular expression pattern that includes all words
extracellular_keywords <- c('extracellular', 'secreted', 'cell surface',
                           'pericellular', 'secretory vesicle membrane',
                           'cytoplasmic vesicle, secretory vesicle membrane',
                           'synaptic vesicle membrane')
extracellular_pattern <- paste(extracellular_keywords, collapse = "|")

molecular_function_keywords <- c("receptor", "extracellular", "secreted",
                                 "cytokine", "neuropeptide", "hormone")
molecular_function_pattern <- paste(molecular_function_keywords, collapse = "|")

immune_expression_keywords <- c("innate immune response", "immune response",
                                "adaptive immune response", "cytokine signaling",
                                "antigen presentation", "t-cell receptor", "b-cell",
                                "nk-cell", "monocyte", "neutrophil", "macrophage",
                                "plasma cell", "dendritic cell", "lymphoid tissue",
                                "spleen", "thymus")
immune_keywords_pattern <- paste(immune_expression_keywords, collapse = "|")

# Note that "immunity" will capture innate immunity, adaptive immunity, and
# immunity
biological_process_keywords <- c("inflammatory response", "immunity",
                                 "host-virus interaction")
biological_process_pattern <- paste(biological_process_keywords, collapse = "|")

# Some genes that are immune involved are left out. The three patterns below
# rescue some of these.
gene_description_keywords <- c("interleukin", "immunoglobulin", "chemokine")
gene_description_pattern <- paste(gene_description_keywords, collapse = "|")

protein_class_keywords <- c("cd marker")
protein_class_pattern <- paste(protein_class_keywords, collapse = "|")

specific_immune_involved_proteins <- c("ATRNL1", "ERBB2", "GPR83", "HTRA3",
                                       "LGALS3BP", "PRSS57", "RIPK4", "SERPINA3",
                                       "SERPINE1", "SERPINE2", "SRPB1", "TYRO3",
                                        "VSTM5")
specific_immune_involved_pattern <- paste(specific_immune_involved_proteins, collapse = "|")

# Read in metadata files and combine --------------------------------------
uniprot <- read_tsv(args$input_uniprot) %>%
  clean_names() %>%
  rename_with(~ paste0("uniprot_", .), everything()) %>%
  rename(uniprot = uniprot_protid)

protein_atlas <- read_tsv(args$input_protein_atlas) %>% 
  clean_names() %>%
  rename_with(~ paste0("proteinatlas_", .), everything()) %>%
  rename(uniprot = proteinatlas_uniprot, ensembl = proteinatlas_ensembl) %>%
  relocate(uniprot, ensembl, .before = proteinatlas_gene)

mouse_phenotypes <- read_csv(args$input_mouse_opentargets) %>%
  rename(opentargets_mouse_phenotype = mouse_phenotype)

metadata <- uniprot %>%
  left_join(protein_atlas, by = "uniprot") %>%
  left_join(mouse_phenotypes, by = "ensembl") %>%
  # make columns lowercase that we'll be searching with regular expressions
  mutate(across(
    c(uniprot_subcellular_location_cc, proteinatlas_secretome_location,
      proteinatlas_molecular_function, proteinatlas_blood_expression_cluster,
      proteinatlas_tissue_expression_cluster,
      proteinatlas_single_cell_expression_cluster,
      proteinatlas_biological_process, proteinatlas_gene_description,
      proteinatlas_protein_class),
    tolower))

# process metadata to pull out information of interest --------------------

# Determine extracellular location
metadata <- metadata %>%
  mutate(
    any_extracellular_location = 
      if_else(
        str_detect(uniprot_subcellular_location_cc, extracellular_pattern) |
          str_detect(proteinatlas_secretome_location, 'secreted') |
          str_detect(proteinatlas_molecular_function, molecular_function_pattern) |
          str_detect(uniprot_signal_peptide,  "SIGNAL"),
        'extracellular', 'not extracellular')
  ) %>%
  mutate(any_extracellular_location = ifelse(is.na(any_extracellular_location), "not extracellular", any_extracellular_location))

# Annotate whether there is evidence of expression of the gene in
# immune-relevant cells.
metadata <- metadata %>%
  mutate(
    immune_involvement = if_else(
      str_detect(proteinatlas_blood_expression_cluster, immune_keywords_pattern) |
        str_detect(proteinatlas_tissue_expression_cluster, immune_keywords_pattern) |
        str_detect(proteinatlas_single_cell_expression_cluster, immune_keywords_pattern) |
        str_detect(proteinatlas_biological_process, biological_process_pattern) |
        str_detect(proteinatlas_gene_description, gene_description_pattern) |
        str_detect(proteinatlas_protein_class, protein_class_pattern) |
        str_detect(uniprot_gene_names_primary, specific_immune_involved_pattern),
      'immune cell/tissue expression or involved in inflammation/immunity',
      'no immune cell/tissue expression or involved in inflammation/immunity')
  ) %>%
  mutate(
    immune_involvement = ifelse(
      is.na(immune_involvement),
      "no immune cell/tissue expression or involved in inflammation/immunity",
      immune_involvement)
  )

# Annotate genes that have an immune phenotype in mice
metadata <- metadata %>%
  mutate(
    mouse_immune_phenotype = if_else(
      str_detect(opentargets_mouse_phenotype, 'immune system phenotype'),
      'mouse immune phenotype', 'no mouse immune phenotype')
    ) %>%
  mutate(mouse_immune_phenotype = ifelse(is.na(mouse_immune_phenotype), "no mouse immune phenotype", mouse_immune_phenotype))

write_csv(metadata, args$output)

filtered_metadata <- metadata %>%
  filter(mouse_immune_phenotype == "mouse immune phenotype") %>%
  filter(immune_involvement == "immune cell/tissue expression or involved in inflammation/immunity") %>%
  filter(any_extracellular_location == "extracellular")

write_csv(filtered_metadata, args$output_filtered)
