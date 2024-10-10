library(openxlsx)
library(tidyverse)
library(janitor)
library(optparse)


option_list <- list(
  make_option(c("--nomburg"), type="character",
              help="Path human viruses TSV file."),
  make_option(c("--human_viruses"), type="character",
              help="Path human viruses TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output NCBI identifier TXT file.")
)

args <- parse_args(OptionParser(option_list=option_list))

# Read in supplementary table from Nomburg, J., Doherty, E.E., Price, N. et al. 
# Birth of protein folds and functions in the virome. Nature 633, 710–717 (2024). 
# https://doi.org/10.1038/s41586-024-07809-y
# Table includes names and taxonimic identifiers for viral proteins in the database.
supp_table1 <- read.xlsx(args$nomburg, sheet = 1)

# Parse viral protein metadata information from protein names.
supp_table1_parsed <- supp_table1 %>%
  # Fix row that causes parsing error.
  mutate(cluster_member = gsub("putative_membrane__secreted_protein", "putative_membrane_secreted_protein", cluster_member)) %>%
  separate(col = cluster_member, into = c("gene_name", "ncbi_id", "virus", "taxon_id"), sep = "__", remove = FALSE)

# Read in metadata for viruses that infect humans.
# Note that 16 viruses in this table are not represented in the Nomburg database.
viralzone_human_viruses <- read_tsv(args$human_viruses)

# Filter to proteins from viruses that infect humans.
# There are likely additional proteins in Nomburg db that originate from viruses 
# that infect humans. However, we ignore these for this first pass to focus on 
# bonafide cases that we think will produces the strongest signals of structural
# mimicry.
supp_table1_parsed_filtered <- supp_table1_parsed %>%
  filter(taxon_id %in% viralzone_human_viruses$genome_assembly_ncbi_taxon_id_aligned_with_nomburg) 

# Write out NCBI accessions for viral proteins that infect humans and that are
# represented in the Nomburg database. This information will be used to retrieve
# UniProt protein IDs and metadata. Note that not all proteins on this list will
# have a UniProt accession, so in requiring UniProt metadata, we are further
# excluding some proteins that are from viruses that infect humans.
write_tsv(supp_table1_parsed_filtered %>% select(ncbi_id),
          file = args$output,
          col_names = FALSE)