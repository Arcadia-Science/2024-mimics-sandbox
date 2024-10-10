library(openxlsx)
library(tidyverse)
library(janitor)

option_list <- list(
  make_option(c("--nomburg"), type="character",
              help="Path human viruses TSV file."),
  make_option(c("--human_viruses"), type="character",
              help="Path human viruses TSV file."),
  make_option(c("--output_structure_filepaths"), type="character",
              help="Path to output TXT file with relative structure filepaths in nomburg structure zip archive."),
  make_option(c("--output_tsv", type = "character",
                help="Path to output TSV file with combined metadata."))
)

args <- parse_args(OptionParser(option_list=option_list))

supp_table1 <- read.xlsx(args$nomburg, sheet = 1)

# Parse viral protein metadata information from protein names.
supp_table1_parsed <- supp_table1 %>%
  # Fix row that causes parsing error.
  mutate(cluster_member = gsub("putative_membrane__secreted_protein", "putative_membrane_secreted_protein", cluster_member)) %>%
  separate(col = cluster_member, into = c("gene_name", "ncbi_id", "virus", "taxon_id"), sep = "__", remove = FALSE)

# Read in metadata for viruses that infect humans.
viralzone_human_viruses <- read_tsv(args$human_viruses) %>%
  rename_with(.cols = everything(), function(x){paste0("viralzone_", x)}) %>%
  rename(viralzone_link = viralzone_viralzone_link)

# Read in and format UniProt metadata 
uniprot_metadata <- read_tsv("sandbox/try_viral_foldseek_db/human_virus_uniprot_features.tsv") %>%
  clean_names() %>%
  # Remove any rows where refseq does not have a value.
  filter(!is.na(ref_seq)) %>%
  # Some UniProt accessions correspond to multiple RefSeq accessions.
  # Separate these out so we can match between the Nomburg identifiers and the
  # UniProt <-> NCBI accession map.
  separate(col = ref_seq, into = c("ref_seq1", "ref_seq2"), sep = ";", remove = F) %>%
  # Some accessions have additional information reported in brackets; remove this.
  mutate(ref_seq1 = str_remove(string = ref_seq1, pattern = " \\[.*\\]")) %>%
  # Nomburg removed accession versions, so remove these here as well.
  mutate(ref_seq1 = str_remove(pattern = "\\.[0-9]$", string = ref_seq1)) %>%
  mutate(ref_seq2 = str_remove(pattern = "\\.[0-9]$", string = ref_seq2)) %>%
  # If there were multiple RefSeq IDs associated with a UniProt ID, choose the
  # one that Nomburg used.
  mutate(ref_seq_join = ifelse(ref_seq1 %in% supp_table1_parsed$ncbi_id, ref_seq1, ref_seq2)) %>%
  select(-ref_seq1, -ref_seq2)

# join information together and output structure filepaths in zip archive
human_viruses_with_structures <- uniprot_metadata %>%
  left_join(supp_table1_parsed, by = c("ref_seq_join" = "ncbi_id")) %>%
  left_join(viralzone_human_viruses, by = c("taxon_id" = "genome_assembly_ncbi_taxon_id_aligned_with_nomburg")) %>%
  # remove TaxonID as it is redundant, and remove subcluster_rep bc it isn't
  # clearly described in Nomburg et al.
  select(-taxonID, -subcluster_rep) %>%
  # Rename the cluster columns to something more intuitive, since it won't be
  # clear what these are outside of the context of the supplementary table.
  rename(nomburg_cluster_id = cluster_ID, nomburg_cluster_representative = cluster_rep,
         nomburg_cluster_count = cluster_count, nomburg_protein_name = cluster_member) %>%
  # generate the filepaths for the structures in the zip archive
  mutate(structure_filepaths = paste0("viral_structures", "/", family, "/", nomburg_protein_name, ".pdb"))

# Write out relative structure filepaths that will be used to selective decompress
# structures from zip archive of all viral structures in db.
write_tsv(human_viruses_with_structures %>% pull(structure_filepaths),
          file = args$output_structure_filepaths,
          col_names = FALSE)

# Write out combined metadata
write_tsv(human_viruses_with_structures,
          file = args$output_tsv)