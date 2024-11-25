# This script selects viral mimics according to criteria that the translation
# pilot team came up with at our October 2024 translation planning meeting.
# Our aim is to select the viral mimics that we think are the best and easiest
# candidates for being a drug or identifying a drug target in the
# anti-inflammatory space. To do this, we focused on identifying hits that are:
# 1. Likely to be real matches with a strong signature of similarity to the 
#    structure of the host protein sequence. These 
#    a. tcov, qcov, evalue: we set these values based on intuition of what we
#       anticipate real hits will looks like.
#    b. alnlen: we set alnlen to 100 to try to capture viral hits that are
#       likely real mimics acquired by HGT from the host (or a closely related
#       organism). The exact length of 100 was based on intuition.
#    c. max_tmscore: We decided to use a cut off of 0.4 after observing that
#       when we use our viral structural mimicry detection approach on viral
#       mimics in the literature, the lowest TM Score for a correct hit was
#       0.49. We buffered this slightly to 0.4 given that we only have a handful
#       of these examples. We define viral mimics in the literature as those
#       with experimentally confirmed functions and often with 1) solved crystal
#       structures that demonstrate mimicry of host proteins and 2) evolutionary
#       analysis of where the mimic originated (HGT).
# 2. High quality matches from high quality structural predictions. We limit the
#    each pLDDT to > 50. AlphaFold labels structural predictions with a pLDDT <
#    50 as low-confidence predictions.
#    a. host_pdb_plddt
#    b. query_pdb_plddt
#    c. lddt
# 3. Not something the virus uses to replicate itself. We think these proteins
#    may hav interesting properties, but for this first application we're
#    focused on we don't think these proteins have the highest chance of
#    helping us find ways to reduce inflammation.
#    a. Right now, we're only eliminating viral proteins that mimic host genes
#       that have the word "polymerase" in their description. This likely leaves
#       behind some replicative machinery, but our other filters seem to remove
#       other proteins in this category well enough for our use case.
# 4. The human protein that is mimicked is easy reachable by a drug (e.g. the
#    protein has an extracellular or membrane annotation of some kind).
#    a. host_signal_peptide: UniProt annotates proteins that contain a signal
#       peptide. Signal peptides target direct proteins to different parts of
#       the cell and function as a secretion signal.
#    b. host_subcellular_location_cc: UniProt annotates the subcellular location
#       of proteins. If the subcellular location is annotated as cell membrane
#       other non-intracellular locations, we keep the protein.
# 5. Selecting the best structural mimic.
#    a. evalue: when a viral protein has matches to multiple host proteins, we
#       keep the viral hit with the lowest evalue. We decided to use this
#       filtering approach after observing that when we use our viral structural
#       mimicry detection approach on viral mimics in the literature, the
#       the correct hit was most frequently the one with the lowest e-value. We
#       define viral mimics in the literature as those with experimentally
#       confirmed functions and often with 1) solved crystal structures that
#       demonstrate mimicry of host proteins and 2) evolutionary analysis of 
#       where the mimic originated (HGT).
#
# We anticipate creating multiple filtering criteria to isolate different kinds
# of viral structural mimics. As we work with these results, we may also refine
# our ideas of criteria that identify traits of matches (quality, etc.).

library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--input"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--output"), type="character",
              help="Path to output CSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

foldseek_results <- read_tsv(args$input, show_col_types = FALSE) 

# Define constants
TCOV <- 0.25
QCOV <- 0.25
ALNLEN <- 100
EVALUE <- 0.0001
HOST_PDB_PLDDT <- 50
QUERY_PDB_PLDDT <- 50
LDDT <- 0.5
MAX_TMSCORE <- 0.4

foldseek_results_filtered <- foldseek_results %>%
  # filter on quality and strength of foldseek match
  filter(tcov > TCOV) %>%
  filter(qcov > QCOV) %>%
  filter(alnlen > ALNLEN) %>%
  filter(evalue < EVALUE) %>%
  filter(host_pdb_plddt > HOST_PDB_PLDDT) %>%
  filter(query_pdb_plddt > QUERY_PDB_PLDDT) %>%
  filter(lddt > LDDT) %>%
  rowwise() %>%
  mutate(min_tmscore = min(alntmscore, qtmscore, ttmscore),
         max_tmscore = max(alntmscore, qtmscore, ttmscore),
         avg_tmscore = sum(alntmscore, qtmscore, ttmscore) / 3 ) %>%
  ungroup() %>%
  filter(max_tmscore > MAX_TMSCORE) %>%
  dplyr::relocate(all_of(c("avg_tmscore", "min_tmscore")),
                  .after = ttmscore) %>%
  # remove viral matches that hit DNA/RNA replicative machinery
  filter(!grepl(pattern = "polymerase", x = host_function_cc)) %>%
  # filter to host proteins that have signal peptides
  filter(grepl(pattern = "SIGNAL", x=host_signal_peptide)) %>%
  # include extracellular and cell membrane targets
  mutate(host_subcellular_location_cc = tolower(host_subcellular_location_cc)) %>%
  filter(grepl(pattern = "cell membrane|cell junction|cell projection|subcellular location: membrane|secreted",
               x = host_subcellular_location_cc)) %>%
  # change order of columns so it's easier to see the host organism
  dplyr::relocate(host_organism, .before = target) %>%
  # select only the top hit for each query
  arrange(query, evalue) %>%
  group_by(query, query_species) %>%
  slice_min(evalue)

write_csv(foldseek_results_filtered, args$output)
