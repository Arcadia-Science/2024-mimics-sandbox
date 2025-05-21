library(tidyverse)
library(ape)
library(phytools)
library(phylolm)
library(optparse)

# This script uses two approaches to detect "dual mimicry" in viruses. We think
# that at certain points during an infection, viruses need to orchestrate a
# coordinated attack on their host to achieve their goals. Because structural
# mimicry points us toward proteins in the host that can be used to manipulate
# host physiology, we think that investigating instances of "dual mimicry" -
# where a virus mimics two host proteins at the same time to achieve its goals -
# may point us toward more effective strategies for therapeutics based on
# simultaneously targeting two host proteins.
#
# So far, we think that dual mimicry can occur in one of two ways.
# 1. Two viral proteins are co-expressed so that they can work in concert.
# 2. A single viral protein mimics two or more host proteins.
#
# This script analyzes our mimicry data to detect these two types of dual
# mimicry.


# Set up argument parsing ----------------------------------------------------

option_list <- list(
  make_option(c("--input_foldseek_results"), type="character",
              help="Path to foldseek results TSV file."),
  make_option(c("--output_taxonomic_plot"), type="character",
              help="Path to taxonimic tree PNG file."),
  make_option(c("--output_cooccurrence"), type="character",
              help="Path to co-occurrence method output CSV file."),
  make_option(c("--output_dual_domain"), type="character",
              help="Path to dual domain method output CSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

# 1. Phylogenetic co-occurrence analysis -------------------------------------

# We hypothesize that when multiple structural mimics are present in a viral
# genome, they work synergistically to modulate host physiology more effectively
# than single gene alone. We aim to discover these co-occurrence patterns to
# point us toward dual-target drugs that might be more effective than
# single-target drugs. This section uses a phylogenetic approach to control for
# taxonomic relatedness of viruses when testing for significant patterns of
# co-occurrence; viruses that are more closely related should have more genes
# that co-occur becaus of descent from a common ancestor. 

# Read in and filter foldseek results and filter to hits we think are real
foldseek_results <- read_tsv(args$input_foldseek_results) %>%
  group_by(query) %>% 
  slice_min(evalue) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(max_tmscore = max(alntmscore, qtmscore, ttmscore)) %>%
  filter(max_tmscore > 0.5) %>%
  filter(alnlen > 100) %>%
  mutate(host_function_cc = tolower(host_function_cc)) %>%
  # Remove hits to polymerase, which we expect to have high similarity to host
  # polymerase but functions mostly for viral replication.
  filter(!grepl(pattern = "polymerase", x = host_function_cc))

# Extract taxonomy information
taxonomy <- foldseek_results %>%
  select(query_species, query_superkingdom, query_phylum, query_class, 
         query_order, query_family, query_genus) %>%
  distinct()

taxonomy <- taxonomy %>%
  mutate(taxonomy_string = paste(query_superkingdom,
                                 query_phylum,
                                 query_class,
                                 query_order,
                                 query_family,
                                 query_genus,
                                 query_species,
                                 sep = ";"))

taxa_matrix <- taxonomy %>%
  select(query_superkingdom, query_phylum, query_class, query_order,
         query_family, query_genus, query_species) %>%
  column_to_rownames(var = "query_species")

# Build a tree from taxonomy information.
taxa_numeric <- taxa_matrix %>%
  mutate_all(as.factor) %>%
  mutate_all(as.numeric)

dist_matrix <- dist(taxa_numeric, method = "manhattan")

hc <- hclust(dist_matrix, method = "average")

phylo_tree <- as.phylo(hc)

pdf(args$output_taxonomic_plot, width = 4, height = 5)
plot(phylo_tree, cex = 0.6, label.offset = 0.5)
dev.off()

# Pull out gene (protein) presence/absence information
presence_absence <- foldseek_results %>%
  distinct(query_species, host_gene_names_primary) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = host_gene_names_primary, 
              values_from = present, 
              values_fill = list(present = 0)) %>%
  column_to_rownames(var = "query_species")

species_in_data <- rownames(presence_absence)
species_in_tree <- phylo_tree$tip.label

# Check for mismatches, should be empty
setdiff(species_in_data, species_in_tree)

presence_absence <- presence_absence[phylo_tree$tip.label, ]

# Use a phylogenetic comparison method to determine pairs of genes (proteins)
# that co-occur more than we would expect them to by chance, given the
# relatedness of the viruses.

# Get the list of genes
genes <- colnames(presence_absence)
gene_pairs <- combn(genes, 2, simplify = FALSE)
phylo_results <- list()

for (pair in gene_pairs) {
  gene1 <- pair[1]
  gene2 <- pair[2]
  
  # Prepare the data
  data_pair <- data.frame(species = rownames(presence_absence),
                          gene1 = presence_absence[[gene1]],
                          gene2 = presence_absence[[gene2]])
  
  rownames(data_pair) <- data_pair$species
  
  # Make sure the species order matches between data and tree.
  data_pair <- data_pair[phylo_tree$tip.label, ]
  
  # Fit phylogenetic logistic regression model: gene1 ~ gene2
  model <- tryCatch(phyloglm(gene1 ~ gene2,
                             phy = phylo_tree,
                             data = data_pair,
                             method = "logistic_MPLE",
                             btol = 50,
                             log.alpha.bound = 4),
                    error = function(e) NULL)
  
  # If the model didn't converge, skip.
  if (is.null(model)) next
  
  p_value <- summary(model)$coefficients["gene2", "p.value"]
  coefficient <- summary(model)$coefficients["gene2", "Estimate"]
  
  phylo_results[[paste(gene1, gene2, sep = "_")]] <- data.frame(
    gene1 = gene1,
    gene2 = gene2,
    coefficient = coefficient,
    p_value = p_value
  )
}

# Check if the results are empty, and if so write an empty data frame. If not,
# process the results and write out.
if(length(phylo_results) == 0) {
  significant_results = data.frame("gene1" = character(),
                                   "gene2" = character(),
                                   "coefficient" = numeric(),
                                   "p_adjusted" = numeric())
  write_csv(significant_results, args$output_cooccurrence)
} else {
  phylo_results_df <- bind_rows(phylo_results)
  
  phylo_results_df <- phylo_results_df %>%
    mutate(p_adjusted = p.adjust(p_value, method = "BH"))
  
  significant_results <- phylo_results_df %>%
    filter(p_adjusted < 0.05)
  
  significant_results <- significant_results %>% 
    arrange(desc(coefficient)) %>%
    mutate(coefficient = round(coefficient, digits = 1),
           p_adjusted = round(p_adjusted, digits = 3)) %>%
    select(gene1, gene2, coefficient, p_adjusted)
  
  write_csv(significant_results, args$output_cooccurrence)
}

# 2. Detect whether a single viral protein mimics two host proteins ------------

foldseek_results <- read_tsv(args$input_foldseek_results) %>%
  rowwise() %>%
  mutate(max_tmscore = max(alntmscore, qtmscore, ttmscore), .after = ttmscore) %>%
  ungroup() %>%
  filter(evalue <= 1e-5)

# Process each viral query protein to determine if there are strong hits that
# don't overlap in alignment interval with the top hit.
overlaps_with_top_hit <- foldseek_results %>%
  arrange(query, evalue) %>%
  group_by(query) %>%
  # Identify the top hit by lowest e-value and save information about hit.
  mutate(top_target = first(target),
         top_gene_name = first(host_gene_names_primary),
         top_host_function_cc = first(host_function_cc),
         top_max_tmscore = first(max_tmscore),
         top_qstart = first(qstart),
         top_qend = first(qend),
         top_evalue = first(evalue)) %>%
  # Exclude the top hit from comparison.
  filter(row_number() > 1) %>%
  # Calculate overlap with the top hit's alignment interval.
  mutate(overlap_length = pmax(0, pmin(qend, top_qend) - pmax(qstart, top_qstart) + 1),
         alignment_length = qend - qstart + 1,
         top_alignment_length = top_qend - top_qstart + 1,
         overlap_prop = overlap_length / alignment_length,
         overlap_prop_top = overlap_length / top_alignment_length) %>%
  ungroup()

nonoverlapping_hits <- overlaps_with_top_hit %>%
  select(query, target, host_gene_names_primary, max_tmscore, host_function_cc,
         evalue, qstart, qend, alignment_length, top_target, top_gene_name,
         top_host_function_cc, top_max_tmscore, top_qstart, top_qend,
         top_evalue, top_alignment_length, overlap_length, overlap_prop,
         overlap_prop_top) %>%
  # Filter to hits where the alignment intervals overlap by less than 50% of 
  # aligned base pairs in the viral query protein. 50% is empiric and can be
  # lowered to be more stringent; in practice there aren't that many hits, so
  # I allowed a more permissive proportion for the user to make subsequent 
  # decisions about what is interesting. Something like 10% is more
  # conservative. Unlike general structural mimicry, I found no positive
  # control examples of documented mimicry of this type, so I selected these
  # thresholds without using "positive control" examples.
  filter(overlap_prop < 0.5) %>%
  # Select top non-overlapping hit.
  group_by(query) %>%
  slice_min(evalue)

write_tsv(nonoverlapping_hits, args$output_dual_domain)
