library(tidyverse)
library(vegan)
setwd("~/github/2024-parasite-human-mimics/")
clusters <- read_tsv("sandbox/viro3d_clusters/viro3d_human_viruses_tmscorethreshold0.7_cluster.tsv",
                     col_names = c("representative", "member")) %>%
  mutate(representative = str_replace(string = representative, pattern = ".*?(CF-|EF-)", replacement = "\\1")) %>%
  mutate(member = str_replace(string = member, pattern = ".*?(CF-|EF-)", replacement = "\\1")) %>%
  distinct()
  
# remove viral proteins that are singlets (only members in their clusters)
singlet_clusters <- clusters %>%
  group_by(representative) %>%
  tally() %>%
  filter(n == 1) %>%
  select(-n)

type1_results_full <- read_tsv("~/Downloads/viro3d_virus_matches_alignmenttype1.tsv.gz", show_col_types = FALSE) %>%
  select(query_virus_name, query, target, qlen, tlen, alnlen, alntmscore, qtmscore, 
         ttmscore, host_gene_names_primary, host_function_cc, lddt, prob, qcov,
         tcov, pident, bits, evalue, qstart, qend, tstart, tend, 
         query_structure_file, query_protlen, query_genome_composition,
         query_kingdom, query_phylum, query_class, query_order, query_family,
         query_genus, query_species, query_virus_names) %>%
  # remove hits that are in singlet clusters
  mutate(query = str_replace(string = query, pattern = ".*?(CF-|EF-)", replacement = "\\1"),
         query = str_remove(string = query, pattern = ".pdb$")) %>% 
  filter(!query %in% singlet_clusters$representative)

type1_results <- type1_results_full %>%
  mutate(qtmscore = ifelse(qtmscore >= 0.5, 1, 0)) %>%
  select(query, target, qtmscore) %>%
  distinct() %>%
  pivot_wider(id_cols = query, names_from = target, values_from = qtmscore, values_fill = 0) %>%
  column_to_rownames("query")

metadata <- 
  rownames_to_column(as.data.frame(type1_results), var = "query") %>%
  left_join(clusters, by = c("query" = "member")) %>%
  mutate(cluster = as.factor(representative)) %>%
  select(query, cluster)

# make sure metadata is ordered the same as the dissimilarity object
metadata <- column_to_rownames(metadata, "query")
all.equal(rownames(metadata), rownames(type1_results))

type1_dissimilarity <- vegan::vegdist(type1_results, method="jaccard")

adonis_res <- adonis2(type1_dissimilarity ~ cluster, data = metadata, 
                      parallel = 6)
print(adonis_res)

adonis_res9 <- adonis_res




adonis_res85 <- adonis_res
adonis_res8 <- adonis_res
adonis_res75 <- adonis_res
adonis_res7 <- adonis_res

adonis_res6 <- adonis_res
adonis_res65 <- adonis_res
