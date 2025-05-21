library(tidyverse)

thresh5 <- read_csv("~/Downloads/per_cluster_gmm_results022825_typ2thres05_022825.csv")  
thresh15 <- read_csv("~/Downloads/per_cluster_gmm_results022825_typ2thres15_fixed.csv")

nrow(thresh5)
nrow(thresh15)

table(thresh5$cluster_member_queries %in% thresh15$cluster_member_queries)
table(thresh15$cluster_member_queries %in% thresh5$cluster_member_queries)
# 539 rows are the same between threshold 5 and 15 by query members

table(thresh5$cluster_members_host_genes %in% thresh15$cluster_members_host_genes)
table(thresh15$cluster_members_host_genes %in% thresh5$cluster_members_host_genes)
# 279 or 281 sets of host genes are the samebetween threshold 5 and 15 by host gene

thresh5 <- thresh5 %>%
  filter(evalue_min < 1e-5) %>%
  filter(evalue_max < 0.01)
write_tsv(thresh5, "~/Downloads/tmp.tsv")
thresh15 <- thresh15 %>%
  filter(evalue_min < 1e-5) %>%
  filter(evalue_max < 0.01)

table(thresh5$cluster_member_queries %in% thresh15$cluster_member_queries)
table(thresh15$cluster_member_queries %in% thresh5$cluster_member_queries)
# 117 rows are the same between threshold 5 and 15 by query members after filtering by evalue

table(thresh5$cluster_members_host_genes %in% thresh15$cluster_members_host_genes)
table(thresh15$cluster_members_host_genes %in% thresh5$cluster_members_host_genes)
# 102 or 101 sets of host genes are the samebetween threshold 5 and 15 by host gene



# split out individual identifiers to see 1:1 matches ---------------------
thresh5_sep <- thresh5 %>%
  separate_rows(cluster_members_host_genes, sep = ",\\s*") %>%
  separate_rows(cluster_member_queries, sep = ",\\s*") %>%
  unite(col = query_target, cluster_member_queries, cluster_members_host_genes, sep = "_", remove = FALSE)

length(unique(thresh5_sep$cluster_member_queries))
length(unique(thresh5_sep$cluster_members_host_genes))
# thresh5 has 449 distinct virus proteins that match 460 distinct human genes


thresh15_sep <- thresh15 %>%
  separate_rows(cluster_members_host_genes, sep = ",\\s*") %>%
  separate_rows(cluster_member_queries, sep = ",\\s*") %>%
  unite(col = query_target, cluster_member_queries, cluster_members_host_genes, sep = "_", remove = FALSE)

length(unique(thresh15_sep$cluster_member_queries))
length(unique(thresh15_sep$cluster_members_host_genes))
# thresh15 has 491 distinct virus proteins that match 486 distinct human genes

# which unique pairs are recovered in one and not the other?
thresh5_unique <- thresh5_sep %>%
  filter(!query_target %in% thresh15_sep$query_target)
nrow(thresh5_unique)
# 223 unique pairs


thresh15_unique <- thresh15_sep %>%
  filter(!query_target %in% thresh5_sep$query_target)
nrow(thresh15_unique)
# 499 unique


#...but there is some overlap in the host proteins, just not 1:1 hits
table(thresh15_unique$cluster_members_host_genes %in% thresh5_unique$cluster_members_host_genes)
table(unique(thresh15_unique$cluster_members_host_genes) %in% unique(thresh5_unique$cluster_members_host_genes))

table(thresh5_unique$cluster_members_host_genes %in% thresh15_unique$cluster_members_host_genes)
table(unique(thresh5_unique$cluster_members_host_genes) %in% unique(thresh15_unique$cluster_members_host_genes))
