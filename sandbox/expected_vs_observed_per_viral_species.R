library(tidyverse)
library(janitor)
setwd("~/github/2024-parasite-human-mimics/")

viral_metadata <- read_tsv("inputs/viral/merged_viral_metadata_human.tsv") %>%
  clean_names() %>%
  mutate(structure_file = str_remove(pattern = ".pdb$", string = structure_file))

viral_count <- viral_metadata %>%
  group_by(family) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(family = reorder(family, -count))

# Create the plot
ggplot(viral_count, aes(x = family, y = count)) +
  geom_col() +
  theme_classic() +
  coord_flip()


details <- read_csv("~/Downloads/best_cluster_rows_030725_detailed.csv")

viral_metadata <- viral_metadata %>%
  mutate(in_best_hit = ifelse(tolower(structure_file) %in% details$query, "best_hit", "not_best_hit"))

viral_count2 <- viral_metadata %>%
  group_by(family, in_best_hit) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(family = reorder(family, -count))
viral_count2$in_best_hit <- factor(viral_count2$in_best_hit, levels = c("not_best_hit", "best_hit"))

# Create the plot
ggplot(viral_count2, aes(x = family, y = count, fill = in_best_hit)) +
  geom_col() +
  theme_classic() +
  coord_flip()


# chi square --------------------------------------------------------------


# Total counts of 'best_hit' and 'not_best_hit' for each family
counts <- table(viral_metadata$family, viral_metadata$in_best_hit)

# Calculate observed proportions of 'best_hit' in each family
observed_proportions <- prop.table(counts, 1)[, "best_hit"]

# Overall 'best_hit' proportion
overall_proportion <- sum(viral_metadata$in_best_hit == "best_hit") / nrow(viral_metadata)

# Assuming genes are distributed equally among families (adjust if not)
total_genes_per_family <- table(viral_metadata$family)

# Expected counts of 'best_hit' per family under randomness
expected_counts <- total_genes_per_family * overall_proportion

observed_counts <- table(viral_metadata$family, viral_metadata$in_best_hit)[, "best_hit"]
observed_counts[is.na(observed_counts)] <- 0

# Perform the test
chi_square_results <- chisq.test(x = observed_counts, p = expected_counts, rescale.p = TRUE)

# Print the results
print(chi_square_results)


# plot --------------------------------------------------------------------


tmp <- data.frame(observed_counts, expected_counts) %>%
  select(family = Var1, observed = observed_counts, expected = Freq) 

# Pivot data to long format using tidyr
data_long <- tmp %>%
  pivot_longer(cols = -family, names_to = "type", values_to = "count")

# Plotting
ggplot(data_long, aes(x = family, y = count, fill = type)) +
  geom_col(position = position_dodge(width = .5), width = .6) +
  theme_classic() +
  coord_flip() 
