library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(DOSE)

#results <- read_csv("~/Downloads/per_cluster_gmm_results030725_with_summary.csv") 
results <- read_csv("~/Downloads/evaluefiltered2.6type2031725.csv")
ggplot(results, aes(x = evalue_max, y = qtmscore_max)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 2.6)


# go enrichment -----------------------------------------------------------

viruses <- results %>% 
  dplyr::filter(evalue_max < 2.6) %>%
  dplyr::select(model_identifier, cluster_members_host_genes) %>%
  separate_rows(cluster_members_host_genes, sep = ", ") %>%
  group_by(model_identifier) %>%
  slice_head(n = 1)

gene_symbols <- viruses$cluster_members_host_genes
gene2entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = c("ENTREZID"))
gene2entrez <- na.omit(gene2entrez)
gene_list <- unique(gene2entrez$ENTREZID)
ego <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL",
                pAdjustMethod= "bonferroni", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                readable = TRUE)

# Now simplify to reduce redundancy
ego_simplified <- simplify(ego, 
                           cutoff = 0.5,    # similarity cutoff
                           by = "p.adjust", # which p-value to use in selecting representative terms
                           select_fun = min # how to select terms among a group
)

head(ego_simplified)

# theme -------------------------------------------------------------------

# font constants
REGULAR_FONT <- "Suisse Int'l"
SEMIBOLD_FONT <- "Suisse Int'l Semi Bold"
MEDIUM_FONT <- "Suisse Int'l Medium"
MONO_FONT <- "Suisse Int'l Mono"

# Use OS default sans fallback if Suisse is not installed
if (!(REGULAR_FONT %in% extrafont::fonts())) {
  FONT <- system("fc-match -f '%{family}' sans", intern = TRUE)
  REGULAR_FONT <- FONT
  SEMIBOLD_FONT <- FONT
  MEDIUM_FONT <- FONT
  MONO_FONT <- FONT
}

# axis label fonts which differ for categorical and numerical data
CATEGORICAL_FONT <- REGULAR_FONT
NUMERICAL_FONT <- MONO_FONT

# main font types
AXIS_TITLE_FONT <- MEDIUM_FONT
KEY_TITLE <- SEMIBOLD_FONT

# background fill options
BACKGROUND_FILL <- "#FDF8F2"

y_axis_type <- "categorical"
x_axis_type <- "numerical"

# font types
x_axis_family <- if (x_axis_type == "categorical") CATEGORICAL_FONT else NUMERICAL_FONT
y_axis_family <- if (y_axis_type == "categorical") CATEGORICAL_FONT else NUMERICAL_FONT
x_axis_label_family <- MEDIUM_FONT
y_axis_label_family <- MEDIUM_FONT
legend_label_family <- SEMIBOLD_FONT
legend_text_family <- if (x_axis_type == "categorical" || y_axis_type == "categorical") CATEGORICAL_FONT else NUMERICAL_FONT # check both because depends on either to change the font

# axis sizes
x_axis_size <- if (x_axis_type == "categorical") 15 else 14.5
y_axis_size <- if (y_axis_type == "categorical") 15 else 14.5

# axis ticks
# size ratio to pts is 2.13, so .75 in pts is .35 here
x_axis_ticks <- if (x_axis_type == "categorical") ggplot2::element_blank() else ggplot2::element_line(color="black", size = 0.35)
# size ratio to pts is 2.13, so .75 in pts is .35 here
y_axis_ticks <- if (y_axis_type == "categorical") ggplot2::element_blank() else ggplot2::element_line(color="black", size = 0.35)

# background fill
background_color <- BACKGROUND_FILL


# plot --------------------------------------------------------------------

dotplot(ego_simplified, showCategory=24, title="GO Enrichment Analysis", label_format = 60) +
  ggplot2::theme(plot.title = ggplot2::element_text(family = REGULAR_FONT, size = 16, face = "bold"),
                 axis.title.x = ggplot2::element_text(family = x_axis_label_family, size = 15, vjust = -1),
                 axis.title.y = ggplot2::element_text(family = y_axis_label_family, size = 15, vjust = +2),
                 axis.text.x = ggplot2::element_text(family = x_axis_family, size = 10, color = "black"),
                 axis.text.y = ggplot2::element_text(family = y_axis_family, size = 10, color = "black"),
                 legend.title = ggplot2::element_text(family = legend_label_family, size = 16),
                 legend.text = ggplot2::element_text(family = legend_text_family, size = 15),
                 # background specifications
                 plot.background = ggplot2::element_rect(fill = background_color, color = NA),
                 panel.background = ggplot2::element_rect(fill = NA, color=NA),
                 panel.border = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 # tick specifications
                 axis.ticks.x = x_axis_ticks,
                 axis.ticks.y = y_axis_ticks,
                 # 5 pixels is roughly 0.07 inches
                 axis.ticks.length.x = ggplot2::unit(0.07, "in"),
                 axis.ticks.length.y = ggplot2::unit(0.07, "in"),
                 # size ratio to pts is 2.13, so .75 in pts is .35 here
                 axis.line = ggplot2::element_line(color="black", size=0.35),
                 
                 # legend specifications
                 legend.background = ggplot2::element_rect(fill = background_color, color = NA))

write_tsv(ego@result, file = "~/Downloads/go_enrichment_results.tsv")
