# The following links/resources were used to create the code in this file:
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/gene-id-annotation_rnaseq_01_ensembl.html
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html
# https://www.r-bloggers.com/2019/05/quick-and-easy-t-sne-analysis-in-r/
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html

# How to clear Rstudio environment: rm(list = ls())
# How to clear Rstudio console: CTRL + L

# Package Installation
install.packages("devtools")
install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Install the Homo sapiens package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(org.Hs.eg.db)

if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

if (!("BiocVersion" %in% installed.packages())) {
  BiocManager::install(version = 'devel')
}

if (!("M3C" %in% installed.packages())) {
  BiocManager::install("M3C")
}

# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(M3C)
library(umap)


# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Create the plots folder if it doesn't exist
if (!dir.exists("plots")) {
  dir.create(plots_dir)
}

# Create the results folder if it doesn't exist
if (!dir.exists("results")) {
  dir.create(results_dir)
}


# Define the file path to the data directory
data_dir <- file.path("data", "SRP018853")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
data_file <- file.path(data_dir, "SRP018853.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
metadata_file <- file.path(data_dir, "metadata_SRP018853.tsv")


file.exists(data_file)
file.exists(metadata_file)



# Part 1.b: Initializing data and converting Ensembl IDs to Hugo gene names

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

library(magrittr)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

library(org.Hs.eg.db)

# Map Ensembl IDs to their associated HUGO gene names
mapped_list <- mapIds(
  org.Hs.eg.db, # Package for humans
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Current ENSEMBL data
  column = "SYMBOL", # Convert to SYMBOL, or HUGO gene names
  multiVals = "list"
)

head(mapped_list)


# Converting list to a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Hugo") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Hugo)

head(mapped_df)



# Part 1.c: Log-scale data, calculate per-gene median expression ranges,
# then make a density plot showing results. 

library(dplyr)

expr_format <- expression_df %>% full_join(mapped_df, by = c("Gene" = "Ensembl")) %>% mutate(Gene = ifelse(!is.na(Hugo), Hugo, Gene))
expr_format = subset(expr_format, select = -Hugo)

head(expr_format)

# What size is your expression matrix? How many genes does it include? 
print(dim(expr_format))
# Expression matrix size: 44293 x 81
# Number of genes: 44293 genes

# Log-scale the data
expr_names = subset(expr_format, select = Gene)
expr_log <- expr_format %>% subset(select = -Gene) %>% mutate(across(everything(), ~log(., base = 2)))
expr_log_full = expr_log
expr_log_full["Gene"] = expr_names
expr_log_full <- expr_log_full %>% relocate(Gene, .before = SRR764979)

# Calculate per-gene median expression ranges
median_ranges2 <- expr_log_full %>%
  rowwise() %>%
  mutate(Median = median(c_across(-Gene))) %>%
  dplyr::select(Gene, Median)

median_ranges <- expr_log %>% rowwise() %>% mutate(Median = median(c_across(), na.rm = TRUE)) %>% subset(select = Median)
median_ranges_full = median_ranges
median_ranges_full["Gene"] = expr_names
median_ranges_full <- median_ranges_full %>% relocate(Gene, .before = Median)

head(median_ranges_full)

# Make a density plot showing those results
densityPlot <- ggplot(data = median_ranges_full, aes(x = Median)) +
  geom_density(fill = "red", color = "black") +
  labs(title = "Density Plot of Median Gene Expression",
       x = "Median Gene Expression at Log2 Scale",
       y = "Density") + ylim(0, .8) + xlim(-.5, 8)

# Save density plot
ggsave(
  filename = file.path("plots", "SRP018853_median_gene_density_plot.png"),
  plot = densityPlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)



# Part 2.a-c: Generate a PCA plot
metadata$TestGroups <- ifelse(stringr::str_starts(metadata$refinebio_title, 'M'), "Healthy", "Pre-T1D")
rounded_df = round(subset(expr_format, select = -Gene))

dds <- DESeqDataSetFromMatrix(countData = rounded_df,
                              colData = metadata,
                              design= ~TestGroups)

# Performance of differential expression analysis
dds <- DESeq(dds)

# Stabilize variance
dds_norm <- vst(dds)

pcaData <- plotPCA(dds_norm, intgroup = "TestGroups", returnData = TRUE)

# metadata$TestGroups <- NULL : To clear and remove TestGroups columns

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=TestGroups)) +
  geom_point(size=3) +
  ggtitle("Principal Component Analysis of Healthy and Pre-T1D Samples")
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
  
  
# Save PCA plot
ggsave(
  filename = file.path("plots", "SRP018853_pca_plot.png"),
  plot = pcaPlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)



# Part 2.d.i: Generate t-SNE plot
set.seed(123)

tsneData <- tsne(rounded_df,labels=as.factor(metadata$TestGroups))

tsnePlot <- last_plot() + 
  ggtitle("t-SNE Plot of Healthy and Pre-T1D Samples")

print(plot)

# Save t-SNE plot
ggsave(
  filename = file.path("plots", "SRP018853_t-sne_plot.png"),
  plot = tsnePlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)



# Part 2.d.ii: Generate UMAP plot
set.seed(246)

normalized_counts <- assay(dds_norm) %>%
  t()

umap_results <- umap::umap(normalized_counts)

umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("refinebio_accession_code") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(metadata, by = "refinebio_accession_code")

umapPlot <- ggplot(umap_plot_df, aes(x = X1, y = X2, color=TestGroups)) +
  geom_point() +
  ggtitle("UMAP Plot of Healthy and Pre-T1D Samples")

# Save UMAP  plot
ggsave(
  filename = file.path("plots", "SRP018853_umap_plot.png"),
  plot = umapPlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)

