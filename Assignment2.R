library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(ggplot2)
if (!("org.Dr.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(org.Hs.eg.db)

# Define the file path to the data directory
# Replace with the path of the folder the files will be in
# !! Update this for wherever you place the repo on your computer !!
data_dir <- file.path("D:/_schoolwork/16th Grade/CGS4144 bioinformatics/repo/CGS4144-Group-Project/Dataset", "SRP018853")
img_dir <- file.path("D:/_schoolwork/16th Grade/CGS4144 bioinformatics/repo/CGS4144-Group-Project", "Graphs")
data_file <- file.path(data_dir, "SRP018853.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP018853.tsv")

metadata <- readr::read_tsv(metadata_file)
dataset <- readr::read_tsv(data_file)

expression_df <- readr::read_tsv(data_file) %>% tibble::column_to_rownames("Gene")
expression_df <- expression_df %>% dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

expression_df <- expression_df %>% tibble::rownames_to_column("Gene")

mapped_list <- mapIds(
  org.Hs.eg.db, # Package for humans
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # We have ENSEMBL
  column = "SYMBOL", # Want SYMBOL
  multiVals = "list"
)
head(mapped_list)

# enframe() makes a `list` column; we will simplify it with unnest()
# This will result in one row of our data frame per list item
mapped_df <- mapped_list %>% tibble::enframe(name = "Ensembl", value = "Symbol") %>% tidyr::unnest(cols = Symbol)

head(mapped_df)
summary(as.factor(mapped_df$Symbol), maxsum = 100)

multi_mapped <- mapped_df %>% dplyr::count(Ensembl, name = "symbol_id_count") %>% dplyr::arrange(desc(symbol_id_count))

head(multi_mapped)


# c. Load the data into your chosen programming language (R or python recommended). 
# What size is your expression matrix? How many genes does it include? 
# How much variation do you see in the data? 

# To answer these questions, log-scale the data, calculate per-gene median expression ranges, 
# then make a density plot showing those results. Summarize your findings.

expr_format <- expression_df %>% full_join(mapped_df, by = c("Gene" = "Ensembl")) %>% mutate(Gene = ifelse(!is.na(Symbol), Symbol, Gene))
expr_format = subset(expr_format, select = -Symbol)

head(expr_format)
# What size is your expression matrix? How many genes does it include? 
print(dim(expr_format))
# 44181 x 81, 44181 genes

# log-scale the data
expr_names = subset(expr_format, select = Gene)
expr_log <- expr_format %>% subset(select = -Gene) %>% mutate(across(everything(), ~log(., base = 2)))
expr_log_full = expr_log
expr_log_full["Gene"] = expr_names
expr_log_full <- expr_log_full %>% relocate(Gene, .before = SRR764979)

# calculate per-gene median expression ranges
median_ranges2 <- expr_log_full %>%
  rowwise() %>%
  mutate(Median = median(c_across(-Gene))) %>%
  dplyr::select(Gene, Median)

median_ranges <- expr_log %>% rowwise() %>% mutate(Median = median(c_across(), na.rm = TRUE)) %>% subset(select = Median)
median_ranges_full = median_ranges
median_ranges_full["Gene"] = expr_names
median_ranges_full <- median_ranges_full %>% relocate(Gene, .before = Median)

head(median_ranges_full)

# make a density plot showing those results
ggplot(data = median_ranges_full, aes(x = Median)) +
  geom_density(fill = "red", color = "black") +
  labs(title = "Density Plot of Median Gene Expression",
       x = "Median Gene Expression at Log2 Scale",
       y = "Density") + ylim(0, .8) + xlim(-.5, 8)
# and save it 
ggsave(
  filename = file.path(img_dir, "SRP018853_median_gene_density_plot.png"),
  plot = last_plot(),
  device = "png",
)

#section 2