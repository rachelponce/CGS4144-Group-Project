# The following links/resources were used to create the code in this file:
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/gene-id-annotation_rnaseq_01_ensembl.html
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html
# https://www.r-bloggers.com/2019/05/quick-and-easy-t-sne-analysis-in-r/
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html
# https://alexslemonade.github.io/refinebio-examples/02-microarray/clustering_microarray_01_heatmap.html#4_Clustering_Heatmap_-_Microarray
# https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html 
# https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/




# How to clear Rstudio environment
rm(list = ls())

# How to clear Rstudio console: CTRL + L
# Delete a data frame: rm(df1)


# Package Installation
install.packages("devtools")
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("gprofiler2")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("apeglm")
BiocManager::install("M3C")
BiocManager::install("ComplexHeatmap")
BiocManager::install("umap")
BiocManager::install("clusterProfiler")
BiocManager::install("msigdbr")
BiocManager::install("AnnotationDbi")
BiocManager::install("topGO")
BiocManager::install("biomaRt")
BiocManager::install("Rgraphviz")


# Install the Homo sapiens package
BiocManager::install("org.Hs.eg.db") 

if (!("BiocVersion" %in% installed.packages())) {
  BiocManager::install(version = 'devel')
}


# Libraries
library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(M3C)
library(umap)
library(ComplexHeatmap)
library(devtools)
library(circlize)
library(EnhancedVolcano)
library(pheatmap)
library(gprofiler2)
library(topGO)
library(biomaRt)
library(Rgraphviz)
library(clusterProfiler)
library(msigdbr)
library(AnnotationDbi)



# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Create the plots folder if it doesn't exist
if (!dir.exists("Assn2-plots")) {
  dir.create("Assn2-plots")
}

# Create the results folder if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
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
metadata <- read_tsv(metadata_file)

# Read in data TSV file
expression_df <- read_tsv(data_file) %>%
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

summary(as.factor(mapped_df$Hugo), maxsum = 10)
# 13,721 NAs: # of Ensembl IDs that did not map to HUGO gene names 

collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
  # Collapse the HUGO gene names `mapped_df` into one column named `all_HUGO_names`
  dplyr::summarize(all_HUGO_names = paste(Hugo, collapse = ";"))

head(collapsed_mapped_df)

final_mapped_df <- data.frame(
  "first_HUGO_name" = mapIds(
    org.Hs.eg.db,
    keys = expression_df$Gene,
    keytype = "ENSEMBL", 
    column = "SYMBOL",
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))

final_mapped_df %>%
  # Filter `final_mapped_df` to rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_HUGO_names, ";")) %>%
  head()

write_tsv(final_mapped_df, file.path("results","Ensembl_HUGO_names.tsv"))




# Part 1.c: Log-scale data, calculate per-gene median expression ranges,
# then make a density plot showing results. 

expr_format <- expression_df %>% full_join(mapped_df, by = c("Gene" = "Ensembl")) %>% mutate(Gene = ifelse(!is.na(Hugo), Hugo, Gene))
expr_format = subset(expr_format, select = -Hugo)

head(expr_format)

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
  filename = file.path("Assn2-plots", "median_gene_density_plot.png"),
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
ddsAnalysis <- DESeq(dds)

# Stabilize variance
dds_norm <- vst(ddsAnalysis)

pcaData <- plotPCA(dds_norm, intgroup = "TestGroups", returnData = TRUE)

# metadata$TestGroups <- NULL : To clear and remove TestGroups columns

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=TestGroups)) +
  geom_point(size=3) +
  ggtitle("Principal Component Analysis of Healthy and Pre-T1D Samples") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
  
  
# Save PCA plot
ggsave(
  filename = file.path("Assn2-plots", "pca_plot.png"),
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

# Save t-SNE plot
ggsave(
  filename = file.path("Assn2-plots", "t-sne_plot.png"),
  plot = tsnePlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)




# Part 2.d.ii: Generate UMAP plot
set.seed(246)

expr_df2 <- read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

expr_df2 <- expr_df2 %>%
  dplyr::select(metadata$refinebio_accession_code)

all.equal(colnames(expr_df2), metadata$refinebio_accession_code)

metadata2 <- metadata %>%
  dplyr::select( # Select only the columns needed for plotting
    refinebio_accession_code,
    TestGroups
  ) %>%
  dplyr::mutate( # Convert the annotation variables into factors
    TestGroups = factor(
      TestGroups,
      levels = c("Healthy", "Pre-T1D")
    )
  )

filtered_expression_df <- expr_df2 %>%
  dplyr::filter(rowSums(.) >= 200)

filtered_expression_df <- round(filtered_expression_df)

dds2 <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df, # Counts values for all samples in dataset
  colData = metadata2, # annotation data for the samples in the counts data frame
  design = ~1
)

dds2_norm <- vst(dds2)

normalized_counts <- assay(dds2_norm) %>%
  t()

umap_results <- umap::umap(normalized_counts)

umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("refinebio_accession_code") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(metadata2, by = "refinebio_accession_code")

umap_plot_df

umapPlot <- ggplot(umap_plot_df, aes(x = X1, y = X2, color=TestGroups)) +
  geom_point() +
  ggtitle("UMAP Plot of Healthy and Pre-T1D Samples")

umapPlot

# Save UMAP  plot
ggsave(
  filename = file.path("Assn2-plots", "umap_plot.png"),
  plot = umapPlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)




# Part 3.a Perform differential analysis
set.seed(369)

exprdf_DA <- expression_df %>%
  tibble::column_to_rownames("Gene")

exprdf_DA <- exprdf_DA %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(exprdf_DA), metadata$refinebio_accession_code)

metadata <- metadata %>%
  dplyr::mutate(
    TestGroups = factor(TestGroups, levels = c("Healthy", "Pre-T1D"))
  )

filter_exprdf_DA <- exprdf_DA %>%
  dplyr::filter(rowSums(.) >= 75)

gene_matrix <- round(filter_exprdf_DA)

ddset <- DESeqDataSetFromMatrix(
  # Non-normalized count data
  countData = gene_matrix,
  # Metadata data frame
  colData = metadata,
  # Supply experimental variable
  design = ~TestGroups
)

deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(
  deseq_object, # Original DESeq2 object after running DESeq()
  coef = 2, # Log fold change coefficient used in DESeq(); default is 2.
  res = deseq_results # Original DESeq2 results table
)

head(deseq_results)

# Convert data set to data frame
deseq_df <- deseq_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.05) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)

write_tsv(deseq_df,file.path("results","diff_expr_results.tsv"))




# Part 3.b: Create a volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01,
)

#Save volcano plot
ggsave(
  file.path("Assn2-plots", "volcano_plot.png"),
  plot = volcano_plot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)




# Part 3.c: Create a table of the top 50 differentially expressed genes
top_50_genes <- deseq_df %>%
  dplyr::arrange(padj, desc(log2FoldChange)) %>%  # Sort by padj and log2FoldChange
  dplyr::slice_head(n = 50)


# Save the top 50 genes to a TSV file
write_tsv(top_50_genes, file.path("results", "top_50_diff_expr_genes.tsv"))




# Part 4: Create a heatmap
sigGenes <- read_tsv("results/top_50_diff_expr_genes.tsv")

df <- read_tsv(data_file) 

top50_geneNames <- sigGenes$Gene

df <- df %>%
  dplyr::filter(Gene %in% top50_geneNames)

df <- data.frame(df) %>% 
  tibble::column_to_rownames("Gene")

annotation_df <- metadata %>%
  # Select the variables that we want for annotating the heatmap
  dplyr::select(
    refinebio_accession_code,
    TestGroups
  ) %>%
  # The `pheatmap()` function requires that the row names of our
  # annotation object matches the column names of our dataset object
  tibble::column_to_rownames("refinebio_accession_code")

heatmap <- pheatmap(
  df,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples)
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_df,
  main = "Annotated Heatmap",
  colorRampPalette(c(
    "green",
    "white",
    "darkmagenta"
  ))(25),
  scale = "row" # Scale values in the direction of genes (rows)
)

ggsave(
  file.path("Assn2-plots", "heatmap.png"),
  plot = heatmap,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)




# Part 5.1 Enrichment analysis using gProfiler2 and gene ontology
diffGeneList <- read_tsv("results/diff_expr_results.tsv")
diffGeneList <- diffGeneList$Gene

gostres <- gost(query = diffGeneList, 
                organism = "hsapiens", 
                ordered_query = FALSE, 
                multi_query = FALSE, 
                significant = TRUE, 
                exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, 
                evcodes = FALSE, 
                user_threshold = 0.05, 
                correction_method = "g_SCS", 
                domain_scope = "annotated", 
                custom_bg = NULL, 
                numeric_ns = "", 
                sources = c("GO:BP"), 
                as_short_link = FALSE, 
                highlight = TRUE)

names(gostres)

gProf_table <- gostres$result

gostplot(gostres, capped = TRUE, interactive = TRUE)
gProfPlot <- gostplot(gostres, capped = FALSE, interactive = FALSE)

ggsave(
  file.path("Assn2-plots", "gProfiler_plot.png"),
  plot = gProfPlot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)

write_tsv(gProf_table, file.path("results", "gProfiler_EA.tsv"))




# Part 5.2 Enrichment analysis using clustProfiler and gene ontology
sigs <- read_tsv("results/diff_expr_results.tsv")
clustGeneList <- sigs$log2FoldChange
names(clustGeneList) <- sigs$Gene

clustProf_results <- gseGO(
  clustGeneList,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  keyType = "ENSEMBL",
  eps = 1e-300
)

clustProf_results_df <- as.data.frame(clustProf_results)
gseaplot(clustProf_results, geneSetID = 1)

write_tsv(clustProf_results_df, file.path("results","clustProfiler_EA.tsv"))




# Part 5.3 Enrichment analysis using topGO and gene ontology
totalGenes <- read_tsv("results/Ensembl_HUGO_names.tsv")
totalGeneList <- na.omit(totalGenes$first_HUGO_name)

length(totalGeneList)
head(totalGeneList)

mapped_diffGenes <- mapIds(
  org.Hs.eg.db, 
  keys = diffGeneList,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "list"
)

head(mapped_diffGenes)

mapped_df_diffGenes <- mapped_diffGenes %>%
  tibble::enframe(name = "Ensembl", value = "Hugo") %>%
  tidyr::unnest(cols = Hugo)

diffGeneList_HUGO <- na.omit(mapped_df_diffGenes$Hugo)

length(diffGeneList_HUGO)
head(diffGeneList_HUGO)

db = useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids = getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=totalGeneList, mart=db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = diffGeneList_HUGO %in% go_ids[,2]
keep =which(keep==TRUE)
diffGeneList_HUGO=diffGeneList_HUGO[keep]

# make named factor showing which genes are of interest
geneList=factor(as.integer(totalGeneList %in% diffGeneList_HUGO))
names(geneList)= totalGeneList

topGO_data <- new("topGOdata", 
                  ontology = "BP",
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = gene_2_GO)

classic_fisher_result = runTest(topGO_data, algorithm='classic', statistic='fisher')

weight_fisher_result = runTest(topGO_data, algorithm='weight01', statistic='fisher') 

allGO=usedGO(topGO_data)
all_res <- GenTable(topGO_data, classicFisher= classic_fisher_result, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_final=cbind(all_res,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]

#get list of significant GO after multiple testing correction
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]

write_tsv(all_res_final,file.path("results","topGO_EA.tsv"))

# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
pdf(file='Assn2-plots/topGO_plot.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(topGO_data, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
dev.off()

myterms = results.table.bh$GO.ID
mygenes = genesInTerm(topGO_data, myterms)

var=c()
for (i in 1:length(myterms))
{
  myterm=myterms[i]
  mygenesforterm= mygenes[myterm][[1]]
  mygenesforterm=paste(mygenesforterm, collapse=',')
  var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
write.table(var,"genetoGOmapping.txt",sep="\t",quote=F)

var <- as.data.frame(var)

write_tsv(var,file.path("results","topGO_genetoGOmapping.tsv"))




# Part 6: Combined table of enrichment analysis results
gprof_results <- read_tsv("results/gProfiler_EA.tsv")
clust_results <- read_tsv("results/clustProfiler_EA.tsv")
topgo_results <- read_tsv("results/topGO_EA.tsv")

gprof_results <- gprof_results %>%
  select(
    GO.ID = term_id,
    term = term_name,
    gProfiler.p_value = p_value,
    significant_gprof = significant
  )

clust_results <- clust_results %>%
  select(
    GO.ID = ID,
    term = Description,
    clustProfiler.p_value = pvalue,
  )%>%
  mutate(
    significant_clustprof = TRUE
  )

topgo_results <- topgo_results %>%
  select(
    GO.ID = GO.ID,
    term = Term,
    topGO.p_value = weightFisher,
    significant_topgo = Significant
  )

topgo_results <- topgo_results %>%
  mutate(
    significant_topgo = ifelse(significant_topgo != 0, "TRUE", "FALSE")
  )

# Merge results by Gene Set / Term
merged_results <- gprof_results %>%
  full_join(clust_results, by = c("GO.ID", "term")) %>%
  full_join(topgo_results, by = c("GO.ID", "term"))

merged_results <- merged_results %>%
  mutate(
    across(starts_with("significant"), ~ as.logical(.)),
    across(ends_with("p_value"), ~ as.numeric(.))
  )

sorted_results <- merged_results %>%
  arrange(gProfiler.p_value, clustProfiler.p_value, topGO.p_value)

enriched_terms <- sorted_results %>%
  rowwise() %>%
  mutate(
    significant_count = sum(c_across(starts_with("significant")) == "TRUE", na.rm = TRUE),
    included_count = sum(!is.na(c_across(ends_with("p_value"))))
  ) %>%
  ungroup()

enriched_terms <- enriched_terms %>% 
  select(-c(significant_gprof, significant_clustprof, significant_topgo))

write_tsv(enriched_terms, "results/merged_enrichment_analysis.tsv")

top10_enriched_terms <- enriched_terms %>% 
  filter(significant_count >= 2) %>%
  filter(included_count == 3) %>%
  slice_head(n = 10)

write_tsv(top10_enriched_terms, "results/top10_terms_EA.tsv")
