# https://www.tidymodels.org/learn/statistics/k-means/index.html 
# https://cran.r-project.org/web/packages/cluster/cluster.pdf 
# https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html


install.packages("readr")
install.packages("ggplot")
install.packages("ggalluvial")
install.packages("cluster")
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("cluster")


library(cluster)
library(readr)
library(ggplot2)
library(tidyr)
library(ggalluvial)
library(readr)
library(cluster)
library(pheatmap)


data <- read_tsv("./data/SRP018853/SRP018853.tsv")
metadata <- read_tsv("./data/SRP018853/metadata_SRP018853.tsv")

variance_result <- apply(data[2, 2:81], 1, var, na.rm = TRUE)
variance_result <- apply(data[2:nrow(data), 2:81], 1, var, na.rm = TRUE)
top_geneNumbers <- order(variance_result, decreasing = TRUE)[1:5000]
top_5000 <- data[top_geneNumbers, ]

data_matrix <- t(data[1:5000, 2:81]) 
data_matrix_10 <- t(data[1:10, 2:81])
data_matrix_100 <- t(data[1:100, 2:81])
data_matrix_1000 <- t(data[1:1000, 2:81])
data_matrix_10000 <- t(data[1:10000, 2:81])
data_matrix_80 <- t(data[1:80, 2:81])



#Unsupervised Analysis using K-means
set.seed(999)

# 2 clusters
kmeans_result_5000_Groups2 <- kmeans(data_matrix, centers = 2)
cluster_assignments <- kmeans_result_5000_Groups2$cluster
pca_5000_Groups2 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups2 <- as.data.frame(pca_5000_Groups2$x)
pca_data_5000_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for k = 2", x = "x1", y = "x2")

# 3 clusters
kmeans_result_5000_Groups3 <- kmeans(data_matrix, centers = 3)
cluster_assignments <- kmeans_result_5000_Groups3$cluster
pca_5000_Groups3 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups3 <- as.data.frame(pca_5000_Groups3$x)
pca_data_5000_Groups3$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups3, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for k = 3", x = "x1", y = "x2")

# 4 clusters
kmeans_result_5000_Groups4 <- kmeans(data_matrix, centers = 4)
cluster_assignments <- kmeans_result_5000_Groups4$cluster
pca_5000_Groups4 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups4 <- as.data.frame(pca_5000_Groups4$x)
pca_data_5000_Groups4$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups4, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for k = 4", x = "x1", y = "x2")

# 5 clusters
kmeans_result_5000_Groups5 <- kmeans(data_matrix, centers = 5)
cluster_assignments <- kmeans_result_5000_Groups5$cluster
pca_5000_Groups5 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups5 <- as.data.frame(pca_5000_Groups5$x)
pca_data_5000_Groups5$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups5, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for k = 5", x = "x1", y = "x2")


# 10 genes
kmeans_result_10_Groups2 <- kmeans(data_matrix_10, centers = 2)
cluster_assignments <- kmeans_result_10_Groups2$cluster
pca_10_Groups2 <- prcomp(data_matrix_10, center = TRUE)
pca_data_10_Groups2 <- as.data.frame(pca_10_Groups2$x)
pca_data_10_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_10_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for n = 10", x = "x1", y = "x2")

# 100 genes
kmeans_result_100_Groups2 <- kmeans(data_matrix_100, centers = 2)
cluster_assignments <- kmeans_result_100_Groups2$cluster
pca_100_Groups2 <- prcomp(data_matrix_100, center = TRUE)
pca_data_100_Groups2 <- as.data.frame(pca_100_Groups2$x)
pca_data_100_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_100_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for n = 100", x = "x1", y = "x2")

# 1,000 genes
kmeans_result_1000_Groups2 <- kmeans(data_matrix_1000, centers = 2)
cluster_assignments <- kmeans_result_1000_Groups2$cluster
pca_1000_Groups2 <- prcomp(data_matrix_1000, center = TRUE, scale. = TRUE)
pca_data_1000_Groups2 <- as.data.frame(pca_1000_Groups2$x)
pca_data_1000_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_1000_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for n = 1000", x = "x1", y = "x2")

# 10,000 genes
kmeans_result_10000_Groups2 <- kmeans(data_matrix_10000, centers = 2)
cluster_assignments <- kmeans_result_10000_Groups2$cluster
pca_10000_Groups2 <- prcomp(data_matrix_10000, center = TRUE, scale. = TRUE)
pca_data_10000_Groups2 <- as.data.frame(pca_10000_Groups2$x)
pca_data_10000_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_10000_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Plot for n = 10000", x = "x1", y = "x2")




#Unsupervised Analysis using PAM clustering
pam_result_10_Groups2 <- pam(data_matrix_10, 2)
plot(pam_result_10_Groups2)

pam_result_10_Groups3 <- pam(data_matrix_10, 3)
plot(pam_result_10_Groups3)

pam_result_10_Groups4 <- pam(data_matrix_10, 4)
plot(pam_result_10_Groups4)

pam_result_10_Groups5 <- pam(data_matrix_10, 5)
plot(pam_result_10_Groups5)

pam_result_80_Groups2 <- pam(data_matrix_80, 2)
plot(pam_result_80_Groups2)




#Unsupervised Analysis using Hierarchical clustering
distance_matrix <- dist(data_matrix)
hclust_result <- hclust(distance_matrix, method = "complete")
plot(hclust_result, labels = FALSE, main = "Hclust dendrogram n = 5000")

# Number of genes will vary, number of clusters is selected for us
distance_matrix_10 <- dist(data_matrix_10)
hclust_result_10 <- hclust(distance_matrix_10, method = "complete")
plot(hclust_result_10, labels = FALSE, main = "Hclust dendrogram n = 10")

distance_matrix_100 <- dist(data_matrix_100)
hclust_result_100 <- hclust(distance_matrix_100, method = "complete")
plot(hclust_result_100, labels = FALSE, main = "Hclust dendrogram n = 100")

distance_matrix_1000 <- dist(data_matrix_1000)
hclust_result_1000 <- hclust(distance_matrix_1000, method = "complete")
plot(hclust_result_1000, labels = FALSE, main = "Hclust dendrogram n = 1000")

distance_matrix_10000 <- dist(data_matrix_10000)
hclust_result_10000 <- hclust(distance_matrix_10000, method = "complete")
plot(hclust_result_10000, labels = FALSE, main = "Hclust dendrogram n = 10000")




#Chi-Squared Test
kmeans_10_100 <- table(kmeans_result_10_Groups2$cluster, kmeans_result_100_Groups2$cluster)
kmeans_10_1000 <- table(kmeans_result_10_Groups2$cluster, kmeans_result_1000_Groups2$cluster)
kmeans_100_1000 <- table(kmeans_result_100_Groups2$cluster, kmeans_result_1000_Groups2$cluster)
kmeans_100_10000 <- table(kmeans_result_100_Groups2$cluster, kmeans_result_10000_Groups2$cluster)

pam_10_80 <- table(pam_result_10_Groups2$clustering, pam_result_80_Groups2$clustering)

hclust_10_100 <- table(hclust_result_10$merge, hclust_result_100$merge)
hclust_10_1000 <- table(hclust_result_10$merge, hclust_result_1000$merge)
hclust_100_1000 <- table(hclust_result_100$merge, hclust_result_1000$merge)
hclust_100_10000 <- table(hclust_result_100$merge, hclust_result_10000$merge)


# K-means
chisq_kmeans_10_100 <- chisq.test(kmeans_10_100)
chisq_kmeans_10_1000 <- chisq.test(kmeans_10_1000)
chisq_kmeans_100_1000 <- chisq.test(kmeans_100_1000)
chisq_kmeans_100_10000 <- chisq.test(kmeans_100_10000)

# PAM
chisq_pam_10_80 <- chisq.test(pam_10_80)

# Hierarchical clustering
chisq_hclust_10_100 <- chisq.test(hclust_10_100)
chisq_hclust_10_1000 <- chisq.test(hclust_10_1000)
chisq_hclust_100_1000 <- chisq.test(hclust_100_1000)
chisq_hclust_100_10000 <- chisq.test(hclust_100_10000)

# Create a data frame to store the results
gene_clust_chisq <- data.frame(
  Method = c("K-means", "K-means", "K-means", "K-means", "PAM", "Hierarchical", "Hierarchical", "Hierarchical", "Hierarchical"),
  Gene_Count_1 = c(10, 10, 100, 100, 10, 10, 10, 100, 100),
  Gene_Count_2 = c(100, 1000, 1000, 10000, 80, 100, 1000, 1000, 10000),
  Chi_Squared_p_value = c(chisq_kmeans_10_100$p.value,
                          chisq_kmeans_10_1000$p.value,
                          chisq_kmeans_100_1000$p.value,
                          chisq_kmeans_100_10000$p.value,
                          chisq_pam_10_80$p.value,
                          chisq_hclust_10_100$p.value,
                          chisq_hclust_10_1000$p.value,
                          chisq_hclust_100_1000$p.value,
                          chisq_hclust_100_10000$p.value)
)

write_tsv(gene_clust_chisq, file.path("results", "gene_clustering_chisquare.tsv"))




#Alluvial diagram

# Create a data frame of cluster assignments
cluster_data <- data.frame(
  Sample = rownames(data_matrix),
  Genes_10 = kmeans_result_10_Groups2$cluster,
  Genes_100 = kmeans_result_100_Groups2$cluster,
  Genes_1000 = kmeans_result_1000_Groups2$cluster,
  Genes_5000 = kmeans_result_5000_Groups2$cluster,
  Genes_10000 = kmeans_result_10000_Groups2$cluster
)

long_data_cluster <- pivot_longer(cluster_data, 
                                  cols = starts_with("Genes"), 
                                  names_to = "Gene_Count", 
                                  values_to = "cluster")

long_data_cluster$cluster <- as.factor(long_data_cluster$cluster)

alluvial_plot <- ggplot(
  long_data_cluster,
  aes(x = Gene_Count, stratum = cluster, alluvium = Sample, fill = cluster, label = cluster)) +
  scale_fill_brewer(type = "qual", palette = "Pastel1") + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Alluvial Diagram of Cluster Memberships Across Gene Subsets")

# Save the plot
ggsave(
  file.path("Assn3 Plots", "alluvial_diagram.png"),
  plot = alluvial_plot,
  device = "png",
  width = 8,
  height = 6,
  units = "in"
)




# Heatmap
metadata$TestGroups <- ifelse(stringr::str_starts(metadata$refinebio_title, 'M'), "Healthy", "Pre-T1D")

annotations <- data.frame(
  kmeans = kmeans_result_5000_Groups2$cluster,
  pam = pam_result_80_Groups2$clustering,
  hierarchical = hclust_result$order,  
  Assn1_TestGroups = metadata$TestGroups
)

data_matrix <- as.matrix(data.matrix)

pheatmap(
  data_matrix,
  annotation_col = annotations,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  scale = "none",
  main = "Heatmap with Clustering Annotations"
)




# Statistics - Chi-Squared Test
# Assignment 1 vs. Clusters


