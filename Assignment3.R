
install.packages("readr")
library(readr) 
data <- read_tsv("./data/SRP018853/SRP018853.tsv")
metadata <- read_tsv("./data/SRP018853/metadata_SRP018853.tsv")
variance_result <- apply(data[2, 2:81], 1, var, na.rm = TRUE)
variance_result <- apply(data[2:nrow(data), 2:81], 1, var, na.rm = TRUE)
top_geneNumbers <- order(variance_result, decreasing = TRUE)[1:5000]
top_5000 <- data[top_geneNumbers, ]
data_matrix <- data[1:5000, 2:81]
data_matrix_10 <- data[1:10, 2:81]
data_matrix_100 <- data[1:100, 2:81]
data_matrix_1000 <- data[1:1000, 2:81]
data_matrix_10000 <- data[1:10000, 2:81]

#K-means
set.seed(999)
kmeans_result <- kmeans(data_matrix, centers = 2)
cluster_assignments <- kmeans_result$cluster
pca <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca$x)
pca_data$cluster <- as.factor(cluster_assignments)
install.packages("ggplot")
library(ggplot2)
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

#PAM
install.packages("cluster")
library(cluster)
pam_result <- pam(data_matrix, 2)
plot(pam_result)

#Hclust
distance_matrix <- dist(data_matrix)
hclust_result <- hclust(distance_matrix, method = "complete")
plot(hclust_result, labels = FALSE, main = "Hclust dendrogram n = 5000")
