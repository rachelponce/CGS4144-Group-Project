
install.packages("readr")
library(readr)
install.packages("ggplot")
library(ggplot2)

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
kmeans_result_5000_Groups2 <- kmeans(data_matrix, centers = 2)
cluster_assignments <- kmeans_result_5000_Groups2$cluster
pca_5000_Groups2 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups2 <- as.data.frame(pca_5000_Groups2$x)
pca_data_5000_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_5000_Groups3 <- kmeans(data_matrix, centers = 3)
cluster_assignments <- kmeans_result_5000_Groups3$cluster
pca_5000_Groups3 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups3 <- as.data.frame(pca_5000_Groups3$x)
pca_data_5000_Groups3$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups3, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_5000_Groups4 <- kmeans(data_matrix, centers = 4)
cluster_assignments <- kmeans_result_5000_Groups4$cluster
pca_5000_Groups4 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups4 <- as.data.frame(pca_5000_Groups4$x)
pca_data_5000_Groups4$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups4, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_5000_Groups5 <- kmeans(data_matrix, centers = 5)
cluster_assignments <- kmeans_result_5000_Groups5$cluster
pca_5000_Groups5 <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
pca_data_5000_Groups5 <- as.data.frame(pca_5000_Groups5$x)
pca_data_5000_Groups5$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_5000_Groups5, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_10_Groups2 <- kmeans(data_matrix_10, centers = 2)
cluster_assignments <- kmeans_result_10_Groups2$cluster
pca_10_Groups2 <- prcomp(data_matrix_10, center = TRUE)
pca_data_10_Groups2 <- as.data.frame(pca_10_Groups2$x)
pca_data_10_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_10_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_100_Groups2 <- kmeans(data_matrix_100, centers = 2)
cluster_assignments <- kmeans_result_100_Groups2$cluster
pca_100_Groups2 <- prcomp(data_matrix_100, center = TRUE)
pca_data_100_Groups2 <- as.data.frame(pca_100_Groups2$x)
pca_data_100_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_100_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_1000_Groups2 <- kmeans(data_matrix_1000, centers = 2)
cluster_assignments <- kmeans_result_1000_Groups2$cluster
pca_1000_Groups2 <- prcomp(data_matrix_1000, center = TRUE, scale. = TRUE)
pca_data_1000_Groups2 <- as.data.frame(pca_1000_Groups2$x)
pca_data_1000_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_1000_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")

kmeans_result_10000_Groups2 <- kmeans(data_matrix_10000, centers = 2)
cluster_assignments <- kmeans_result_10000_Groups2$cluster
pca_10000_Groups2 <- prcomp(data_matrix_10000, center = TRUE, scale. = TRUE)
pca_data_10000_Groups2 <- as.data.frame(pca_10000_Groups2$x)
pca_data_10000_Groups2$cluster <- as.factor(cluster_assignments)
ggplot(pca_data_10000_Groups2, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    labs(title = "K-means plot", x = "x1", y = "x2")
#PAM
install.packages("cluster")
library(cluster)
pam_result_5000_Groups2 <- pam(data_matrix, 2)
plot(pam_result_5000_Groups2)

pam_result_5000_Groups3 <- pam(data_matrix, 3)
plot(pam_result_5000_Groups3)

pam_result_5000_Groups4 <- pam(data_matrix, 4)
plot(pam_result_5000_Groups4)

pam_result_5000_Groups5 <- pam(data_matrix, 5)
plot(pam_result_5000_Groups5)

pam_result_1000_Groups2 <- pam(data_matrix_1000, 2)
plot(pam_result_1000_Groups2)

pam_result_10000_Groups2 <- pam(data_matrix_10000, 2)
plot(pam_result_10000_Groups2)
#Hclust
distance_matrix <- dist(data_matrix)
hclust_result <- hclust(distance_matrix, method = "complete")
plot(hclust_result, labels = FALSE, main = "Hclust dendrogram n = 5000")

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