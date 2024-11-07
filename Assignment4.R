# The following links/resources were used to create the code in this file:
# https://seandavi.github.io/ITR/machine_learning_mlr3.html#Background


install.packages("readr")
install.packages("ggplot2")
install.packages("mlr3verse")
install.packages("mlr3learners")
install.packages("pROC")
install.packages("ranger")
install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("magrittr")
install.packages("dplyr")
install.packages("ConsensusClusterPlus")
install.packages("cluster")
install.packages("factoextra")
install.packages("pheatmap")
install.packages("ggaluvial")
install.packages("tidymodels")
BiocManager::install("GEOquery")
install.packages("kknn")
install.packages("e1071")



library(DESeq2)
library(readr)
library(magrittr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(mlr3verse)
library(GEOquery)
library(mlr3learners) # for K nearest neighbors
library(ranger) # for Random forest
library(pROC)
library(factoextra)
library(randomForest)
library(tidymodels)
library(kknn)
library(e1071)


data <- read_tsv("./data/SRP018853/SRP018853.tsv")
metadata <- read_tsv("./data/SRP018853/metadata_SRP018853.tsv")

variance_result <- apply(data[2, 2:81], 1, var, na.rm = TRUE)
variance_result <- apply(data[2:nrow(data), 2:81], 1, var, na.rm = TRUE)
top_geneNumbers <- order(variance_result, decreasing = TRUE)[1:5000]
top_5000 <- data[top_geneNumbers, ]

metadata$TestGroups <- ifelse(stringr::str_starts(metadata$refinebio_title, 'M'), "Healthy", "Pre-T1D")

# Transpose data
transp_5000 <- as.data.frame(t(top_5000[ , -1]))

colnames(transp_5000) <- top_5000$Gene
transp_5000$sample_id <- colnames(data)[-1]

# Select only the necessary columns from the metadata
relevant_columns = c("refinebio_accession_code", "TestGroups")  # Adjust based on the actual names
clean_metadata = metadata[, relevant_columns]

# Merge gene expression data with metadata using sample IDs
merged_data <- merge(transp_5000, clean_metadata, by.x = "sample_id", by.y = "refinebio_accession_code")

# Set health_status column as factor for classification task
merged_data$TestGroups <- as.factor(merged_data$TestGroups)

# Remove 'sample_id' column from merged_data
merged_data <- merged_data[, !(names(merged_data) == "sample_id")]

# Creating the task
task = as_task_classif(merged_data,target='TestGroups')

set.seed(24)
train_set = sample(task$row_ids, 0.67 * task$nrow)
test_set = setdiff(task$row_ids, train_set)



# K-nearest neighbors (Test Groups from Assignment 1)
learner = lrn("classif.kknn")
learner$train(task, row_ids = train_set)
learner$model

pred_train = learner$predict(task, row_ids=train_set)
pred_test = learner$predict(task, row_ids=test_set)

pred_train$confusion

measures = msrs(c('classif.acc'))
pred_train$score(measures)

pred_test$confusion
pred_test$score(measures)



# Random forest (Test Groups from Assignment 1)
learner = lrn("classif.ranger", importance = "impurity")
learner$train(task, row_ids = train_set)
learner$model

pred_train = learner$predict(task, row_ids=train_set)
pred_test = learner$predict(task, row_ids=test_set)

pred_train$confusion

measures = msrs(c('classif.acc'))
pred_train$score(measures)

pred_test$confusion
pred_test$score(measures)

variab_filter = flt("importance", learner = learner)
variab_filter$calculate(task)
head(as.data.table(variab_filter), 10)



# Support vector machine (Test Groups from Assignment 1)
learner = lrn("classif.svm")
learner$train(task, row_ids = train_set)
learner$model

pred_train = learner$predict(task, row_ids=train_set)
pred_test = learner$predict(task, row_ids=test_set)

pred_train$confusion

measures = msrs(c('classif.acc'))
pred_train$score(measures)

pred_test$confusion
pred_test$score(measures)



pred_matrix <- matrix(NA, nrow = length(test_set), ncol = 3)  # 3 for KNN, Random Forest, and SVM
colnames(pred_matrix) <- c("K-NN", "RandomForest", "SVM")
