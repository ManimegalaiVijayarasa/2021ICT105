# Leukemia Classification Using Gene Expression Data
# Dataset: Golub Leukemia Dataset
# Analysis: Exploratory Data Analysis and SVM Classification

# 1. Load Required Packages

library(golubEsets)
library(tidyverse)
library(pheatmap)
library(e1071)
library(Biobase)




# 2. Load Golub Leukemia Dataset

# Load the dataset
data(Golub_Merge)

golub_data <- Golub_Merge

# Extract gene expression matrix
expression_matrix <- exprs(golub_data)

# Extract sample information
sample_info <- as.data.frame(pData(golub_data))

# Match sample information with expression matrix columns
sample_info <- sample_info[
  match(
    colnames(expression_matrix),
    rownames(sample_info)
  ),
  ,
  drop = FALSE
]

# Extract leukemia classification labels
group <- sample_info$ALL.AML

# Check dataset dimensions
cat("Dataset Dimensions:\n")
cat("Number of genes:", nrow(expression_matrix), "\n")
cat("Number of samples:", ncol(expression_matrix), "\n")

# Display sample information
cat("\nSample Information:\n")
head(sample_info)

# Display leukemia classification
cat("\nLeukemia Type Distribution:\n")
print(table(group))

# Check sample matching
cat("\nSample matching successful:",
    identical(
      colnames(expression_matrix),
      rownames(sample_info)
    ),
    "\n")




# 3. Data Quality Checking

# Check dataset dimensions
cat("Number of genes:", nrow(expression_matrix), "\n")
cat("Number of samples:", ncol(expression_matrix), "\n")

# Check ALL and AML sample counts
cat("\nLeukemia Type Distribution:\n")
print(table(sample_info$ALL.AML))

# Calculate class percentages
cat("\nClass Percentages:\n")
print(prop.table(table(sample_info$ALL.AML)) * 100)

# Check missing values in sample information
cat("\nMissing Values in Sample Information:\n")
print(colSums(is.na(sample_info)))

# Check missing values in gene expression matrix
cat("\nMissing Values in Expression Matrix:\n")
print(sum(is.na(expression_matrix)))

# Check expression data range
cat("\nExpression Value Range:\n")
print(range(expression_matrix, na.rm = TRUE))




# 4. Exploratory Data Analysis (EDA)


# 4.1 Distribution of ALL and AML Samples

sample_counts <- as.data.frame(
  table(sample_info$ALL.AML)
)

colnames(sample_counts) <- c(
  "Leukemia_Type",
  "Count"
)

# Create bar chart
ggplot(
  sample_counts,
  aes(
    x = Leukemia_Type,
    y = Count,
    fill = Leukemia_Type
  )
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = Count),
    vjust = -0.5,
    size = 5
  ) +
  labs(
    title = "Distribution of ALL and AML Samples",
    x = "Leukemia Type",
    y = "Number of Samples"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

# Save the plot
ggsave(
  "figures/leukemia_sample_distribution.png",
  width = 7,
  height = 5,
  dpi = 300
)





# 5. Gene Expression Distribution Analysis


# Convert gene expression matrix into a data frame

expression_df <- data.frame(
  Sample = rep(
    colnames(expression_matrix),
    each = nrow(expression_matrix)
  ),
  Expression = as.vector(expression_matrix)
)

# Match each sample with its leukemia type

sample_info$Sample_ID <- rownames(sample_info)

expression_df$Leukemia_Type <- sample_info$ALL.AML[
  match(
    expression_df$Sample,
    sample_info$Sample_ID
  )
]

# Create gene expression distribution boxplot

ggplot(
  expression_df,
  aes(
    x = Sample,
    y = Expression,
    fill = Leukemia_Type
  )
) +
  geom_boxplot(outlier.size = 0.2) +
  labs(
    title = "Gene Expression Distribution Across Samples",
    x = "Samples",
    y = "Gene Expression Values",
    fill = "Leukemia Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 5,
      hjust = 1
    )
  )

# Save the plot

ggsave(
  "figures/gene_expression_distribution.png",
  width = 14,
  height = 7,
  dpi = 300
)





# 6. Top 50 Variable Genes Heatmap


# Calculate variance for each gene

gene_variance <- apply(
  expression_matrix,
  1,
  var
)

# Sort genes by variance (highest first)

gene_variance_sorted <- sort(
  gene_variance,
  decreasing = TRUE
)

# Select top 50 most variable genes

top_genes <- names(
  gene_variance_sorted
)[1:50]

# Extract expression values for top 50 genes

top_gene_matrix <- expression_matrix[
  top_genes,
  ,
  drop = FALSE
]

# Standardise gene expression values

scaled_matrix <- t(
  scale(t(top_gene_matrix))
)

# Prepare sample annotations

annotation_col <- data.frame(
  Leukemia_Type = sample_info$ALL.AML
)

rownames(annotation_col) <- colnames(expression_matrix)

# Create heatmap

pheatmap(
  scaled_matrix,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 50 Most Variable Genes in ALL and AML",
  filename = "figures/top_50_variable_genes_heatmap.png",
  width = 12,
  height = 10
)






# 7. Differential Gene Expression Analysis


# Calculate mean expression for ALL samples

mean_ALL <- rowMeans(
  expression_matrix[, group == "ALL", drop = FALSE]
)

# Calculate mean expression for AML samples

mean_AML <- rowMeans(
  expression_matrix[, group == "AML", drop = FALSE]
)

# Perform Welch's t-test for each gene

p_values <- apply(
  expression_matrix,
  1,
  function(x) {
    t.test(
      x[group == "ALL"],
      x[group == "AML"]
    )$p.value
  }
)

# Create differential expression results

differential_results <- data.frame(
  Gene = rownames(expression_matrix),
  Mean_ALL = mean_ALL,
  Mean_AML = mean_AML,
  Difference = mean_AML - mean_ALL,
  P_value = p_values
)

# Adjust p-values using Benjamini-Hochberg method

differential_results$Adjusted_P_value <- p.adjust(
  differential_results$P_value,
  method = "BH"
)

# Sort genes by adjusted p-value

differential_results <- differential_results[
  order(differential_results$Adjusted_P_value),
]

# Display top 20 genes

head(differential_results, 20)

# Identify statistically significant genes

significant_genes <- subset(
  differential_results,
  Adjusted_P_value < 0.05
)

# Display number of significant genes

cat(
  "Number of significant genes:",
  nrow(significant_genes),
  "\n"
)

# Display first 10 significant genes

head(significant_genes, 10)

# Save differential expression results

write.csv(
  differential_results,
  "output/differential_gene_results.csv",
  row.names = FALSE
)





# 8. Volcano Plot


# Classify genes based on adjusted p-values

differential_results$Significance <- ifelse(
  differential_results$Adjusted_P_value < 0.05,
  "Significant",
  "Not Significant"
)

# Create volcano plot

ggplot(
  differential_results,
  aes(
    x = Difference,
    y = -log10(
      pmax(
        Adjusted_P_value,
        .Machine$double.xmin
      )
    ),
    color = Significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.2
  ) +
  labs(
    title = "Volcano Plot: ALL vs AML",
    x = "Mean Expression Difference (AML - ALL)",
    y = "-Log10 Adjusted P-value",
    color = "Gene Significance"
  ) +
  theme_minimal()

# Save volcano plot

ggsave(
  "figures/volcano_plot.png",
  width = 9,
  height = 6,
  dpi = 300
)







# 9. SVM Classification


# 9.1 Split dataset into training and testing sets

set.seed(123)

# Select training samples separately for each class
train_ALL <- sample(
  which(group == "ALL"),
  floor(0.7 * sum(group == "ALL"))
)

train_AML <- sample(
  which(group == "AML"),
  floor(0.7 * sum(group == "AML"))
)

# Combine training samples
train_index <- c(train_ALL, train_AML)

# Remaining samples for testing
test_index <- setdiff(
  seq_along(group),
  train_index
)

# Check sample distribution
cat("Training sample distribution:\n")
print(table(group[train_index]))

cat("\nTesting sample distribution:\n")
print(table(group[test_index]))


# 9.2 Select the top 50 genes using training data only

training_p_values <- apply(
  expression_matrix[, train_index, drop = FALSE],
  1,
  function(x) {
    t.test(
      x[group[train_index] == "ALL"],
      x[group[train_index] == "AML"]
    )$p.value
  }
)

# Select 50 genes with the smallest p-values
top_genes_ml <- names(
  sort(training_p_values)
)[1:50]

cat("\nNumber of selected genes:", length(top_genes_ml), "\n")


# 9.3 Prepare training and testing data

train_data <- t(
  expression_matrix[
    top_genes_ml,
    train_index,
    drop = FALSE
  ]
)

test_data <- t(
  expression_matrix[
    top_genes_ml,
    test_index,
    drop = FALSE
  ]
)

# Convert labels into factors
train_labels <- factor(
  group[train_index],
  levels = c("ALL", "AML")
)

test_labels <- factor(
  group[test_index],
  levels = c("ALL", "AML")
)


# 9.4 Train SVM model

svm_model <- svm(
  x = train_data,
  y = train_labels,
  kernel = "radial",
  scale = TRUE
)

# Display model
print(svm_model)


# 9.5 Predict test samples

predictions <- predict(
  svm_model,
  newdata = test_data
)

# Compare actual and predicted classes
results <- data.frame(
  Actual = test_labels,
  Predicted = predictions
)

print(results)


# 9.6 Create confusion matrix

conf_matrix <- table(
  Actual = test_labels,
  Predicted = predictions
)

print(conf_matrix)


# 9.7 Calculate accuracy

accuracy <- sum(diag(conf_matrix)) /
  sum(conf_matrix) * 100

cat(
  "\nModel Accuracy:",
  round(accuracy, 2),
  "%\n"
)


# 9.8 Save prediction results

write.csv(
  results,
  "output/svm_prediction_results.csv",
  row.names = FALSE
)


# 9.9 Calculate precision, recall and F1-score

precision <- diag(conf_matrix) /
  colSums(conf_matrix)

recall <- diag(conf_matrix) /
  rowSums(conf_matrix)

f1_score <- 2 * (precision * recall) /
  (precision + recall)

performance <- data.frame(
  Class = names(precision),
  Precision = as.numeric(precision),
  Recall = as.numeric(recall),
  F1_Score = as.numeric(f1_score)
)

print(performance)


# 9.10 Calculate balanced accuracy

balanced_accuracy <- mean(recall) * 100

cat(
  "\nBalanced Accuracy:",
  round(balanced_accuracy, 2),
  "%\n"
)


# Save performance metrics

write.csv(
  performance,
  "output/svm_performance_metrics.csv",
  row.names = FALSE
)






# 10. 5-Fold Cross-Validation


# 10.1 Create stratified folds

set.seed(123)

fold_id <- integer(length(group))

# Assign ALL samples to five folds
all_samples <- sample(which(group == "ALL"))

fold_id[all_samples] <- rep(
  1:5,
  length.out = length(all_samples)
)

# Assign AML samples to five folds
aml_samples <- sample(which(group == "AML"))

fold_id[aml_samples] <- rep(
  1:5,
  length.out = length(aml_samples)
)

# Check class distribution in each fold
for (k in 1:5) {
  cat("\nFold", k, "\n")
  print(table(group[fold_id == k]))
}


# 10.2 Perform 5-fold cross-validation

cv_predictions <- rep(
  NA_character_,
  length(group)
)

fold_accuracy <- numeric(5)

for (k in 1:5) {
  
  # Separate training and testing samples
  train_idx <- which(fold_id != k)
  test_idx <- which(fold_id == k)
  
  # Select genes using training samples only
  training_p_values <- apply(
    expression_matrix[, train_idx, drop = FALSE],
    1,
    function(x) {
      t.test(
        x[group[train_idx] == "ALL"],
        x[group[train_idx] == "AML"]
      )$p.value
    }
  )
  
  # Select top 50 genes
  selected_genes <- names(
    sort(training_p_values)
  )[1:50]
  
  # Prepare training and testing data
  x_train <- t(
    expression_matrix[
      selected_genes,
      train_idx,
      drop = FALSE
    ]
  )
  
  x_test <- t(
    expression_matrix[
      selected_genes,
      test_idx,
      drop = FALSE
    ]
  )
  
  # Prepare class labels
  y_train <- factor(
    group[train_idx],
    levels = c("ALL", "AML")
  )
  
  y_test <- factor(
    group[test_idx],
    levels = c("ALL", "AML")
  )
  
  # Train SVM model
  model <- svm(
    x = x_train,
    y = y_train,
    kernel = "radial",
    scale = TRUE
  )
  
  # Predict testing samples
  predicted <- predict(
    model,
    newdata = x_test
  )
  
  # Store predictions
  cv_predictions[test_idx] <- as.character(predicted)
  
  # Calculate fold accuracy
  fold_accuracy[k] <- mean(
    predicted == y_test
  ) * 100
  
  cat(
    "Fold", k,
    "- Accuracy:",
    round(fold_accuracy[k], 2),
    "%\n"
  )
}


# 10.3 Create overall confusion matrix

cv_conf_matrix <- table(
  Actual = factor(group, levels = c("ALL", "AML")),
  Predicted = factor(
    cv_predictions,
    levels = c("ALL", "AML")
  )
)

print(cv_conf_matrix)


# 10.4 Calculate overall accuracy

cv_accuracy <- sum(
  diag(cv_conf_matrix)
) / sum(cv_conf_matrix) * 100

cat(
  "\n5-Fold Cross-Validation Accuracy:",
  round(cv_accuracy, 2),
  "%\n"
)


# 10.5 Calculate average fold accuracy

cat(
  "Average Fold Accuracy:",
  round(mean(fold_accuracy), 2),
  "%\n"
)

# Calculate standard deviation

cat(
  "Accuracy Standard Deviation:",
  round(sd(fold_accuracy), 2),
  "%\n"
)


# 10.6 Calculate precision, recall and F1-score

cv_precision <- diag(cv_conf_matrix) /
  colSums(cv_conf_matrix)

cv_recall <- diag(cv_conf_matrix) /
  rowSums(cv_conf_matrix)

cv_f1 <- 2 * (cv_precision * cv_recall) /
  (cv_precision + cv_recall)

cv_performance <- data.frame(
  Class = names(cv_precision),
  Precision = as.numeric(cv_precision),
  Recall = as.numeric(cv_recall),
  F1_Score = as.numeric(cv_f1)
)

print(cv_performance)


# 10.7 Calculate balanced accuracy

cv_balanced_accuracy <- mean(cv_recall) * 100

cat(
  "Balanced Accuracy:",
  round(cv_balanced_accuracy, 2),
  "%\n"
)


# 10.8 Save cross-validation results

cv_results <- data.frame(
  Actual = group,
  Predicted = cv_predictions,
  Fold = fold_id
)

write.csv(
  cv_results,
  "output/svm_5fold_cv_results.csv",
  row.names = FALSE
)

# Save fold accuracy

write.csv(
  data.frame(
    Fold = 1:5,
    Accuracy = fold_accuracy
  ),
  "output/svm_5fold_accuracy.csv",
  row.names = FALSE
)

# Save performance metrics

write.csv(
  cv_performance,
  "output/svm_cv_performance_metrics.csv",
  row.names = FALSE
)

