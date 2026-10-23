# 02_train_and_validate.R
# Train tumor vs normal classifiers on the 113 matched pairs without leakage:
#   - split by patient, so a patient's tumor and normal never end up on different sides
#   - normalize each sample on its own (log-CPM)
#   - choose genes using training samples only
# Then estimate performance with a held-out test set, nested cross-validation,
# and a permutation test.

library(SummarizedExperiment)
library(caret)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

data   <- readRDS("data/tcga_brca_paired.rds")
counts <- assays(data)[["unstranded"]]

label   <- factor(ifelse(data$sample_type == "Primary Tumor", "Tumor", "Normal"),
                  levels = c("Normal", "Tumor"))
patient <- substr(colnames(data), 1, 12)
print(table(label))

# ---- Normalization: log2 counts per million (per sample, nothing shared) --
cpm    <- t(t(counts) / colSums(counts)) * 1e6
logcpm <- log2(cpm + 1)

# Gene selection: expressed (CPM > 1 in >= 20% of samples), then top 1,000 by variance.
# Always called on training samples only.
select_genes <- function(idx, n = 1000) {
  expressed <- rowMeans(cpm[, idx] > 1) >= 0.2
  v <- apply(logcpm[expressed, idx], 1, var)
  names(sort(v, decreasing = TRUE))[1:n]
}
as_df <- function(genes, idx) data.frame(t(logcpm[genes, idx]), label = label[idx])

# ---- Patient-level 70/30 split -------------------------------------------
set.seed(123)
train_patients <- sample(unique(patient), size = round(0.7 * length(unique(patient))))
is_train <- patient %in% train_patients
print(table(is_train, label))

genes    <- select_genes(is_train)
train_df <- as_df(genes, is_train)
test_df  <- as_df(genes, !is_train)

# ---- PCA (fit on training, test projected) --------------------------------
pca      <- prcomp(train_df[, -ncol(train_df)], scale. = TRUE)
test_pcs <- predict(pca, newdata = test_df[, -ncol(test_df)])
var_expl <- summary(pca)$importance[2, 1:3]
print(var_expl)

pca_df <- rbind(
  data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], label = train_df$label, set = "Train"),
  data.frame(PC1 = test_pcs[, 1], PC2 = test_pcs[, 2], label = test_df$label, set = "Test")
)
p <- ggplot(pca_df, aes(PC1, PC2, color = label, shape = set)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "PCA of top 1,000 genes (fit on training patients only)",
       x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2])) +
  theme_minimal()
ggsave("figures/pca_train_test.png", p, width = 7, height = 5, dpi = 150)

# ---- Models ---------------------------------------------------------------
ctrl     <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
svm_grid <- expand.grid(sigma = c(0.0005, 0.001, 0.005), C = c(0.5, 1, 2))

set.seed(123)
svm_fit <- train(label ~ ., data = train_df, method = "svmRadial",
                 trControl = ctrl, tuneGrid = svm_grid, metric = "ROC")
print(svm_fit$bestTune)
svm_pred <- predict(svm_fit, test_df)

# Logistic regression on the top 30 PCs (1,000 genes is too many predictors for 158 samples)
pc_train <- data.frame(pca$x[, 1:30], label = train_df$label)
pc_test  <- data.frame(test_pcs[, 1:30], label = test_df$label)
log_fit  <- glm(label ~ ., data = pc_train, family = binomial)
log_pred <- factor(ifelse(predict(log_fit, pc_test, type = "response") > 0.5, "Tumor", "Normal"),
                   levels = levels(label))

set.seed(123)
rf_fit  <- train(label ~ ., data = train_df, method = "rf",
                 trControl = ctrl, tuneLength = 3, metric = "ROC")
rf_pred <- predict(rf_fit, test_df)

summarize <- function(model, pred) {
  cm <- confusionMatrix(pred, test_df$label, positive = "Tumor")
  data.frame(model = model,
             correct = sum(pred == test_df$label), n = nrow(test_df),
             accuracy = round(cm$overall["Accuracy"], 3),
             sensitivity_tumor = round(cm$byClass["Sensitivity"], 3),
             specificity = round(cm$byClass["Specificity"], 3),
             false_neg = cm$table["Normal", "Tumor"], false_pos = cm$table["Tumor", "Normal"])
}
test_results <- rbind(summarize("SVM (radial)", svm_pred),
                      summarize("Logistic regression (30 PCs)", log_pred),
                      summarize("Random forest", rf_pred))
rownames(test_results) <- NULL
print(test_results)
write.csv(test_results, "results/test_set_results.csv", row.names = FALSE)

# ---- Nested cross-validation, grouped by patient ---------------------------
# Gene selection and tuning are repeated inside each outer fold.
set.seed(123)
folds  <- groupKFold(patient, k = 5)
nested <- data.frame(fold = seq_along(folds), accuracy = NA)
for (i in seq_along(folds)) {
  tr <- folds[[i]]
  te <- setdiff(seq_along(label), tr)
  g  <- select_genes(tr)
  set.seed(123)
  fit <- train(label ~ ., data = as_df(g, tr), method = "svmRadial",
               trControl = ctrl, tuneGrid = svm_grid, metric = "ROC")
  nested$accuracy[i] <- mean(predict(fit, as_df(g, te)) == label[te])
  cat("Fold", i, "accuracy:", round(nested$accuracy[i], 3), "\n")
}
cat("Nested CV mean accuracy:", round(mean(nested$accuracy), 3), "\n")
write.csv(nested, "results/nested_cv.csv", row.names = FALSE)

# ---- Permutation test: shuffle training labels, keep tuned settings ----------
set.seed(123)
fixed    <- trainControl(method = "none", classProbs = TRUE)
perm_acc <- replicate(50, {
  shuffled <- train_df
  shuffled$label <- sample(shuffled$label)
  fit <- train(label ~ ., data = shuffled, method = "svmRadial",
               trControl = fixed, tuneGrid = svm_fit$bestTune)
  mean(predict(fit, test_df) == test_df$label)
})
real_acc <- mean(svm_pred == test_df$label)
cat("Permutation mean:", round(mean(perm_acc), 3),
    " range:", round(range(perm_acc), 3),
    " permutations >= real accuracy:", sum(perm_acc >= real_acc), "of 50\n")
write.csv(data.frame(permutation = 1:50, accuracy = perm_acc), "results/permutation_test.csv", row.names = FALSE)

png("figures/permutation_test.png", width = 900, height = 650, res = 150)
hist(perm_acc, breaks = seq(0, 1, 0.1), xlim = c(0, 1), col = "grey80",
     main = "SVM accuracy with shuffled labels (50 permutations)", xlab = "Test accuracy")
abline(v = real_acc, col = "red", lwd = 2)
legend("topleft", legend = "Real labels", col = "red", lwd = 2, bty = "n")
dev.off()
