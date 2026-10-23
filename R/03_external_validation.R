# 03_external_validation.R
# Train the final SVM on all 226 paired samples, then test it on tumors from
# patients it has never seen. Since every external sample is a tumor, the
# metric here is sensitivity (share called "Tumor"), broken down by tissue
# source site and by PAM50 subtype.

library(SummarizedExperiment)
library(TCGAbiolinks)
library(caret)
library(ggplot2)

paired   <- readRDS("data/tcga_brca_paired.rds")
external <- readRDS("data/tcga_brca_external_tumors.rds")

log_cpm <- function(counts) log2(t(t(counts) / colSums(counts)) * 1e6 + 1)

counts <- assays(paired)[["unstranded"]]
cpm    <- t(t(counts) / colSums(counts)) * 1e6
logcpm <- log2(cpm + 1)
label  <- factor(ifelse(paired$sample_type == "Primary Tumor", "Tumor", "Normal"),
                 levels = c("Normal", "Tumor"))

# ---- Final model on all paired samples -----------------------------------
expressed <- rowMeans(cpm > 1) >= 0.2
v         <- apply(logcpm[expressed, ], 1, var)
genes     <- names(sort(v, decreasing = TRUE))[1:1000]

ctrl     <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
svm_grid <- expand.grid(sigma = c(0.0005, 0.001, 0.005), C = c(0.5, 1, 2))
set.seed(123)
final_fit <- train(label ~ ., data = data.frame(t(logcpm[genes, ]), label = label),
                   method = "svmRadial", trControl = ctrl, tuneGrid = svm_grid, metric = "ROC")
print(final_fit$bestTune)

# ---- Predict external tumors (normalized the same way, same 1,000 genes) ----
ext_logcpm <- log_cpm(assays(external)[["unstranded"]])
ext_pred   <- predict(final_fit, newdata = data.frame(t(ext_logcpm[genes, ])))

normal_sites <- unique(substr(colnames(paired)[label == "Normal"], 6, 7))
res <- data.frame(barcode = colnames(external),
                  patient = substr(colnames(external), 1, 12),
                  site    = substr(colnames(external), 6, 7),
                  correct = ext_pred == "Tumor")
res$site_had_normals <- res$site %in% normal_sites

cat("External tumors:", nrow(res), " called Tumor:", sum(res$correct),
    " sensitivity:", round(mean(res$correct), 3), "\n")

by_normals <- aggregate(correct ~ site_had_normals, data = res,
                        FUN = function(x) c(n = length(x), sensitivity = round(mean(x), 3)))
by_normals <- do.call(data.frame, by_normals)
print(by_normals)
write.csv(by_normals, "results/external_by_normal_site.csv", row.names = FALSE)

by_site <- aggregate(correct ~ site, data = res,
                     FUN = function(x) c(n = length(x), sensitivity = round(mean(x), 3)))
by_site <- do.call(data.frame, by_site)
names(by_site) <- c("site", "n", "sensitivity")
by_site <- by_site[order(-by_site$n), ]
print(by_site[by_site$n >= 20, ], row.names = FALSE)
write.csv(by_site, "results/external_by_site.csv", row.names = FALSE)

# ---- Which tumors were missed? PAM50 subtype -----------------------------
subtypes  <- TCGAquery_subtype(tumor = "BRCA")
res$pam50 <- subtypes$BRCA_Subtype_PAM50[match(res$patient, subtypes$patient)]
by_pam50 <- aggregate(correct ~ pam50, data = res,
                      FUN = function(x) c(n = length(x), called_tumor = sum(x), sensitivity = round(mean(x), 3)))
by_pam50 <- do.call(data.frame, by_pam50)
names(by_pam50) <- c("pam50", "n", "called_tumor", "sensitivity")
print(by_pam50, row.names = FALSE)
write.csv(by_pam50, "results/external_by_pam50.csv", row.names = FALSE)

p <- ggplot(by_site[by_site$n >= 20, ], aes(reorder(site, sensitivity), sensitivity)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip(ylim = c(0.8, 1.02)) +
  labs(title = "External tumors called Tumor, by tissue source site",
       subtitle = "Sites with 20+ tumors; label = number of tumors",
       x = "Tissue source site", y = "Sensitivity") +
  theme_minimal()
ggsave("figures/external_by_site.png", p, width = 7, height = 5, dpi = 150)

write.csv(res, "results/external_predictions.csv", row.names = FALSE)
