# 01_get_data.R
# Select samples from TCGA-BRCA, check where they come from, build matched
# tumor-normal pairs, and download/prepare the RNA-seq counts.
#
# Outputs (not tracked in git, recreated by this script):
#   data/tcga_brca_paired.rds            113 matched tumor-normal pairs (226 samples)
#   data/tcga_brca_external_tumors.rds   one tumor per remaining patient (external test set)

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- Find all BRCA RNA-seq samples ---------------------------------------
query_all <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
cases <- getResults(query_all)$cases

all_tumors  <- TCGAquery_SampleTypes(barcode = cases, typesample = "TP")  # primary tumor (01)
all_normals <- TCGAquery_SampleTypes(barcode = cases, typesample = "NT")  # solid tissue normal (11)
cat("Tumor samples:", length(all_tumors), " Normal samples:", length(all_normals), "\n")

# Barcode parts: TCGA-SS-PPPP-TTV...  SS = tissue source site, first 12 chars = patient
site    <- function(x) substr(x, 6, 7)
patient <- function(x) substr(x, 1, 12)

# ---- Design check: where do tumors and normals come from? ----------------
# The original analysis used the first 113 tumors vs all 113 normals.
first_113 <- all_tumors[1:113]
site_check <- data.frame(
  group     = c("First 113 tumors", "All 113 normals"),
  n_sites   = c(length(unique(site(first_113))), length(unique(site(all_normals)))),
  top_site  = c(names(which.max(table(site(first_113)))), names(which.max(table(site(all_normals))))),
  top_share = c(max(table(site(first_113))), max(table(site(all_normals)))) / 113
)
print(site_check)
write.csv(site_check, "results/site_check_original_design.csv", row.names = FALSE)

# ---- Matched pairs: one tumor and one normal from the same patient --------
pick_one <- function(barcodes, p) {
  b <- barcodes[patient(barcodes) == p]
  a <- b[substr(b, 14, 16) == "01A"]           # prefer the first tumor portion
  if (length(a) > 0) a[1] else b[1]
}
paired <- intersect(patient(all_normals), patient(all_tumors))
cat("Patients with both a tumor and a normal sample:", length(paired), "\n")

tumor_pick  <- unname(sapply(paired, function(p) pick_one(all_tumors, p)))
normal_pick <- unname(sapply(paired, function(p) all_normals[patient(all_normals) == p][1]))
stopifnot(identical(sort(site(tumor_pick)), sort(site(normal_pick))))  # same site mix in both groups

# ---- External test set: tumors from patients NOT in the pairs -------------
unused <- all_tumors[!patient(all_tumors) %in% paired]
external_pick <- unname(sapply(unique(patient(unused)), function(p) pick_one(unused, p)))
cat("External tumors:", length(external_pick),
    "from", length(unique(site(external_pick))), "tissue source sites\n")

# ---- Download and prepare -------------------------------------------------
get_counts <- function(barcodes, file) {
  if (file.exists(file)) return(invisible(NULL))
  q <- GDCquery(project = "TCGA-BRCA",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = "STAR - Counts",
                barcode = barcodes)
  GDCdownload(q, files.per.chunk = 50)
  saveRDS(GDCprepare(q), file)
}
get_counts(c(tumor_pick, normal_pick), "data/tcga_brca_paired.rds")
get_counts(external_pick,              "data/tcga_brca_external_tumors.rds")
