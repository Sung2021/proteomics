# ============================================================================
# 06_dromics_bmd.R  [Branch B only: dose-response]
# Fit dose-response models and calculate Benchmark Dose (BMD)
# Input : data/processed/protein_matrix.csv
#         data/processed/responding_proteins.csv  (from 05b)
# Output: results/tables/bmd_results.csv
#         results/figures/bmd_heatmap.png
# ============================================================================

library(DRomics)
library(yaml)

cfg <- yaml::read_yaml("config.yaml")
cat("\n=== Step 6: DROmics BMD calculation ===\n")

# --- 1. Load protein matrix & responding proteins ---------------------------
protein_matrix   <- read.csv(cfg$paths$protein_matrix,    row.names = 1)
responding_df    <- read.csv(cfg$paths$responding_proteins)
responding_ids   <- responding_df$protein_id

cat("All proteins    :", nrow(protein_matrix), "\n")
cat("Responding      :", length(responding_ids), "\n")

# Subset to responding proteins only
protein_subset <- protein_matrix[rownames(protein_matrix) %in% responding_ids, ]
cat("Subset rows     :", nrow(protein_subset), "\n")

# --- 2. Build DROmics-format input ------------------------------------------
# DROmics expects: first column = item IDs, remaining columns = dose replicates
# Column names of protein_matrix must encode dose info (e.g. "dose0_1", "dose10_2")
# Dose is extracted from sample metadata in config.yaml
sample_meta <- read.csv(cfg$paths$sample_metadata)
dose_vec    <- sample_meta[[cfg$study$dose_col]][
  match(colnames(protein_subset), sample_meta[[cfg$limpa$sample_col]])
]

# Write temporary DROmics-format file
drdata_df <- cbind(item = rownames(protein_subset), protein_subset)
colnames(drdata_df)[-1] <- dose_vec
tmp_file <- tempfile(fileext = ".csv")
write.csv(drdata_df, tmp_file, row.names = FALSE)

# --- 3. Import data into DROmics --------------------------------------------
o <- continuousomicdata(tmp_file)
cat("DROmics object created\n")
print(o)

# --- 4. Select dose-responding items ----------------------------------------
s <- itemselect(o,
                select.method = cfg$dromics$select_method,
                FDR           = cfg$dromics$fdr)
cat("Items selected:", nrow(s$informativeHits), "\n")

# --- 5. Fit dose-response models --------------------------------------------
cat("Fitting dose-response models (this may take a few minutes)...\n")
f <- drcfit(s, progressbar = FALSE)
cat("Items with convergence:", sum(!is.na(f$fitres$model)), "\n")

# --- 6. Calculate BMD -------------------------------------------------------
bmd_res <- bmdcalc(f,
                   z      = cfg$dromics$bmd_z,
                   x      = cfg$dromics$bmd_x)

bmd_df <- bmd_res$res
cat("\n=== BMD Summary ===\n")
print(summary(bmd_df$BMD.zSD))

dir.create(dirname(cfg$paths$bmd_results), recursive = TRUE, showWarnings = FALSE)
write.csv(bmd_df, cfg$paths$bmd_results, row.names = FALSE)
cat("\n✓ BMD results saved to:", cfg$paths$bmd_results, "\n")

# --- 7. BMD distribution plot -----------------------------------------------
dir.create(dirname(cfg$paths$bmd_heatmap), recursive = TRUE, showWarnings = FALSE)
png(cfg$paths$bmd_heatmap, width = 1000, height = 800, res = 120)
bmdplot(bmd_res, BMDtype = "zSD",
        main = "BMD Distribution (z × SD threshold)")
dev.off()
cat("✓ BMD plot saved to:", cfg$paths$bmd_heatmap, "\n")
