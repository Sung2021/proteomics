# ============================================================================
# 05b_limpa_prefilter.R  [Branch B: dose-response — FILTER STEP]
# Identify proteins that respond to dose treatment via limma DEA.
# Responding proteins are passed to 06_dromics_bmd.R.
# Output: data/processed/responding_proteins.csv
# ============================================================================

library(limpa)
library(limma)
library(yaml)

cfg <- yaml::read_yaml("config.yaml")
cat("\n=== Step 5b: Prefilter — dose-responding proteins ===\n")

# --- 1. Load peptide data ---------------------------------------------------
peptide_dt <- read.csv(cfg$paths$peptide_data)
cat("Peptide rows:", nrow(peptide_dt), "\n")

# --- 2. Build limpa object --------------------------------------------------
L <- limpa::newLimpa(
  peptide_dt,
  peptide_id   = cfg$limpa$peptide_col,
  protein_id   = cfg$limpa$protein_col,
  sample_id    = cfg$limpa$sample_col,
  intensity    = cfg$limpa$intensity_col
)

# In dose-response design use dose_col groups for the linear model
# (i.e., treat each dose as a separate group to call responding proteins)
sample_map <- unique(peptide_dt[, c(cfg$limpa$sample_col, cfg$study$group_col)])
L <- limpa::addTargets(L, targets = sample_map,
                       sample_id    = cfg$limpa$sample_col,
                       condition_id = cfg$study$group_col)

# --- 3. DPC & protein quantification ----------------------------------------
cat("Estimating DPC...\n")
L <- limpa::estimateDPC(L)
L <- limpa::DPCquant(L)

protein_expr <- limpa::proteinExpression(L)
cat("Proteins quantified:", nrow(protein_expr), "\n")

# Export full protein matrix for DROmics
dir.create(dirname(cfg$paths$protein_matrix), recursive = TRUE, showWarnings = FALSE)
write.csv(protein_expr, cfg$paths$protein_matrix, row.names = TRUE)
cat("✓ Protein matrix saved to:", cfg$paths$protein_matrix, "\n")

# --- 4. Filter responding proteins (any dose vs control) --------------------
targets <- limpa::targets(L)
design  <- model.matrix(~ 0 + targets[[cfg$study$group_col]])
colnames(design) <- levels(factor(targets[[cfg$study$group_col]]))

fit <- limpa::limmaFit(L, design = design)
fit <- limma::eBayes(fit)

# Build one contrast per non-reference group
ref <- cfg$study$reference_group
non_ref <- setdiff(colnames(design), ref)
contrast_list <- paste(non_ref, "-", ref)
contrast_mat  <- limma::makeContrasts(contrasts = contrast_list, levels = design)
fit2 <- limma::contrasts.fit(fit, contrast_mat)
fit2 <- limma::eBayes(fit2)

# A protein is "responding" if significant in ANY contrast
all_res <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
responding <- rownames(subset(all_res,
                               adj.P.Val < cfg$thresholds$prefilter_adj_pval))

cat("Responding proteins:", length(responding),
    "/ Total:", nrow(protein_expr), "\n")

if (length(responding) == 0)
  warning("No responding proteins found. Consider relaxing prefilter_adj_pval in config.yaml.")

# --- 5. Save responding protein list ----------------------------------------
dir.create(dirname(cfg$paths$responding_proteins), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(protein_id = responding),
          cfg$paths$responding_proteins, row.names = FALSE)
cat("✓ Responding proteins saved to:", cfg$paths$responding_proteins, "\n")
cat("→ Pass this list to 06_dromics_bmd.R\n")
