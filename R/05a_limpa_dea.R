# ============================================================================
# 05a_limpa_dea.R  [Branch A: group comparison — FINAL OUTPUT]
# Peptide → protein quantification (DPC-quant) + limma DEA
# Output: results/tables/dea_results.csv, results/figures/volcano.png
# ============================================================================

library(limpa)
library(limma)
library(yaml)

cfg <- yaml::read_yaml("config.yaml")
cat("\n=== Step 5a: DEA (group comparison) ===\n")

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

sample_map <- unique(peptide_dt[, c(cfg$limpa$sample_col, cfg$study$group_col)])
L <- limpa::addTargets(L, targets = sample_map,
                       sample_id    = cfg$limpa$sample_col,
                       condition_id = cfg$study$group_col)

# --- 3. DPC estimation & protein quantification ----------------------------
cat("Estimating DPC...\n")
L <- limpa::estimateDPC(L)
L <- limpa::DPCquant(L)
cat("Proteins quantified:", nrow(limpa::proteinExpression(L)), "\n")

# --- 4. limma DEA -----------------------------------------------------------
targets <- limpa::targets(L)
design  <- model.matrix(~ 0 + targets[[cfg$study$group_col]])
colnames(design) <- levels(factor(targets[[cfg$study$group_col]]))

fit  <- limpa::limmaFit(L, design = design)
fit  <- limma::eBayes(fit)

contrast_str <- paste(cfg$study$treatment_group, "-", cfg$study$reference_group)
contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)
fit2 <- limma::contrasts.fit(fit, contrast_mat)
fit2 <- limma::eBayes(fit2)

all_res <- limma::topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# --- 5. Save results --------------------------------------------------------
dir.create(dirname(cfg$paths$dea_results), recursive = TRUE, showWarnings = FALSE)
write.csv(all_res, cfg$paths$dea_results, row.names = TRUE)
cat("\n✓ DEA results saved to:", cfg$paths$dea_results, "\n")

# --- 6. Volcano plot --------------------------------------------------------
sig <- cfg$thresholds$adj_pval
fc  <- cfg$thresholds$log2fc
sig_proteins <- subset(all_res, adj.P.Val < sig & abs(logFC) > fc)

dir.create(dirname(cfg$paths$volcano_plot), recursive = TRUE, showWarnings = FALSE)
png(cfg$paths$volcano_plot, width = 900, height = 800, res = 120)
plot(all_res$logFC, -log10(all_res$adj.P.Val),
     main = paste0("Volcano Plot (", cfg$study$treatment_group,
                   " vs ", cfg$study$reference_group, ")"),
     xlab = "Log2 Fold Change", ylab = "-log10(adj.P.Val)",
     pch = 20, cex = 0.6, col = "grey60")
points(sig_proteins$logFC, -log10(sig_proteins$adj.P.Val),
       col = "firebrick", pch = 20, cex = 0.8)
abline(h = -log10(sig), v = c(-fc, fc), col = "steelblue", lty = 2)
legend("topright", legend = c(
  paste0("Significant (n=", nrow(sig_proteins), ")"),
  paste0("adj.P.Val < ", sig, "  |logFC| > ", fc)
), bty = "n", text.col = c("firebrick", "steelblue"))
dev.off()

cat("✓ Volcano plot saved to:", cfg$paths$volcano_plot, "\n")
cat("Significant proteins:", nrow(sig_proteins), "\n")
