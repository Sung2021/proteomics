# ============================================================================
# 03_qc.R - Quality Control Plots
# ============================================================================

cat("\n=== Step 3: Quality Control ===\n")

library(Spectra)

# Load clean spectra
sp <- readRDS("output/processed/clean_spectra.rds")
cat("Loaded", length(sp), "spectra\n")

# Extract MS1 spectra
ms1 <- filterMsLevel(sp, 1)
cat("MS1 spectra:", length(ms1), "\n")

# Extract TIC and retention time
tic <- spectraData(ms1)$totIonCurrent
rt  <- spectraData(ms1)$rtime

# Create TIC plot
cat("\nGenerating TIC plot...\n")
png("output/qc_plots/tic.png", width = 1200, height = 800, res = 120)
plot(rt, tic, 
     type = "l", 
     main = "Total Ion Current (TIC)",
     xlab = "Retention Time (s)",
     ylab = "Total Ion Current",
     col = "blue",
     lwd = 2)
grid()
dev.off()

cat("✓ TIC plot saved to: output/qc_plots/tic.png\n")

# Additional QC metrics
cat("\n=== QC Metrics ===\n")
cat("TIC range:", format(range(tic), scientific = TRUE), "\n")
cat("RT range:", round(range(rt), 2), "seconds\n")
cat("Mean TIC:", format(mean(tic), scientific = TRUE), "\n")
cat("Median TIC:", format(median(tic), scientific = TRUE), "\n")

# Optional: MS level distribution plot
ms_levels <- table(msLevel(sp))
png("output/qc_plots/ms_level_distribution.png", 
    width = 800, height = 600, res = 120)
barplot(ms_levels,
        main = "MS Level Distribution",
        xlab = "MS Level",
        ylab = "Number of Spectra",
        col = c("lightblue", "lightcoral"),
        border = "black")
dev.off()

cat("✓ MS level distribution plot saved\n")
