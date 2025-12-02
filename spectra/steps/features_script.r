# ============================================================================
# 04_features.R - Extract MS2 Features
# ============================================================================

cat("\n=== Step 4: Feature Extraction ===\n")

library(Spectra)

# Load clean spectra
sp <- readRDS("output/processed/clean_spectra.rds")
cat("Loaded", length(sp), "spectra\n")

# Extract MS2 spectra
ms2 <- filterMsLevel(sp, 2)
cat("MS2 spectra:", length(ms2), "\n")

if (length(ms2) == 0) {
  stop("No MS2 spectra found. Check your data.")
}

# Extract features
cat("\nExtracting features...\n")

# Base peak intensity (max intensity per spectrum)
bpi <- sapply(intensity(ms2), max)

# Total ion current (sum of intensities)
tic2 <- sapply(intensity(ms2), sum)

# Number of peaks per spectrum
npeaks <- sapply(peaksData(ms2), nrow)

# Precursor information
precursor_mz <- precursorMz(ms2)
precursor_int <- precursorIntensity(ms2)

# Create feature table
df <- data.frame(
  file            = dataOrigin(ms2),
  scan_number     = acquisitionNum(ms2),
  rtime           = rtime(ms2),
  precursor_mz    = precursor_mz,
  precursor_int   = precursor_int,
  bpi             = bpi,
  tic             = tic2,
  npeaks          = npeaks,
  stringsAsFactors = FALSE
)

# Summary statistics
cat("\n=== Feature Summary ===\n")
cat("Number of MS2 spectra:", nrow(df), "\n")
cat("Precursor m/z range:", round(range(precursor_mz, na.rm = TRUE), 2), "\n")
cat("Number of peaks per spectrum:\n")
cat(" - Mean:", round(mean(npeaks), 1), "\n")
cat(" - Median:", median(npeaks), "\n")
cat(" - Range:", range(npeaks), "\n")

# Save
output_file <- "output/features/features.csv"
write.csv(df, output_file, row.names = FALSE)
cat("\n✓ Features saved to:", output_file, "\n")

# Optional: Create feature distribution plots
cat("\nGenerating feature distribution plots...\n")

png("output/qc_plots/feature_distributions.png", 
    width = 1200, height = 800, res = 120)
par(mfrow = c(2, 2))

hist(precursor_mz, breaks = 50, main = "Precursor m/z Distribution",
     xlab = "m/z", col = "lightblue", border = "black")

hist(log10(tic2), breaks = 50, main = "TIC Distribution (log10)",
     xlab = "log10(TIC)", col = "lightcoral", border = "black")

hist(npeaks, breaks = 50, main = "Peaks per Spectrum",
     xlab = "Number of Peaks", col = "lightgreen", border = "black")

plot(rtime(ms2), precursor_mz, pch = 16, cex = 0.5,
     main = "Precursor m/z vs Retention Time",
     xlab = "Retention Time (s)", ylab = "Precursor m/z")

dev.off()

cat("✓ Feature distribution plots saved\n")
