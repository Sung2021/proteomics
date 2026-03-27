# ============================================================================
# 04_feature_extraction.R  [Common Step]
# Extract MS2 features → protein quantification matrix
# Output: data/processed/protein_matrix.csv
# ============================================================================

library(Spectra)
library(yaml)

cfg <- yaml::read_yaml("config.yaml")
cat("\n=== Step 4: Feature Extraction ===\n")

sp  <- readRDS(cfg$paths$clean_spectra)
ms2 <- filterMsLevel(sp, 2)
cat("MS2 spectra:", length(ms2), "\n")
if (length(ms2) == 0) stop("No MS2 spectra found. Check your data.")

bpi          <- sapply(intensity(ms2),  max)
tic2         <- sapply(intensity(ms2),  sum)
npeaks       <- sapply(peaksData(ms2),  nrow)
precursor_mz <- precursorMz(ms2)

feature_df <- data.frame(
  file          = dataOrigin(ms2),
  scan_number   = acquisitionNum(ms2),
  rtime         = rtime(ms2),
  precursor_mz  = precursor_mz,
  bpi           = bpi,
  tic           = tic2,
  npeaks        = npeaks,
  stringsAsFactors = FALSE
)

cat("\n=== Feature Summary ===\n")
cat("MS2 spectra          :", nrow(feature_df), "\n")
cat("Precursor m/z range  :", round(range(precursor_mz, na.rm = TRUE), 2), "\n")
cat("Peaks / spectrum     : mean =", round(mean(npeaks), 1),
    "| median =", median(npeaks), "\n")

dir.create(dirname(cfg$paths$features),    recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$qc_plots,             recursive = TRUE, showWarnings = FALSE)

write.csv(feature_df, cfg$paths$features, row.names = FALSE)
cat("\n✓ Features saved to:", cfg$paths$features, "\n")

# Distribution plots
plt_out <- file.path(cfg$paths$qc_plots, "feature_distributions.png")
png(plt_out, width = 1200, height = 800, res = 120)
par(mfrow = c(2, 2))
hist(precursor_mz, breaks = 50, main = "Precursor m/z",
     xlab = "m/z", col = "lightblue")
hist(log10(tic2), breaks = 50, main = "TIC (log10)",
     xlab = "log10(TIC)", col = "lightcoral")
hist(npeaks, breaks = 50, main = "Peaks per Spectrum",
     xlab = "# Peaks", col = "lightgreen")
plot(rtime(ms2), precursor_mz, pch = 16, cex = 0.4,
     main = "m/z vs RT", xlab = "RT (s)", ylab = "Precursor m/z")
dev.off()
cat("✓ Feature plots:", plt_out, "\n")
