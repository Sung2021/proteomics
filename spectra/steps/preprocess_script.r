# ============================================================================
# 02_preprocess.R - Preprocess spectra
# ============================================================================

cat("\n=== Step 2: Preprocessing ===\n")

library(Spectra)

# Load raw spectra
sp <- readRDS("output/processed/raw_spectra.rds")
cat("Loaded", length(sp), "raw spectra\n")

# Count before filtering
n_before <- length(sp)
ms1_before <- sum(msLevel(sp) == 1)
ms2_before <- sum(msLevel(sp) == 2)

cat("\nBefore filtering:\n")
cat(" - Total:", n_before, "\n")
cat(" - MS1:", ms1_before, "\n")
cat(" - MS2:", ms2_before, "\n")

# Apply filters
cat("\nApplying filters...\n")
sp <- sp |>
  filterEmptySpectra() |>
  filterIntensity(c(100, Inf)) |>
  filterMzRange(c(300, 2000))

# Count after filtering
n_after <- length(sp)
ms1_after <- sum(msLevel(sp) == 1)
ms2_after <- sum(msLevel(sp) == 2)

cat("\nAfter filtering:\n")
cat(" - Total:", n_after, "(", round(100*n_after/n_before, 1), "% retained)\n")
cat(" - MS1:", ms1_after, "(", round(100*ms1_after/ms1_before, 1), "% retained)\n")
cat(" - MS2:", ms2_after, "(", round(100*ms2_after/ms2_before, 1), "% retained)\n")

# Save
output_file <- "output/processed/clean_spectra.rds"
saveRDS(sp, output_file)
cat("\n✓ Saved to:", output_file, "\n")
