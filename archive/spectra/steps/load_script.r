# ============================================================================
# 01_load.R - Load mzML files
# ============================================================================

cat("\n=== Step 1: Loading Data ===\n")

library(Spectra)

# Find all mzML files
files <- list.files("raw", pattern = "mzML$", full.names = TRUE)

if (length(files) == 0) {
  stop("No mzML files found in 'raw/' directory. Please add data files.")
}

cat("Found", length(files), "mzML files:\n")
for (f in files) {
  cat(" -", basename(f), "\n")
}

# Load spectra
cat("\nLoading spectra...\n")
sp <- Spectra(files)

# Summary
cat("\n=== Summary ===\n")
cat("Total spectra:", length(sp), "\n")
cat("MS levels:", toString(unique(msLevel(sp))), "\n")
cat("Retention time range:", round(range(rtime(sp)), 2), "seconds\n")

# Save
output_file <- "output/processed/raw_spectra.rds"
saveRDS(sp, output_file)
cat("\n✓ Saved to:", output_file, "\n")
