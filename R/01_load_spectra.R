# ============================================================================
# 01_load_spectra.R  [Common Step]
# Load mzML files using the Spectra package
# ============================================================================

library(Spectra)
library(yaml)

cfg  <- yaml::read_yaml("config.yaml")
cat("\n=== Step 1: Loading mzML files ===\n")

files <- list.files(cfg$paths$raw_mzml, pattern = "\\.mzML$", full.names = TRUE)
if (length(files) == 0) stop("No mzML files found in: ", cfg$paths$raw_mzml)

cat("Found", length(files), "mzML file(s):\n")
for (f in files) cat(" -", basename(f), "\n")

sp <- Spectra(files)

cat("\n=== Summary ===\n")
cat("Total spectra   :", length(sp), "\n")
cat("MS levels       :", toString(unique(msLevel(sp))), "\n")
cat("RT range        :", round(range(rtime(sp)), 2), "seconds\n")

dir.create(dirname(cfg$paths$raw_spectra), recursive = TRUE, showWarnings = FALSE)
saveRDS(sp, cfg$paths$raw_spectra)
cat("\n✓ Saved to:", cfg$paths$raw_spectra, "\n")
