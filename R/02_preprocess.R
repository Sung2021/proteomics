# ============================================================================
# 02_preprocess.R  [Common Step]
# Filter and clean raw spectra
# ============================================================================

library(Spectra)
library(yaml)

cfg <- yaml::read_yaml("config.yaml")
cat("\n=== Step 2: Preprocessing ===\n")

sp       <- readRDS(cfg$paths$raw_spectra)
n_before <- length(sp)
cat("Loaded", n_before, "raw spectra\n")

sp <- sp |>
  filterEmptySpectra() |>
  filterIntensity(c(cfg$preprocessing$min_intensity, Inf)) |>
  filterMzRange(cfg$preprocessing$mz_range)

n_after <- length(sp)
cat(sprintf("\nRetained %d / %d spectra (%.1f%%)\n",
            n_after, n_before, 100 * n_after / n_before))
cat("  MS1:", sum(msLevel(sp) == 1), "\n")
cat("  MS2:", sum(msLevel(sp) == 2), "\n")

dir.create(dirname(cfg$paths$clean_spectra), recursive = TRUE, showWarnings = FALSE)
saveRDS(sp, cfg$paths$clean_spectra)
cat("\n✓ Saved to:", cfg$paths$clean_spectra, "\n")
