# ============================================================================
# 03_qc.R  [Common Step]
# Quality control plots from clean spectra
# ============================================================================

library(Spectra)
library(yaml)

cfg <- yaml::read_yaml("config.yaml")
cat("\n=== Step 3: Quality Control ===\n")

sp  <- readRDS(cfg$paths$clean_spectra)
ms1 <- filterMsLevel(sp, 1)
cat("MS1 spectra:", length(ms1), "\n")

dir.create(cfg$paths$qc_plots, recursive = TRUE, showWarnings = FALSE)

# TIC plot
tic <- spectraData(ms1)$totIonCurrent
rt  <- spectraData(ms1)$rtime

tic_out <- file.path(cfg$paths$qc_plots, "tic.png")
png(tic_out, width = 1200, height = 800, res = 120)
plot(rt, tic, type = "l",
     main = "Total Ion Current (TIC)",
     xlab = "Retention Time (s)", ylab = "Total Ion Current",
     col = "steelblue", lwd = 2)
grid()
dev.off()
cat("✓ TIC plot:", tic_out, "\n")

# MS level distribution
ms_out <- file.path(cfg$paths$qc_plots, "ms_level_distribution.png")
png(ms_out, width = 800, height = 600, res = 120)
barplot(table(msLevel(sp)),
        main = "MS Level Distribution",
        xlab = "MS Level", ylab = "Number of Spectra",
        col = c("lightblue", "lightcoral"))
dev.off()
cat("✓ MS level distribution:", ms_out, "\n")

cat("\n=== QC Metrics ===\n")
cat("TIC range  :", format(range(tic), scientific = TRUE), "\n")
cat("RT range   :", round(range(rt), 2), "seconds\n")
cat("Mean TIC   :", format(mean(tic),   scientific = TRUE), "\n")
cat("Median TIC :", format(median(tic), scientific = TRUE), "\n")
