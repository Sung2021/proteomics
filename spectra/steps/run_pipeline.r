# ============================================================================
# run_pipeline.R - Master Pipeline Script
# ============================================================================

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║       Proteomics MS Data Processing Pipeline              ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# Record start time
start_time <- Sys.time()

# Function to run script with error handling
run_step <- function(script_path, step_name) {
  cat("\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("Running:", step_name, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  tryCatch({
    source(script_path)
    cat("\n✓", step_name, "completed successfully\n")
    return(TRUE)
  }, error = function(e) {
    cat("\n✗ Error in", step_name, ":\n")
    cat(e$message, "\n")
    return(FALSE)
  })
}

# Run pipeline steps
steps <- list(
  list(path = "R/01_load.R",        name = "Step 1: Load Data"),
  list(path = "R/02_preprocess.R",  name = "Step 2: Preprocess"),
  list(path = "R/03_qc.R",          name = "Step 3: Quality Control"),
  list(path = "R/04_features.R",    name = "Step 4: Feature Extraction")
)

results <- sapply(steps, function(step) {
  run_step(step$path, step$name)
})

# Summary
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "secs")

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("PIPELINE SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

for (i in seq_along(steps)) {
  status <- if(results[i]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("%-40s %s\n", steps[[i]]$name, status))
}

cat(rep("-", 60), "\n", sep = "")
cat(sprintf("Total time: %.1f seconds\n", as.numeric(elapsed)))
cat(rep("=", 60), "\n", sep = "")

if (all(results)) {
  cat("\n🎉 Pipeline completed successfully!\n")
  cat("\nOutput files:\n")
  cat(" - output/processed/raw_spectra.rds\n")
  cat(" - output/processed/clean_spectra.rds\n")
  cat(" - output/qc_plots/tic.png\n")
  cat(" - output/qc_plots/ms_level_distribution.png\n")
  cat(" - output/qc_plots/feature_distributions.png\n")
  cat(" - output/features/features.csv\n")
} else {
  cat("\n⚠️  Pipeline completed with errors. Check output above.\n")
}
