# ============================================================================
# setup.R - Initial Setup and Package Installation
# ============================================================================

cat("Starting proteomics pipeline setup...\n")

# Create directory structure
dirs <- c(
  "raw",
  "R",
  "output",
  "output/qc_plots",
  "output/processed",
  "output/features"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("Created directory:", d, "\n")
  } else {
    cat("Directory already exists:", d, "\n")
  }
}

# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Required packages
packages <- c(
  "Spectra",
  "S4Vectors",
  "msdata"  # Optional: for testing with example data
)

# Install packages
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(pkg, "is already installed.\n")
  }
}

# Verify installation
cat("\n=== Verification ===\n")
for (pkg in packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗", pkg, "failed to load\n")
  }
}

cat("\n=== Setup Complete ===\n")
cat("Place your mzML files in the 'raw/' directory.\n")
cat("Then run: source('run_pipeline.R')\n")
