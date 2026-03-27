# proteomics

> End-to-end proteomics analysis: raw MS spectra → differential expression → dose-response modeling

---

## Pipeline Overview

```
mzML files
    │
    ▼
[spectra]  Load → Preprocess → QC → Feature Extraction
    │
    ▼
[limpa]    Peptide → Protein Quantification → Differential Expression (limma)
    │
    ▼
[DROmics]  Dose-Response Modeling → BMD Calculation → Pathway Sensitivity
```

---

## Modules

### 1. `spectra/` — Raw MS Data Processing

**What:** Ingests raw `.mzML` mass spectrometry files and extracts a clean feature matrix ready for downstream statistical analysis.

**Why:** `MsBackendMzR` / `Spectra` (Bioconductor) provides a memory-efficient, lazy-loading interface to mzML — essential for large MS datasets where loading everything into RAM is impractical. The modular step design (load → preprocess → QC → features) makes it easy to swap or re-run individual stages.

**Key Steps:**
| Script | Purpose |
|---|---|
| `steps/load_script.r` | Read mzML → create `Spectra` object |
| `steps/preprocess_script.r` | Filter empty spectra, remove low-intensity signals (< 100), restrict m/z 300–2000 |
| `steps/qc_script.r` | TIC distribution, peak count QC plots |
| `steps/features_script.r` | Bin spectra, extract intensity feature matrix |
| `steps/run_pipeline.r` | Master script: runs all steps in order |

**Result:** A sample × feature intensity matrix (`output/features/`) and QC plots (`output/qc_plots/`) confirming data quality before statistical testing.

---

### 2. `limpa/` — Peptide-to-Protein Quantification & Differential Expression

**What:** Aggregates peptide-level intensity measurements to protein level and identifies differentially expressed proteins between experimental groups using a precision-weighted linear model.

**Why:** Standard protein summarization (e.g., mean/median) ignores the fact that missing peptide values are not random — they are more likely to occur for low-abundance peptides (Missing Not At Random, MNAR). `limpa` explicitly models this via a **Detection Probability Curve (DPC)**, then propagates the resulting uncertainty as precision weights into `limma`'s linear model. This reduces false positives caused by imputation artifacts.

**Key Steps:**
1. Build `limpa` object from peptide-level data
2. Estimate DPC — models missingness as a function of intensity
3. `DPCquant()` — summarize peptides → proteins with precision weights
4. `limmaFit()` + `eBayes()` — differential expression with empirical Bayes variance stabilization
5. Volcano plot — visualize logFC vs. −log₁₀(adj. p-value)

**Result:** Ranked differential expression table (logFC, adj.P.Val) and a volcano plot highlighting significantly changed proteins (|logFC| > 1, FDR < 0.05).

---

### 3. `DROmics/` — Dose-Response Modeling & BMD Calculation

**What:** Fits dose-response curves to each protein's intensity profile across a dose gradient and calculates the **Benchmark Dose (BMD)** — the dose at which a biologically meaningful response begins.

**Why:** Classical pairwise DEA (treated vs. control) only tests whether a change exists at a single dose. BMD analysis identifies *at what dose* each protein starts to change, enabling prioritization of the most sensitive molecular targets. This is critical in toxicology and drug safety studies where dose-response shape matters as much as the presence of an effect.

**Key Steps:**
1. Import continuous proteomics data (`continuousomicdata()`)
2. Select dose-responsive proteins (`itemselect()`, FDR < 0.05, quadratic trend test)
3. Fit 5 candidate models (linear, exponential, Hill, Gauss-probit, log-Gauss-probit) — best model selected by AIC
4. Calculate BMD-zSD with 95% bootstrap CI (`bmdboot()`, n = 1000)
5. Pathway-level sensitivity analysis (`sensitivityplot()`, `trendplot()`)

**Result:** Per-protein BMD estimates with confidence intervals; pathway-level ECDF plots revealing which biological processes are most sensitive to the dose gradient.

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("Spectra", "MsBackendMzR", "limpa", "limma", "DESeq2"))

# CRAN
install.packages(c("DRomics", "data.table", "ggplot2"))
```

## Data

All modules use **built-in example datasets** from the respective packages:
- `spectra`: synthetic `.mzML` files generated via `MsBackendMzR`
- `limpa`: `limpa::testData()` — peptide-level MS intensity data
- `DROmics`: `DRomics::proteinomics` — yeast proteomics dose-response data (included in package)

To use your own data, replace the data loading step in each module's script with your file path.
