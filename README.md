# proteomics-dea-bmd-pipeline

> End-to-end proteomics analysis pipeline: raw MS spectra → differential expression → dose-response modeling

---

## Pipeline Overview

```
mzML files
    │
    ├─ 01_load_spectra.R       [Common] Load mzML → Spectra object
    ├─ 02_preprocess.R         [Common] Filter empty / low-intensity spectra
    ├─ 03_qc.R                 [Common] TIC + MS-level QC plots
    ├─ 04_feature_extraction.R [Common] Extract MS2 feature matrix
    │
    ├──── study.mode: group ─────────────────────────────────┐
    │                                                         │
    │  05a_limpa_dea.R   DPC-quant + limma DEA               │
    │                    → dea_results.csv                    │
    │                    → volcano.png            [FINAL]     │
    │                                                         │
    └──── study.mode: dose_response ─────────────────────────┘
         │
         ├─ 05b_limpa_prefilter.R   DPC-quant + limma (all doses vs ctrl)
         │                          → responding_proteins.csv  [FILTER]
         │
         └─ 06_dromics_bmd.R        Dose-response fit + BMD calc
                                    → bmd_results.csv
                                    → bmd_heatmap.png          [FINAL]
```

Switch branches by editing `study.mode` in `config.yaml`, or via CLI:

```bash
python python/run_pipeline.py --mode group          # Branch A
python python/run_pipeline.py --mode dose_response  # Branch B
```

---

## Repository Structure

```
proteomics/
├── R/
│   ├── 01_load_spectra.R         # Load mzML files (Spectra)
│   ├── 02_preprocess.R           # Filter and clean spectra
│   ├── 03_qc.R                   # TIC / MS-level QC plots
│   ├── 04_feature_extraction.R   # Extract MS2 feature matrix
│   ├── 05a_limpa_dea.R           # [Branch A] DPC-quant + DEA → final output
│   ├── 05b_limpa_prefilter.R     # [Branch B] DEA → responding protein filter
│   └── 06_dromics_bmd.R          # [Branch B] BMD calculation
├── config.yaml                   # Centralized paths and parameters
├── python/run_pipeline.py        # Pipeline orchestrator
└── archive/                      # Original tutorial scripts (personal reference)
```

---

## Scripts

### Common steps (01–04): Raw MS processing

**What:** Load `.mzML` files → remove low-quality spectra → QC visualization → extract MS2 feature matrix.

**Why:** `Spectra` (Bioconductor) processes mzML via lazy-loading, making it memory-efficient for large MS datasets. The stepwise design (`01 → 04`) makes it straightforward to re-run or swap individual stages.

| Script | Output |
|---|---|
| `01_load_spectra.R` | `data/processed/raw_spectra.rds` |
| `02_preprocess.R` | `data/processed/clean_spectra.rds` |
| `03_qc.R` | `results/qc/tic.png`, `ms_level_distribution.png` |
| `04_feature_extraction.R` | `data/processed/features.csv` |

---

### Branch A — `05a_limpa_dea.R`: Group comparison (final output)

**What:** Quantify peptide-level intensities to the protein level and detect differentially expressed proteins (DEPs) between experimental groups.

**Why:** `limpa`'s DPC-quant explicitly models the Missing Not At Random (MNAR) structure — where low-abundance peptides are more likely to be missing — and propagates that uncertainty as precision weights into `limma`. This reduces false positives compared to standard imputation.

**Result:** `dea_results.csv` (logFC, adj.P.Val per protein), `volcano.png`

---

### Branch B — `05b_limpa_prefilter.R` + `06_dromics_bmd.R`: Dose-response (final output)

**What (05b):** Run DEA across all doses vs. control to identify dose-responding proteins, then export that list as input to DROmics. DEA here is a **filter tool**, not the final output.

**What (06):** Fit five dose-response models (linear, exponential, Hill, Gauss-probit, log-Gauss-probit) to each responding protein and calculate the Benchmark Dose (BMD-zSD).

**Why:** Simple group comparison (treated vs. control) only answers "does a change exist?" BMD answers "at what dose does the change begin?" This is essential for prioritizing the most sensitive molecular targets in toxicology and drug safety studies.

**Result:** `bmd_results.csv` (per-protein BMD + 95% CI), `bmd_heatmap.png`

---

## Quick Start

```bash
# 1. Set study.mode in config.yaml ("group" or "dose_response")

# 2. Run the full pipeline
python python/run_pipeline.py

# 3. Run specific steps only
python python/run_pipeline.py --steps 1 2 3

# 4. Override branch explicitly
python python/run_pipeline.py --mode dose_response
```

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("Spectra", "MsBackendMzR", "limpa", "limma", "DESeq2"))

# CRAN
install.packages(c("DRomics", "yaml", "data.table", "ggplot2"))
```

```python
pip install pyyaml
```
