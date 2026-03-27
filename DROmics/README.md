# DROmics — Dose-Response Modeling & Benchmark Dose Calculation

## What

Fits dose-response curves to each protein's intensity profile across a dose gradient and calculates the **Benchmark Dose (BMD)** — the dose at which a biologically meaningful change first occurs. Supports pathway-level sensitivity analysis to identify which biological processes are most dose-sensitive.

## Why

Classical pairwise DEA (treated vs. control) answers *whether* a protein changes. BMD analysis answers *at what dose* it starts to change — a fundamentally different and more informative question for toxicology, pharmacology, and drug safety studies.

By fitting continuous dose-response curves and deriving a BMD per protein, this analysis enables:
- **Prioritization** of the most sensitive molecular targets
- **Risk assessment** at doses below the tested range (via model extrapolation)
- **Pathway-level integration**: grouping proteins by biological function reveals which pathways are engaged earliest along the dose axis

`DRomics` automatically selects the best-fitting model from 5 candidates (linear, exponential, Hill, Gauss-probit, log-Gauss-probit) using AIC, removing the need for manual model selection.

## Workflow

```
Protein intensity table (samples × dose)
    │
    ▼
continuousomicdata()   Import and validate data format
                       (log-scale intensities expected)
    │
    ▼
itemselect()           Screen for dose-responsive proteins
                       Method: quadratic trend test, FDR < 0.05
    │
    ▼
drcfit()               Fit 5 candidate dose-response models per protein
                       Automatic best-model selection by AIC
                       Outputs: fitted curves + trend classification
                       (inc / dec / U-shape / bell-shape)
    │
    ▼
bmdcalc()              Calculate BMD-zSD per protein
                       BMD = dose where response deviates by z×SD from baseline
    │
    ▼
bmdboot()              Bootstrap confidence intervals (n = 1000)
bmdfilter()            Retain proteins with well-defined CI
    │
    ▼
Pathway annotation     Merge with GO/pathway annotation file
sensitivityplot()      Pathway-level BMD distribution (1st quartile)
trendplot()            Dominant response trend per pathway
curvesplot()           All fitted curves colored by trend
```

## Scripts

| File | Purpose |
|---|---|
| `01.Tutorial_comprehensive.txt` | Full workflow with all parameter options documented |
| `02.Multi_levels_and_Network.txt` | Multi-omics integration and network-level analysis |

> **Note:** These files are plain-text tutorials. To run them, copy the R code blocks into an R script or RMarkdown document.

## Key Parameters

| Parameter | Value | Rationale |
|---|---|---|
| `select.method` | `"quadratic"` | Detects both monotonic and biphasic responders |
| `FDR` | 0.05 | False discovery rate for protein selection |
| `z` (BMD-zSD) | 1 | Response must deviate ≥ 1 SD from baseline |
| `niter` (bootstrap) | 1000 | For stable 95% CI estimation |

## Result

- Per-protein BMD estimates with 95% bootstrap confidence intervals
- ECDF plot of BMD distribution across all responsive proteins
- Pathway sensitivity plot: which biological processes respond at the lowest doses
- Trend classification (increasing / decreasing / U-shape / bell-shape) per protein

## Data

Uses `DRomics::proteinomics` — a yeast proteomics dose-response dataset included with the package (Larras et al. 2020, *J. Hazard. Mater.*).

## Dependencies

```r
install.packages("DRomics")
BiocManager::install(c("DESeq2", "limma"))  # required by DRomics internally
```
