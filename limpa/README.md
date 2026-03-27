# limpa — Peptide-to-Protein Quantification & Differential Expression

## What

Aggregates peptide-level LC-MS/MS intensity data to protein level and identifies differentially expressed proteins between experimental groups using a precision-weighted `limma` linear model.

## Why

In label-free proteomics, missing peptide values are **not random** — low-abundance peptides are systematically more likely to go undetected (Missing Not At Random, MNAR). Imputing these zeros or excluding incomplete cases introduces bias into differential expression results.

`limpa` resolves this by modeling the missingness explicitly via a **Detection Probability Curve (DPC)**: a sigmoid function fitted to the relationship between peptide intensity and the probability of being observed. The DPC-derived uncertainty is then passed as **precision weights** into `limma`, so the linear model accounts for the reliability of each protein measurement. This reduces false positives compared to naive imputation strategies.

## Workflow

```
Peptide-level intensity table
    │
    ▼
newLimpa()       Initialize limpa object
addTargets()     Attach sample group metadata
    │
    ▼
estimateDPC()    Fit Detection Probability Curve per sample
                 Models P(observed | intensity) as sigmoid
    │
    ▼
DPCquant()       Summarize peptides → proteins
                 Outputs: expression matrix + precision weight matrix
    │
    ▼
limmaFit()       Fit weighted linear model (design matrix)
eBayes()         Empirical Bayes variance stabilization
contrasts.fit()  Define pairwise comparison (e.g., Treated vs Control)
    │
    ▼
topTable()       Differential expression results (logFC, adj.P.Val)
Volcano plot     Visualize significance vs. effect size
```

## Scripts

| File | Purpose |
|---|---|
| `01.Tutorial.R` | Core workflow: DPC estimation → quantification → DEA → volcano plot |
| `01.Tutorials_extended.R` | Extended examples: multi-group contrasts, additional QC |

## Key Parameters

| Parameter | Value | Rationale |
|---|---|---|
| `sig_threshold` | 0.05 | FDR (BH-adjusted p-value) cutoff |
| `fc_threshold` | 1.0 | \|log₂FC\| cutoff — 2-fold change |

## Result

- Ranked differential expression table with logFC and adjusted p-values
- Volcano plot highlighting proteins with |logFC| > 1 and FDR < 0.05
- Precision-weighted analysis reduces false positives from MNAR missingness

## Dependencies

```r
BiocManager::install(c("limpa", "limma"))
install.packages("data.table")
```
