# proteomics

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
│   ├── 01_load_spectra.R         # mzML 로딩 (Spectra)
│   ├── 02_preprocess.R           # 스펙트라 필터링
│   ├── 03_qc.R                   # TIC / MS-level QC
│   ├── 04_feature_extraction.R   # MS2 feature matrix 추출
│   ├── 05a_limpa_dea.R           # [Branch A] DPC-quant + DEA → final output
│   ├── 05b_limpa_prefilter.R     # [Branch B] DEA → responding protein filter
│   └── 06_dromics_bmd.R          # [Branch B] BMD 계산
├── config.yaml                   # 경로 · 파라미터 중앙화
├── python/run_pipeline.py        # 파이프라인 오케스트레이터
└── archive/                      # 원본 튜토리얼 스크립트 보관
```

---

## Scripts

### Common steps (01–04): Raw MS processing

**What:** `.mzML` 파일 로딩 → 불량 스펙트라 제거 → QC 시각화 → MS2 feature matrix 생성.

**Why:** `Spectra` (Bioconductor) 패키지는 mzML을 lazy-loading 방식으로 처리해 대용량 MS 데이터를 메모리 효율적으로 다룬다. 단계별 분리 구조(`01 → 04`)로 중간 단계 재실행이 용이하다.

| Script | Output |
|---|---|
| `01_load_spectra.R` | `data/processed/raw_spectra.rds` |
| `02_preprocess.R` | `data/processed/clean_spectra.rds` |
| `03_qc.R` | `results/qc/tic.png`, `ms_level_distribution.png` |
| `04_feature_extraction.R` | `data/processed/features.csv` |

---

### Branch A — `05a_limpa_dea.R`: Group comparison (final output)

**What:** 펩타이드 수준 데이터를 단백질로 정량화하고 그룹 간 차등발현단백질(DEP)을 검출한다.

**Why:** `limpa`의 DPC-quant는 결측값이 intensity에 의존하는 MNAR 구조를 모델링해 불확실성을 `limma` 가중치로 전달한다. 단순 imputation 대비 위양성 감소.

**Result:** `dea_results.csv` (logFC, adj.P.Val), `volcano.png`

---

### Branch B — `05b_limpa_prefilter.R` + `06_dromics_bmd.R`: Dose-response (final output)

**What (05b):** DEA로 dose에 반응하는 단백질을 필터링 → DROmics 입력 리스트 생성. DEA 자체가 최종 목적이 아니라 **BMD 계산을 위한 필터 도구**로 사용된다.

**What (06):** 반응 단백질에 5가지 dose-response 모델(linear, exponential, Hill, Gauss-probit, log-Gauss-probit)을 피팅하고 BMD-zSD를 계산한다.

**Why:** 단순 그룹 비교(treated vs. control)는 "변화하는가"만 답한다. BMD는 "어느 농도부터 변하기 시작하는가"를 답한다. 독성학·약물 안전성 연구에서 민감도 순위 산정에 필수적이다.

**Result:** `bmd_results.csv` (per-protein BMD + CI), `bmd_heatmap.png`

---

## Quick Start

```bash
# 1. config.yaml에서 study.mode 확인 (group 또는 dose_response)

# 2. 전체 파이프라인 실행
python python/run_pipeline.py

# 3. 특정 단계만 실행
python python/run_pipeline.py --steps 1 2 3

# 4. Branch 명시적 지정
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
