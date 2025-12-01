# ==============================================================================
# limpa 패키지를 이용한 질량 분석 단백질체학 데이터 분석 워크플로우
# ------------------------------------------------------------------------------
# 이 스크립트는 펩타이드 수준 데이터에서 단백질 정량화 및 limma 기반의
# 차등 발현 분석(Differential Expression Analysis)을 수행합니다.
# ==============================================================================

# 1. 환경 설정 및 패키지 설치/로드 -----------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# limpa 및 limma 패키지 설치 (주석 해제 후 실행)
# BiocManager::install("limpa") 

# 필수 패키지 로드
library(limpa)
library(limma)
library(data.table)

# 2. 데이터 준비 및 limpa 객체 생성 -------------------------------------------

# limpa::testData()를 사용하여 예시 데이터를 로드합니다.
# 실제 분석에서는 사용자 펩타이드 수준 데이터프레임 (peptide_dt)을 로드하여 사용합니다.
# e.g., peptide_dt <- read.csv("your_peptide_data.csv")
peptide_dt <- limpa::testData()

# limpa 객체 초기화 (newLimpa)
# 데이터프레임 내의 핵심 컬럼 이름(Peptide, Protein, Sample, Intensity)을 지정합니다.
L <- limpa::newLimpa(peptide_dt, 
                     peptide_id = "Peptide", 
                     protein_id = "Protein", 
                     sample_id = "Sample", 
                     intensity = "Intensity")

# 실험 설계 정보(Targets) 추가
# testData에는 Sample 이름과 그룹 정보("Group")가 이미 포함되어 있습니다.
sample_map <- unique(peptide_dt[, c("Sample", "Group")])
L <- limpa::addTargets(L, targets = sample_map, sample_id = "Sample", condition_id = "Group")

# limpa 객체 요약 정보 출력
print("--- limpa 객체 초기화 완료 ---")
print(L)


# 3. 검출 확률 곡선 (DPC) 추정 및 시각화 --------------------------------------

# DPC는 강도에 따른 결측값 발생 확률을 모델링합니다.
L <- limpa::estimateDPC(L)
print("--- DPC 추정 완료 ---")

# DPC 시각화 (선택 사항: 모델이 데이터의 결측 패턴을 잘 포착했는지 검증)
# limpa::plotDPC(L) 


# 4. DPC-quant를 이용한 단백질 정량화 ------------------------------------------

# DPC-quant를 사용하여 펩타이드 정보를 단백질 수준으로 요약하고,
# 결측값의 불확실성을 포함한 정밀도 가중치(Precision Weights)를 계산합니다.
L <- limpa::DPCquant(L)
print("--- DPC-quant 정량화 완료 ---")

# 정량화된 단백질 발현 행렬 확인 (결측값 없이 완전한 행렬)
protein_expression <- limpa::proteinExpression(L)
print("--- 정량화된 단백질 발현 행렬 (일부) ---")
head(protein_expression)

# 정밀도 가중치 확인 (limma 분석에 사용됨)
# protein_weights <- limpa::proteinWeights(L)


# 5. limma를 사용한 차등 발현 분석 --------------------------------------------

# limpa 객체에서 limma 분석에 필요한 정보 추출
targets <- limpa::targets(L)

# 실험 설계 행렬 (Design Matrix) 생성
# 'Group' 컬럼을 기반으로 Control과 Treated 그룹을 비교하는 설계 행렬을 만듭니다.
design <- model.matrix(~0 + targets$Group)
colnames(design) <- gsub("targets\\$Group", "", colnames(design))
print("--- 실험 설계 행렬 ---")
print(design)

# limpa 객체와 설계 행렬을 사용하여 limma 모델 적합
# limpa::limmaFit는 정량화된 발현 값과 가중치를 자동으로 사용합니다.
fit <- limpa::limmaFit(L, design = design)

# 분산 추정치 안정화 (경험적 베이즈 적용)
fit <- limma::eBayes(fit)

# 비교하고자 하는 그룹 간의 대비(Contrast) 정의 및 적용
contrast_matrix <- limma::makeContrasts(Treated_vs_Control = Treated - Control, 
                                        levels = design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

# 차등 발현 결과 테이블 생성 (상위 20개 단백질을 FDR 보정)
results <- limma::topTable(fit2, coef="Treated_vs_Control", number=20, adjust.method="BH")

print("--- 상위 20개 차등 발현 단백질 결과 ---")
print(results)


# 6. 결과 시각화 (Volcano Plot) ------------------------------------

# 모든 단백질의 결과 테이블을 가져옵니다.
all_results <- limma::topTable(fit2, coef="Treated_vs_Control", number=Inf, adjust.method="BH")

# 유의미성 및 Fold Change 임계값 설정
sig_threshold <- 0.05   # 보정된 p-value (adj.P.Val) 임계값
fc_threshold <- 1       # Log2 Fold Change (|logFC|) 임계값

# 유의미한 단백질 분류
significant_proteins <- subset(all_results, 
                               all_results$adj.P.Val < sig_threshold & 
                               abs(all_results$logFC) > fc_threshold)

# 시각화를 위한 R base plot 생성
print("--- 화산 플롯 생성 시작 ---")
plot(all_results$logFC, -log10(all_results$adj.P.Val), 
     main="Volcano Plot (Treated vs Control)", 
     xlab="Log2 Fold Change", 
     ylab="-log10(Adjusted P-value)",
     pch=20,
     cex=0.6,
     col="grey50") # 기본 단백질은 회색으로 표시

# 유의미한 단백질 강조 (빨간색)
points(significant_proteins$logFC, -log10(significant_proteins$adj.P.Val), 
       col="red", pch=20, cex=0.8)

# 유의미성 경계선 추가
abline(h = -log10(sig_threshold), col="blue", lty=2)
abline(v = c(-fc_threshold, fc_threshold), col="blue", lty=2)

print(paste("유의미한 단백질 개수:", nrow(significant_proteins)))
print("--- 화산 플롯 생성 완료 ---")

# (선택 사항) 결과를 CSV 파일로 저장
# write.csv(all_results, "limpa_de_results.csv", row.names = TRUE)
