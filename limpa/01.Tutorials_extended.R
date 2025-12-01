# ==============================================================================
# limpa 패키지를 이용한 질량 분석 단백질체학 데이터 분석 워크플로우 (확장 버전)
# ------------------------------------------------------------------------------
# 펩타이드 수준 데이터에서 DPC 기반 정량화, limma 차등 발현, 심층 시각화를 수행합니다.
# 시각화는 ggplot2를 사용합니다.
# ==============================================================================

# 1. 환경 설정 및 패키지 로드 -----------------------------------------------

# limpa는 Bioconductor 패키지입니다.
# BiocManager::install(c("limpa", "ComplexHeatmap", "circlize")) 

library(limpa)
library(limma)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap) # 히트맵을 위한 고급 패키지
library(circlize)       # 히트맵 색상 팔레트 관리를 위한 패키지


# 2. 데이터 로드 및 limpa 객체 생성 -------------------------------------------

# limpa::testData()를 사용하여 예시 데이터를 로드합니다.
peptide_dt <- limpa::testData()

# limpa 객체 초기화
L <- limpa::newLimpa(peptide_dt, 
                     peptide_id = "Peptide", 
                     protein_id = "Protein", 
                     sample_id = "Sample", 
                     intensity = "Intensity")

# 실험 설계 정보(Targets) 추가
sample_map <- unique(peptide_dt[, c("Sample", "Group")])
L <- limpa::addTargets(L, targets = sample_map, sample_id = "Sample", condition_id = "Group")

print("--- limpa 객체 초기화 완료 ---")
print(L)


# 3. 데이터 탐색 (EDA): 결측값 패턴 시각화 ----------------------------------

# limpa 객체에서 원시 데이터 행렬을 추출하여 결측값 패턴을 확인합니다.
raw_expression <- limpa::getExpression(L)

# 결측값 (NA) 정보를 추출하여 시각화 데이터 준비
na_data <- raw_expression %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Intensity") %>%
  mutate(is_missing = is.na(Intensity)) %>%
  left_join(limpa::targets(L), by = "Sample")

# 샘플별 결측값 비율 시각화
missing_value_plot <- na_data %>%
  group_by(Sample, Group) %>%
  summarise(Missing_Prop = mean(is_missing), .groups = 'drop') %>%
  ggplot(aes(x = reorder(Sample, Missing_Prop), y = Missing_Prop, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample-wise Missing Value Proportion (Raw Data)",
       x = "Sample",
       y = "Proportion of Missing Protein Intensities")

print(missing_value_plot)
# 


# 4. 검출 확률 곡선 (DPC) 추정 및 시각화 --------------------------------------

# DPC 추정
L <- limpa::estimateDPC(L)
print("--- DPC 추정 완료 ---")

# DPC 시각화: 강도와 검출 확률의 관계를 확인
print("--- DPC 플롯 생성 (새 창에 출력) ---")
limpa::plotDPC(L)
# 


# 5. DPC-quant를 이용한 단백질 정량화 및 정밀도 가중치 계산 ------------------------

# DPC를 기반으로 정량화 및 정밀도 가중치를 계산합니다.
L <- limpa::DPCquant(L)
print("--- DPC-quant 정량화 완료 ---")

# 정량화된 단백질 발현 행렬 확인
protein_expression <- limpa::proteinExpression(L)


# 6. 품질 관리 (QC) 분석: PCA/MDS ------------------------------------------

# 정량화된 단백질 발현 행렬을 사용하여 PCA를 수행하여 샘플 간의 군집화를 확인합니다.
pca_result <- prcomp(t(protein_expression), scale. = TRUE)
pca_data <- as.data.frame(pca_result$x) %>%
  rownames_to_column("Sample") %>%
  left_join(limpa::targets(L), by = "Sample")

# 분산 설명력 계산
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ggplot2 PCA 플롯 생성
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Quantified Protein Expression (QC)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme(legend.position = "bottom")

print(pca_plot)
# 


# 7. limma를 사용한 차등 발현 분석 --------------------------------------------

targets <- limpa::targets(L)
design <- model.matrix(~0 + targets$Group)
colnames(design) <- gsub("targets\\$Group", "", colnames(design))

# limma 모델 적합 (가중치 자동 적용)
fit <- limpa::limmaFit(L, design = design)
fit <- limma::eBayes(fit)

# 대비 정의 및 적용
contrast_matrix <- limma::makeContrasts(Treated_vs_Control = Treated - Control, levels = design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

# 차등 발현 결과 테이블 생성
all_results <- limma::topTable(fit2, coef="Treated_vs_Control", number=Inf, adjust.method="BH")
print("--- 상위 10개 차등 발현 단백질 결과 ---")
print(head(all_results, 10))


# 8. ggplot2를 이용한 결과 시각화: 화산 플롯 (Volcano Plot) --------------------

# 유의미성 및 Fold Change 임계값 설정
sig_threshold <- 0.05   
fc_threshold <- 1       

plot_data <- all_results %>%
  mutate(
    neg_log10_P = -log10(adj.P.Val),
    # 유의미한 단백질 그룹 분류
    Significance = case_when(
      adj.P.Val < sig_threshold & logFC > fc_threshold ~ "Up-regulated",
      adj.P.Val < sig_threshold & logFC < -fc_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot <- ggplot(plot_data, aes(x = logFC, y = neg_log10_P, color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = "blue", 
                                "Up-regulated" = "red", 
                                "Not significant" = "gray")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Treated vs Control Differential Protein Expression",
    x = expression(Log[2]~"Fold Change"),
    y = expression("-Log"[10]~"(Adjusted P-value)"),
    caption = paste("Adj P-value <", sig_threshold, "and |Log2FC| >", fc_threshold)
  ) +
  theme(legend.position = "bottom")

print(volcano_plot)
# 


# 9. 추가 분석: 단백질 가중치 분포 시각화 및 히트맵 ------------------------------------

# 9.1. 단백질 가중치 분포 히스토그램
weights_dt <- data.frame(Weight = limpa::proteinWeights(L)[,1])

weights_plot <- ggplot(weights_dt, aes(x = Weight)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Protein Precision Weights (DPC-quant)",
       x = "Precision Weight",
       y = "Count of Proteins")

print(weights_plot)
# 


# 9.2. 상위 차등 발현 단백질 (DEP) 히트맵 시각화
top_n_proteins <- 50 # 상위 50개 단백질 선택
top_dep <- all_results %>%
  top_n(top_n_proteins, wt = -adj.P.Val) %>% # 가장 낮은 p-value를 가진 단백질 선택
  rownames()

# 히트맵을 위한 데이터 준비: Z-score 표준화
heatmap_data <- protein_expression[top_dep, ]
heatmap_scaled <- t(scale(t(heatmap_data))) # 행(단백질)별 Z-score 스케일링

# 샘플 정보(Group)를 주석으로 추가
col_ha <- HeatmapAnnotation(Group = targets$Group,
                            col = list(Group = c("Control" = "orange", "Treated" = "purple")),
                            annotation_name_side = "left")

# 히트맵 그리기
print("--- 상위 차등 발현 단백질 히트맵 생성 (새 창에 출력) ---")
heatmap_plot <- Heatmap(heatmap_scaled, 
                        name = "Z-Score",
                        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                        top_annotation = col_ha,
                        show_row_names = FALSE,
                        column_title = paste("Top", top_n_proteins, "DE Proteins Heatmap"),
                        row_names_gp = gpar(fontsize = 6),
                        column_names_gp = gpar(fontsize = 8))

# print(heatmap_plot) # ComplexHeatmap 객체 출력



print(paste("총 유의미하게 발현된 단백질 개수:", sum(plot_data$Significance != "Not significant")))
# write.csv(all_results, "limpa_de_results.csv", row.names = TRUE)
