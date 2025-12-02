# **R for Mass Spectrometry – Raw Data 튜토리얼 **

## **0. 목적**

Spectra 패키지를 이용한 raw MS(mzML/mzXML/netCDF/mgf) 데이터 로딩, 시각화, 필터링, 백엔드 사용법을 전체 흐름으로 학습한다.

---

# **1. Raw MS 데이터 개념**

## **1.1 Raw data 추상화**

R에서는 raw MS data를 두 가지 구성으로 본다:

* **metadata table**: spectrum-level annotation (msLevel, rtime 등)
* **spectra list**: 각 스캔의 (m/z, intensity) peak 데이터

Spectra 클래스를 사용해 이 두 요소를 통합 관리한다.

---

# **2. Spectra 클래스 소개**

## **2.1 DataFrame → Spectra 변환**

```r
library(Spectra)
library(S4Vectors)

spd <- DataFrame(
  msLevel = c(1L, 2L),
  rtime   = c(1.1, 1.2)
)
spd$mz        <- list(c(100,103.2,104.3,106.5), c(45.6,120.4,190.2))
spd$intensity <- list(c(200,400,34.2,17),       c(12.3,15.2,6.8))

sp0 <- Spectra(spd)
```

### **핵심 accessor**

* `spectraVariables()`
* `spectraData()`
* `peaksData()`
* `[` (subset)

---

# **3. mzML 파일 로딩**

```r
sp <- Spectra("file.mzML")
length(sp)    # 스캔 수
```

metadata 접근:

```r
msLevel(sp)
rtime(sp)
mz(sp[1])
intensity(sp[1])
```

---

# **4. 핵심 spectra variables**

Spectra 객체가 보장하는 코어 변수:

* **acquisitionNum**
* **centroided**
* **collisionEnergy**
* **dataOrigin**
* **dataStorage**
* **intensity**
* **isolationWindowLower/Upper/TargetMz**
* **msLevel**
* **mz**
* **polarity**
* **precScanNum**
* **precursorCharge**
* **precursorIntensity**
* **precursorMz**
* **rtime**
* **scanIndex**
* **smoothed**

추가 메타데이터는 `$` 로 추가 가능:

```r
sp$rtime_minute <- rtime(sp) / 60
```

---

# **5. 단일 spectrum·MS level 분석 실습**

## **5.1 MS level 분포 확인**

```r
table(msLevel(sp))
```

## **5.2 MS2에서 가장 높은 base peak intensity 스캔 찾기**

```r
idx <- which.max(precursorIntensity(filterMsLevel(sp, 2)))
```

## **5.3 centroided 여부 확인**

```r
centroided(sp)[1:5]
plotSpectra(sp[i])       # 시각적 확인
```

---

# **6. 예제: 특정 스캔의 MS1–MS2 관계 분석**

### **MS1:2807 스캔 및 해당 MS2 children 추출**

```r
ms_2 <- filterPrecursorScan(sp, 2807)
```

### **총 이온 전류(TIC) 플롯**

```r
with(spectraData(filterMsLevel(sp, 1)),
     plot(rtime, totIonCurrent, type="l"))
abline(v = rtime(sp)[2807], col="red")
```

### **MS1 스펙트럼에서 선택된 precursor 표시**

```r
plotSpectra(sp[2807], xlim=c(400,1000))
abline(v = precursorMz(ms_2)[-1], col="grey")
abline(v = precursorMz(ms_2)[2],  col="red")
```

### **isotopic envelope zoom-in**

```r
plotSpectra(sp[2807], xlim=c(521.2,522.5), type="l")
```

### **MS2 10개 스펙트럼 일괄 플롯**

```r
plotSpectra(ms_2[-1])
```

### **MS2 피크 labeling**

```r
mzLabel <- function(z) {
  z <- peaksData(z)[[1]]
  lbl <- format(z[,"mz"], digits=4)
  lbl[z[,"intensity"] < 1e5] <- ""
  lbl
}
plotSpectra(ms_2[7], xlim=c(126,132), labels=mzLabel)
```

---

# **7. 다중 파일 로딩 및 dataOrigin 활용**

```r
fls <- dir(system.file("sciex", package="msdata"), full.names=TRUE)
sp_sciex <- Spectra(fls)
table(dataOrigin(sp_sciex))
```

---

# **8. Spectra backend 시스템**

## **8.1 기본 backends**

1. **MsBackendMzR**

   * peaks는 필요할 때만 디스크에서 읽음(on-disk)
2. **MsBackendMemory**

   * 모든 peak 데이터를 메모리에 올림 (빠름, 메모리 많이 사용)
3. **MsBackendHdf5Peaks**

   * peak만 HDF5에 저장 (scalable)

### **backend 변경**

```r
sp_sciex <- setBackend(sp_sciex, MsBackendMemory())
sp_hdf5  <- setBackend(sp_sciex, MsBackendHdf5Peaks(), hdf5path=tempdir())
```

---

# **9. 기본 filtering**

## **스펙트럼 단위 필터**

* `filterAcquisitionNum`
* `filterDataOrigin`
* `filterEmptySpectra`
* `filterMzRange`
* `filterMsLevel`
* `filterPrecursorMzRange`
* `filterRt`
* 등등

## **peak 단위 필터**

* `filterIntensity`
* `filterMzValues`
* `combinePeaks`
* `deisotopeSpectra`
* `reduceSpectra`
* 등등

### **예제: 특정 파일의 특정 RT 구간만 추출**

```r
fls <- unique(dataOrigin(sp_sciex))
file_2 <- filterDataOrigin(sp_sciex, fls[2])
sps_sub <- filterRt(file_2, c(175,189))
```

---

# **10. Profile → Centroid 변환**

## **pickPeaks() 사용**

```r
plotSpectra(sp[2807], xlim=c(521.2,522.5))
pickPeaks(sp[2807]) |> filterIntensity(1e7) |> 
  plotSpectra(xlim=c(521.2,522.5))
```

---

# **11. 효율성**

### **백엔드별 속도·메모리 특성**

* in-memory: 빠름, 메모리 크게 사용
* on-disk: 느림, 메모리 적게 사용
* HDF5: scalable

---

# **12. Parallel Processing**

Spectra는 BiocParallel 기반 자동 병렬 처리 사용.
파일 단위로 병렬화됨 (`dataStorage` 기반).

---

# **13. Lazy Evaluation**

필터 호출 시 즉시 peak data를 수정하지 않고 **처리 큐(ProcessingQueue)** 에 축적됨.

```r
sp_sciex <- filterIntensity(sp_sciex, c(10, Inf))
sp_sciex@processingQueue     # lazy steps
reset(sp_sciex)              # 모든 lazy 처리 초기화
```

---

# **14. Interactive Visualization - SpectraVis**

* `browseSpectra(sp)` — Shiny 기반 스펙트럼 브라우저
* `plotlySpectra(sp[i])` — 인터랙티브 줌/팬

---

