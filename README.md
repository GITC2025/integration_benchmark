# Integration benchmark scRNA seq analysis

_figuring out the Frontenac HPC, babysitting the downloads, and benchmarking integration • HSCI591 W2026 BHSC, Queen's Uni_

* a first foray into HPC scRNA analysis 
* start date: mid March 2026 after lit review complete (Jan to March 2026) 
* IDE: Positron • Environment: Frontenac Cluster HPC (up to May), then Alliance Canada Narval (May onwards)
* [**wiki - concepts, protocols, etc.**](https://github.com/GITC2025/integration_benchmark/wiki)
* [**custom sctk 2.18.0 + R 4.5.3 container**](https://github.com/GITC2025/integration_benchmark/blob/main/custom_sctk_container.md)

# Workflow
<img width="1000" height="550" alt="Screenshot 2026-04-19 115316" src="https://github.com/user-attachments/assets/373dba5a-1ea8-4fff-a82c-ea01e5a84311" />
<br>when building the containerised pipelines: the pink ovals are logic gates to continue or pause the pipeline as appropriate.

## Methodology
```r
Part One: Data acquisition, alignment and QC
├──1. download fastq.gz from ENA database + md5 checksum
├──2. FastQC
├──3. Alignment: zcat check chemistry, Cellranger dataset 1 + Celescope dataset 2 (2 runs)
	├── Cellranger QC check (done) + celescope QC check (done)
	├── dataset 2 requires background removal (done)
├──4. QC with sctk apptainer or traditional multi-step QC
	├── build custom campbio singlecellTK container (done)
	├── Seurat V5 QC Dataset 1 (done)
	├── Seurat V5 QC Dataset 2 (ongoing)

Bridge
├──5. Normalization and Feature selection (HVGs): SCTransform (predictive modeling method) (ongoing)

Part Two: Core Analysis
├──6. Dimension reduction
	├── PCA: scree plot to check PCs cover 70-90% of variance
	├── choose PCs for correction (cover 90-95% variance)

├──7. Integration benchmarking
	├── control vs scANVI vs Harmony etc. benchmark across diff categories of tools
	├── evaluate latent space with scIB

├──8. Clustering 
	├── Leiden clustering: adjust resolution 
	├── UMAP as qualitative benchmark
		├── set seeds, try diff k values for stability
	├── overlay Leiden cluster labels on UMAP
	├── ID outliers for final QC

├──9. Findmarkers Seurat to ID celltypes via gene enrichment (unsupervised)
	├── reference with database
```
# Dataset 1
Lai, H., Cheng, X., Liu, Q., Luo, W., Liu, M., Zhang, M., Miao, J., Ji, Z., Lin, G. N., Song, W., Zhang, L., Bo, J., Yang, G., Wang, J., & Gao, W. Q. (2021). Single-cell RNA sequencing reveals the epithelial cell heterogeneity and invasive subpopulation in human bladder cancer. International journal of cancer, 149(12), 2099–2115. https://doi.org/10.1002/ijc.33794

* https://www.ebi.ac.uk/ena/browser/view/SRP217277 
* 7 primary tumour samples + 1 paracancerous ('normal') sample 

# Dataset 2
Luo, Y., Tao, T., Tao, R., Huang, G., & Wu, S. (2022). Single-Cell Transcriptome Comparison of Bladder Cancer Reveals Its Ecosystem. Frontiers in oncology, 12, 818147. https://doi.org/10.3389/fonc.2022.818147
* https://www.ebi.ac.uk/ena/browser/view/SRP351272
* 4 samples total: 2x pri BCa, 1 recurrent BCa and one cystitis glandularis

| Run | disease_state | tissue |
| :--- | :--- | :--- |
| SRR17259462 | primary bladder cancer | bladder cancer |
| SRR17259463 | recurrent bladder cancer | bladder cancer |
| SRR17259464 | primary bladder cancer | bladder cancer |
| SRR17259465 | cystitis glandularis | cystitis glandularis |

## [1. download, md5 checksum scripts](https://github.com/GITC2025/integration_benchmark/blob/main/1_download_script.md) 
## [2. FastQC and zcat chemistry checks](https://github.com/GITC2025/integration_benchmark/blob/main/2_fastqc_zcat.md) 
## [3a. Cellranger alignment Dataset 1](https://github.com/GITC2025/integration_benchmark/blob/main/3a_alignment_dataset1.md)
## [3b. Celescope alignment Dataset 2](https://github.com/GITC2025/integration_benchmark/blob/main/3b_alignment_dataset2.md)

***
# QC Dataset 1
## build containers with Campbio singlecellTK and R
get the docker image with apptainer pull
```bash
apptainer pull /global/scratch/$USER/sctk_qc_latest.sif docker://campbio/sctk_qc:latest
```
```bash
# check versions in the latest container
apptainer exec sctk_qc_latest.sif R -e "packageVersion('singleCellTK')"

R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
...
> packageVersion('singleCellTK')
[1] ‘2.14.0’
```
* this is outdated, so we build a custom container with latest R + sctk [here](https://github.com/GITC2025/integration_benchmark/blob/main/custom_sctk_container.md)
* result: a custom container called sctk_v2.18.0_R_v4.5.3.sif
* after much shimming and debugging, our custom container is 95% functional, but needs more refinement
* we will now try out traditional QC with Seurat V5, doublet removal etc, or the Teichman Lab SCTK
* out custom container has all the needed R tools so we'll use it
* [container-building saga here](https://github.com/GITC2025/integration_benchmark/blob/main/stck_container_test.md)

## QC with Seurat V5

* https://satijalab.org/seurat/articles/get_started_v5_new
* QC metrics write up: UMI threshold, gene count, mito (and ribo) %
* working node: 8CPU, 128GB RAM
* there are various strategies and resolutions of QC, we will QC the dataset as a whole, use a mainstream approach here and explore alternative methods
  
```bash
apptainer exec sctk_v2.18.0_R_v4.5.3.sif Rscript -e "write.table(installed.packages()[c('Seurat', 'SeuratObject'), c('Package', 'Version')], col.names=F, row.names=F, quote=F)"

Seurat 5.4.0
SeuratObject 5.4.0
```
enter the apptainer to use R (use batch scripts if you have a lot of datasets)
```bash
apptainer exec \
--bind /global/scratch/$USER:/global/scratch/$USER \
/global/scratch/$USER/sctk_v2.18.0_R_v4.5.3.sif \
R
```
verify our working folder is visible
```R
> list.files("/global/scratch/username/cellranger_outs_dataset1")

[1] "SRR12539462" "SRR12539463" "SRR14615558" "SRR9897621"  "SRR9897622"
[6] "SRR9897623"  "SRR9897624"  "SRR9897625"
```
process all 8 SRRs at once 
```R
library(dplyr)
library(Seurat)
library(patchwork)
library(parallel)

base_dir <- "/global/scratch/username/cellranger_outs_dataset1"
srr_ids <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

process_srr <- function(srr_id) {
matrix_path <- file.path(base_dir, srr_id, "filtered_feature_bc_matrix")
counts <- Read10X(data.dir = matrix_path)
obj <- CreateSeuratObject(counts = counts, project = srr_id, min.cells = 3, min.features = 200)
return(obj)
}

dataset1_list <- mclapply(srr_ids, process_srr, mc.cores = 8)
names(dataset1_list) <- srr_ids
```
verify 8 objects
```r
library(stringr)
dataset1_list <- dataset1_list[str_sort(names(dataset1_list), numeric = TRUE)]
length(dataset1_list)
names(dataset1_list)
dataset1_list
saveRDS(dataset1_list, file = "/global/scratch/username/dataset1_list_raw.rds")
```
good match with our cellranger QC reports for estimated cell counts
```r
# output
[1] 8
[1] "SRR9897621"  "SRR9897622"  "SRR9897623"  "SRR9897624"  "SRR9897625"
[6] "SRR12539462" "SRR12539463" "SRR14615558"
$SRR9897621
An object of class Seurat
26480 features across 8281 samples within 1 assay
Active assay: RNA (26480 features, 0 variable features)
 1 layer present: counts

$SRR9897622
An object of class Seurat
26869 features across 11423 samples within 1 assay
Active assay: RNA (26869 features, 0 variable features)
 1 layer present: counts

$SRR9897623
An object of class Seurat
25640 features across 27978 samples within 1 assay
Active assay: RNA (25640 features, 0 variable features)
 1 layer present: counts

$SRR9897624
An object of class Seurat
26153 features across 37537 samples within 1 assay
Active assay: RNA (26153 features, 0 variable features)
 1 layer present: counts

$SRR9897625
An object of class Seurat
25506 features across 10554 samples within 1 assay
Active assay: RNA (25506 features, 0 variable features)
 1 layer present: counts

$SRR12539462
An object of class Seurat
24722 features across 7846 samples within 1 assay
Active assay: RNA (24722 features, 0 variable features)
 1 layer present: counts

$SRR12539463
An object of class Seurat
25989 features across 19821 samples within 1 assay
Active assay: RNA (25989 features, 0 variable features)
 1 layer present: counts

$SRR14615558
An object of class Seurat
26474 features across 17190 samples within 1 assay
Active assay: RNA (26474 features, 0 variable features)
 1 layer present: counts
 ```
show QC metrics for first 5 cells of first SRR object
```R
head(dataset1_list[[1]]@meta.data, 5)
```
```txt
                   orig.ident nCount_RNA nFeature_RNA
AAACCTGAGAAGAAGC-1 SRR9897621      34803         6673
AAACCTGAGACCGGAT-1 SRR9897621       1628          899
AAACCTGAGCTACCTA-1 SRR9897621      14588         3728
AAACCTGAGGTTACCT-1 SRR9897621       3554         1601
AAACCTGCAAGACGTG-1 SRR9897621       6109         2949
```
* calculate mito and ribo percentages
* in general, high mito + low ribo = dying/stressed cells
* you can also use qs library for parallel save into a qs object
```bash
library(parallel)
library(Seurat)

dataset1_list <- mclapply(dataset1_list, function(x) {
x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
x <- PercentageFeatureSet(x, pattern = "^RP[SL]", col.name = "percent.ribo")
return(x)
}, mc.cores = 8)

head(dataset1_list[[1]]@meta.data, 5)

save_path <- "/global/scratch/hpc6140/R_output_dataset1/dataset1_allQCmetrics.rds"
saveRDS(dataset1_list, file = save_path)

cat("saved to ", save_path, "\n")
```
```txt
# Verify the new columns exist
head(dataset1_list[[1]]@meta.data, 5)
                   orig.ident nCount_RNA nFeature_RNA percent.mt percent.ribo
AAACCTGAGAAGAAGC-1 SRR9897621      34803         6673  4.7610838    19.127661
AAACCTGAGACCGGAT-1 SRR9897621       1628          899  2.4570025    19.656020
AAACCTGAGCTACCTA-1 SRR9897621      14588         3728  2.5911708    20.496298
AAACCTGAGGTTACCT-1 SRR9897621       3554         1601  0.7597074     7.765898
AAACCTGCAAGACGTG-1 SRR9897621       6109         2949  2.2589622    15.027009
```
## visualise QC metrics per sample and globally 
### [dataset 1 metrics and plots](https://github.com/GITC2025/integration_benchmark/blob/main/dataset1_QC_code_metrics.md)
* we have a very healthy dataset 
* ribo percentage threholds use data distribution to remove outliers, rather than strict cutoffs
<img width="1200" height="600" alt="dataset1_QC_by_SRR" src="https://github.com/user-attachments/assets/e34df537-e6e3-4695-b044-e2c91fb4ea7b" />
<br>
<img width="600" height="400" alt="dataset1_QC_Merged_Global" src="https://github.com/user-attachments/assets/625fb3d9-4531-4685-bfc2-61b89b642853" />

### Analysis Seurat V5 QC violin plots
- SRR14615558 is the only paracancerous ('normal') sample here (see metadata)
- long upper tails in feature and counts indicate presence of multiplets, to be filtered
  - ranges from 5 to 11% in our partial output from our sctk run scDblFinder
- density in lower ranges feature and count indicate empty droplets/low quality cells to be filtered
- inter-sample variability in feature/counts/ribo% could indicate batch fx and/or bio variations (need downstream integration)
- mito % is in a fairly tight range, indicating high-quality dataset
  - mito% should correlate inversely with ribo% (visible in violin plots, to confirm via scatter plot)
  - stressed/dying cells hi mito% and low ribo%, and vice versa for high quality cells

### rank and quantiles
<img width="700" height="700" alt="dataset1_Global_Knee_Plots_Sandwich_" src="https://github.com/user-attachments/assets/0f86b81f-731b-478a-9537-ee4bfedd5757" /><br>
* this confirms that standard hard floors and ceilings will safely remove a small subset of outliers
* high resolution quantile info in QC metrics summaries
* calculating our ecdf for our chosen floors/ceilings, we get this
* ribo% - in this bladder cancer tumour cells context: the majority of cells range btw 15% to 30% - characteristic S shape
```r
floors
count < 500: 0.0476%
feat < 200:  0.0014%

ceilings
count > 40k: 0.534%
feat > 8k:  0.1643%
mt > 20%:    1.8374%

dataset attrition estimate
total estimated loss: ~2.42%
```
### [QC metrics summaries](https://github.com/GITC2025/integration_benchmark/blob/main/dataset1_QC_code_metrics.md#per-sample-qc-metrics)

## plot pairwise relationships for QC metrics

<img width="700" height="700" alt="dataset1_QC_Comprehensive_Scatters" src="https://github.com/user-attachments/assets/7aed9c50-376b-47c3-b76e-c02dd73914f9" />

### Analysis Seurat V5 QC pairwise scatter
- feature vs RNA count, R = 0.93
  - expected strong linear r/s - successful library capture
  - upper right quadrant: multiplets
  - lower left quadrant: empty droplets/low Q cells 
- mt % vs RNA count, R = -0.12
  - expected L shape correlation
  - tall left side - stressed/dying cells w/ high mt % and low RNA count
  - dynamic thresholds + higher upper mito limit accounts for metabolic context of cancer cells 
- mt % vs ribo %, R= -0.32
  - inverse correlation expected
  - top left quadrant: high mt, low ribo to filter out
- ribo % vs RNA count, R = 0.15
  - bio variance and/or batch fx across SRRs 
  - moderate correlation of high mito with low ribo - the mito ceiling would filter out low ribo cells
  - high ribo% in cancer cells could indicate high proliferation, this is of biological interest

### consider metadata characteristics
- consider age and staging as confounding covariates
- Metadata from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA558456

| Run | AGE | tumor_stage | subject_status | sex |
| :--- | :--- | :--- | :--- | :--- |
| SRR14615558 | 80 | pT2 | paracancerous tissue ('normal') | male |
| SRR12539462 | 58 | pTa | patient with bladder cancer | male |
| SRR12539463 | 66 | pTa | patient with bladder cancer | male |
| SRR9897621 | 67 | Ta | patient with bladder cancer | male |
| SRR9897622 | 67 | T1 | patient with bladder cancer | male |
| SRR9897623 | 38 | T1 | patient with bladder cancer | male |
| SRR9897624 | 80 | T2 | patient with bladder cancer | male |
| SRR9897625 | 81 | T3 | patient with bladder cancer | male |
***
## Hard filter + scDblFinder before cluster QC
* based on our [quantile summary](https://github.com/GITC2025/integration_benchmark/blob/main/dataset1_QC_code_metrics.md#high-resolution-quantiles-to-justify-standard-floorsceilings) we can safely set standard hard filters
* ignore ribo % floors/ceilings for now as it leads to over-filtering and there are no established threholds - ribo% is of bio interest
* **to remove doublets effectively, we only set mito 20% ceiling, feature and count floors, and retain the feature and count ceilings**
* this allows scDblFinder to learn effectively from normal cells, without being distorted by low feature/count cells, and work on the feature/count ceiling space
* any high feature/count cells missed by scDblFinder can then be passed through the hard ceiling filter again
* see [scDblfinder 1.24.10 vignette ](https://www.bioconductor.org/packages//release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html#should-i-run-qc-cell-filtering-before-or-after-doublet-detection) for this workflow
* since we're using the filtered matrix which is the output from the cellranger emptydroplets algorithm, we can skip emptydrops
* use scDblFinder based on our prelim results from our [partial sctk run](https://github.com/GITC2025/integration_benchmark/blob/main/stck_container_test.md#partial-results-from-custom-sctk-container-on-dataset-1-4-out-of-8) and existing literature
* from the prelim results, background contamination was extremely low, skip DecontX to avoid over-correction, there's additional clustering QC downstream 
```r
floors
count < 500: 0.0476%
feat < 200:  0.0014%

ceilings
count > 40k: 0.534%
feat > 8k:  0.1643%
mt > 20%:    1.8374%

dataset attrition estimate
total estimated loss: ~2.42%
```
### pre filtered subset for scDblfinder
```r
count_floor <- 500
feat_floor <- 200
mt_ceil <- 20

dataset1_pre_dblfilter_list <- lapply(dataset1_list, function(obj) {
subset(obj, subset = nCount_RNA >= count_floor & 
nFeature_RNA >= feat_floor & 
percent.mt < mt_ceil)
```
## run scDblfinder
we use default parameters for random approach
```
pre dblfinder filter summary
      sample raw_cells pre_dblfilter_cells removed removed_pct
  SRR9897621      8281                8169     112        1.35
  SRR9897622     11423               11216     207        1.81
  SRR9897623     27978               27646     332        1.19
  SRR9897624     37537               37378     159        0.42
  SRR9897625     10554               10368     186        1.76
 SRR12539462      7846                7697     149        1.90
 SRR12539463     19821               18861     960        4.84
 SRR14615558     17190               16701     489        2.84
```
verify that our filters worked
```
orig.ident   n_cells min_count median_count min_feat median_feat max_mt mean_mt mean_ribo
SRR9897621      8169       500       7327.0      254        2636  19.98    2.68     16.45
SRR9897622     11216       500       5788.5      217        2536  19.95    5.14     20.72
SRR9897623     27646       500       1267.0      274         747  19.98    3.49     24.05
SRR9897624     37378       500       1333.0      258         767  19.96    6.58     21.36
SRR9897625     10368       501       8734.0      242        2743  19.89    4.05     17.14
SRR12539462     7697       500      12658.0      202        3136  19.99    3.56     37.21
SRR12539463    18861       500       2867.0      260        1284  19.98    6.64     24.17
SRR14615558    16701       500       2282.0      200        1041  19.97    4.66     25.64
```
<img width="700" height="460" alt="dataset1_pre_dblfilter_QC_violins" src="https://github.com/user-attachments/assets/2c06af46-9c2c-4e7d-b63d-5fa31804a4a1" /><br>

### analysis of pre-Dblfinder filtered subset
- ncountRNA is log transformed
- significant bio/batch variance across samples 
- very strong correlation btw nCount and nFeature violin shapes
- mito % range is narrow - good quality dataset
- some low ribo% cells filtered out with the mito% ceiling
- ribo% variance remains
  - SRR...462 shows ribo% upregulation - could indicate ↑ proliferation
  - SRR...558 (paracancerous 'normal') shows a nice distribution up to 40%

## run scDblFinder per sample 
* https://www.bioconductor.org/packages//release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
* we use default parameters for random approach
* histogram diagnostic per sample (comparison is not exact across samples with diff cell counts)
* follow default parameters of 0.8% doublet rate per 1000 cells, matching our 10x v2 chemistry

```
sample      raw_cells pre_dblfilter_cells final_singlets removed_stage2 removed_total_pct
SRR9897621       8281                8169           6885           1284             16.86
SRR9897622      11423               11216           9946           1270             12.93
SRR9897623      27978               27646          23508           4138             15.98
SRR9897624      37537               37378          31020           6358             17.36
SRR9897625      10554               10368           9354           1014             11.37
SRR12539462      7846                7697           7331            366              6.56
SRR12539463     19821               18861          15907           2954             19.75
SRR14615558     17190               16701          14253           2448             17.09

## additional calculations from the output table (assuming 0.8% doublets per 1000 cells per 10x genomics v2)
sample      raw_cells expected_dbl_pct pre_dblfilter_cells removed_dbl_pct final_singlets
SRR9897621  8281      6.62             8169                15.72           6885
SRR9897622  11423     9.14             11216               11.32           9946
SRR9897623  27978     22.38            27646               14.97           23508
SRR9897624  37537     30.03            37378               17.01           31020
SRR9897625  10554     8.44             10368               9.78            9354
SRR12539462 7846      6.28             7697                4.76            7331
SRR12539463 19821     15.86            18861               15.66           15907
SRR14615558 17190     13.75            16701               14.66           14253
```
### post scDblFinder analysis
- SRR...623 (27,978 cells) and SRR...624 (37,537 cells) exceed the standard range of droplet-based microfluidics (10k-20k max usually)
- possibility of under-removal in these overloaded samples - removal % far below expected dbl %
	- overloaded samples may cause multiplets
- scDblfinder focuses on heterotypic doublets
- SRR...621 actual removal % far exceeds expected dbl % 
  - batch/bio factors: incomplete dissoc due to adhesion; multinucleation/polyploidy; singlets with high ambient RNA contamination due to lysis of necrotic cells 
  - dataset1 contamination levels are low according to DecontX from our partial sctk run, so we won't worry about the last factor
- we have removed approx 16% of the initial dataset so far

### scDblFinder histogram diagnostics
- a bimodal distribution is expected, according to documentation
- default settings work well for this dataset
- doublet score of 0.0 = singlets, 1.0 = doublets - in btw are areas of uncertainty
<img width="1937" height="983" alt="scDblFinder_hist_123" src="https://github.com/user-attachments/assets/3776994c-1794-43ed-81f3-bb17f0908c10" />

### post scDblFinder dataset QC metrics
```
orig.ident      n_cells max_count median_count max_feat median_feat max_mt mean_mt mean_ribo
SRR9897621       6885      39986          5639     7730      2217.0  19.98    2.63     16.25
SRR9897622       9946      39937          5183     7981      2361.0  19.95    5.17     20.75
SRR9897623      23508      39605          1104     7099       669.0  19.98    3.48     24.05
SRR9897624      31020      39445          1141     6527       683.0  19.96    6.65     21.26
SRR9897625       9354      38806          8326     6628      2654.5  19.89    4.09     17.07
SRR12539462      7331      39975         12423     6802      3095.0  19.99    3.57     37.12
SRR12539463     15907      39783          2234     7418      1074.0  19.98    6.60     24.12
SRR14615558     14253      39422          1710     6749       853.0  19.97    4.62     25.72
```
### 30 April - all HPCs draining
* [Incident 1583](https://status.alliancecan.ca/view_incident?incident=1583): All compute nodes were down/draining today due to a security patch
* We're migrating to Narval on Alliance Canada after today 
***
## MASS(MCD) + GMM Clustering QC
* this is a clustering based QC to identify outliers after scDblfinder
* the first pass hard filter removes obvious outliers that will distort this step
* intro to MASS + GMM [here](https://github.com/GITC2025/integration_benchmark/wiki/QC_MASS_GMM)
* class 1 is healthy core, class 2 is an extension of the healthy core further away, and class 3 are outliers
* to see how our MASS GMM model works as a filter - we run some diagnostics first
* all 4 QC metrics used - very aggressive filter - 15% cells in class 3 from post scDblfinder dataset
	* full [metrics](https://github.com/GITC2025/integration_benchmark/blob/main/dataset1_QC_code_metrics.md#mass-gmm-class-metrics)

<img width="800" height="400" alt="dataset1_GMM_QC_diagnostics" src="https://github.com/user-attachments/assets/f5433fe9-bc1a-426f-bafa-46a52dfdaa62" />

```
GMM Class 3 Proportions Per Sample and Overall
  orig.ident gmm_distance_class     n percentage
 SRR12539462                  3  2286      31.18
 SRR12539463                  3   976       6.14
 SRR14615558                  3   103       0.72
  SRR9897622                  3  3448      34.67
  SRR9897623                  3  4365      18.57
  SRR9897624                  3  3976      12.82
  SRR9897625                  3  2765      29.56
     OVERALL                  3 17919      15.16
```
```
QC Metric Ranges per GMM Class
 gmm_distance_class nFeature_RNA_min nFeature_RNA_max nCount_RNA_min nCount_RNA_max percent.mt_min percent.mt_max percent.ribo_min percent.ribo_max
                  1              313             7730            500          39986              0          16.31             4.81            47.88
                  2              254             7499            500          39885              0          19.98             1.59            53.36
                  3              200             7981            500          39975              0          19.99             0.91            59.55
```
### analysis of MASS GMM diagnostic
- 558 (paracancerous normal) has the best density of class 1 and class 2 (cells to keep), and only 0.72% class 3
- relatively homogeneous cell pop'n in 558 - stronger genomic instability with most cells in the class 1 and 2
- our MASS GMM class differentiations are clear and align with biological expectations 
  - mt% max increases with each class
  - ribo min decreases (correlated with mt% increase)
  - however, the mito filter here is not appropriate as increasing ribo % has biological significance in cancer context
- with 15% of cells belonging in class 3, this is a fairly aggressive filter
- we can probably adjust this by setting a more lenient mcd (larger viable core), but we need more justification to do so
- **for now, we err on the side of underfiltering and stop at the post scDblfinder dataset, with 16% of cells filtered out**
  - underfiltering is safer as you can preserve rare cell types and QC outliers downstream
  - for variable metrics and seq depth etc, normalization and integration can address this 
- you can rerun this with only [3 QC metrics](https://github.com/GITC2025/integration_benchmark/blob/main/MASS_GMM_3QC%20metrics.md), excluding ribo%, but this makes for an even more aggressive filter, 20.7% of cells in class 3
  - ribo% distribution acts as a stable anchor across cell populations
  - removing ribo% as a variate removes this anchor and results in more outliers in class 3

### median +/- 3MAD
- as a sanity check, you can check if this classic dynamic filter does any good and how it compares with MASS GMM
- still very aggressive, so we'll skip this as well
- additionally, this statistical filter is somewhat blind to the biological state of each sample, whereas the MASS GMM clustering shows some biological differentiation (or batch fx)
- [full metrics](https://github.com/GITC2025/integration_benchmark/blob/main/dataset1_QC_code_metrics.md#median-3mad-filter)
```
3MAD Outlier Proportions Per Sample and Overall
  orig.ident total_cells n_outliers percent_outliers
 SRR12539462        7331       1416            19.32
 SRR12539463       15907       3381            21.25
 SRR14615558       14253       2936            20.60
  SRR9897621        6885       1058            15.37
  SRR9897622        9946       1281            12.88
  SRR9897623       23508       5072            21.58
  SRR9897624       31020       4354            14.04
  SRR9897625        9354        939            10.04
     OVERALL      118204      20437            17.29
```
**no clustering QC was applied for dataset 1**

# Normalization Dataset 1
after removing 16% of cells via hard filters and scDblfinder from original dataset (cellranger filtered matrix), we will run scTransform
***
# Dataset 2 Seurat QC
- we run a similar QC for dataset2, except we need decontX
- so we process both raw and filtered matrices
- [code and metrics here](https://github.com/GITC2025/integration_benchmark/blob/main/dataset2_QC_code_metrics.md)

<img width="800" height="400" alt="dataset2_filtered_QC_by_SRR_clean" src="https://github.com/user-attachments/assets/9f05e69d-9854-4a77-860d-4f451a2115b3" />
<br>
<img width="600" height="400" alt="dataset2_filtered_QC_Merged_Global_clean" src="https://github.com/user-attachments/assets/f46c2cde-af26-4819-8267-24b28e9002f1" />
<br>
<img width="700" height="700" alt="dataset2_QC_Comprehensive_Scatters" src="https://github.com/user-attachments/assets/4480cc13-7c30-4f05-abaf-0a930efaadd3" />
<br>
<img width="700" height="700" alt="dataset2_Global_Knee_Plots_Sandwich" src="https://github.com/user-attachments/assets/67885557-dba1-4261-91ab-99d549caaf0d" />
<br>

# run deContX as first QC step
- we run at [contam <0.2](https://github.com/GITC2025/integration_benchmark/blob/main/dataset2_QC_code_metrics.md#run-decontx-at-contam-02-strict) and [contam <0.6](https://github.com/GITC2025/integration_benchmark/blob/main/dataset2_QC_code_metrics.md#run-decontx-at-contam-06-lenient) to compare attrition rate 
- more info [here](https://github.com/GITC2025/integration_benchmark/wiki/deContX-on-singleron-dataset-2) about deContX for singleron datasets
- contam <0.6 is a suitable threshold, yielding 26% attrition of original filtered matrices
- on the decontamninated dataset, we will then run mito ceilings and other standard filters

<img width="800" height="600" alt="dataset2_decontX_comparison_60_threshold_SRR_fixed - Copy" src="https://github.com/user-attachments/assets/5f80270b-3a46-4caa-9484-2bb65c998afd" />
