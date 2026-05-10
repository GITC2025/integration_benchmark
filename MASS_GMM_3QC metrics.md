- A second run with the MASS GMM QC clustering on post scDblfinder dataset with only 3 QC metrics (mito, count, feature) and without ribo%
- leads to more aggressive filtering due to the loss of ribo% as anchor for the core

```r
library(Seurat)
library(MASS)
library(mclust)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ggridges)
library(parallel)

output_dir <- "/lustre07/scratch/delphine/R_output_dataset1"

hardfiltered <- readRDS(file.path(output_dir, "dataset1_dblrandom_final_list.rds"))

if (inherits(hardfiltered, "list")) {
hardfiltered <- merge(hardfiltered[[1]], y = hardfiltered[-1], add.cell.ids = names(hardfiltered))
}

meta_df <- hardfiltered@meta.data %>% rownames_to_column("cell_id")
meta_split <- split(meta_df, meta_df$orig.ident)

meta_processed <- mclapply(meta_split, function(.x) {

.x <- .x %>% mutate(
log_nCount = log1p(nCount_RNA),
log_nFeature = log1p(nFeature_RNA),
log_mt = log1p(percent.mt)
)

# Matrix construction (3D: nCount, nFeature, mt)
mat <- as.matrix(.x[, c("log_nCount", "log_nFeature", "log_mt")])

# MCD for core
mcd <- cov.rob(mat, method = "mcd")

# Calculate multivariate distance
.x$mahal_dist <- mahalanobis(mat, center = mcd$center, cov = mcd$cov)
.x$log_mahal_dist <- log1p(.x$mahal_dist)

# GMM fitting 
gmm_fit <- Mclust(.x$log_mahal_dist, G = 2:3)

.x$gmm_distance_class <- as.character(gmm_fit$classification)
.x$gmm_uncertainty <- gmm_fit$uncertainty

return(.x)

}, mc.cores = 8)

new_meta <- bind_rows(meta_processed) %>% column_to_rownames("cell_id")
hardfiltered <- AddMetaData(hardfiltered, metadata = new_meta)

saveRDS(hardfiltered, file = file.path(output_dir, "dataset1_scDbl_GMM_3D_annotated.rds"))

plot_gmm_density <- ggplot(hardfiltered@meta.data, aes(x = log_mahal_dist, y = orig.ident, fill = gmm_distance_class)) +
geom_density_ridges(alpha = 0.7, scale = 0.95) +
theme_minimal() +
labs(title = "3D GMM classification", x = "log(mahalanobis distance)", y = "sample")

plot_uncertainty <- FeatureScatter(hardfiltered, 
feature1 = "log_mahal_dist", 
feature2 = "gmm_uncertainty", 
group.by = "gmm_distance_class", 
pt.size = 0.5, 
raster = FALSE) +
labs(title = "3D GMM uncertainty", 
x = "log(mahalanobis distance)", 
y = "uncertainty score")

png(filename = file.path(output_dir, "dataset1_GMM_3D_QC_diagnostics.png"), 
width = 2400, height = 1200, res = 300)

(plot_gmm_density | plot_uncertainty) +
plot_annotation(title = "Multivariate QC: 3D MASS (MCD) + GMM Classification") &
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
legend.text = element_text(size = 6),
axis.title = element_text(size = 8),
axis.text = element_text(size = 7))

dev.off()

# Summary Statistics
class_stats_sample <- hardfiltered@meta.data %>%
group_by(orig.ident, gmm_distance_class) %>%
summarise(n = n(), .groups = "drop_last") %>%
mutate(percentage = round((n / sum(n)) * 100, 2))

class_stats_overall <- hardfiltered@meta.data %>%
group_by(gmm_distance_class) %>%
summarise(n = n(), .groups = "drop") %>%
mutate(orig.ident = "OVERALL", percentage = round((n / sum(n)) * 100, 2))

class_summary <- bind_rows(class_stats_sample, class_stats_overall)

qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
qc_ranges <- hardfiltered@meta.data %>%
group_by(gmm_distance_class) %>%
summarise(across(all_of(qc_metrics), list(
min = ~ min(.),
max = ~ max(.)
), .names = "{.col}_{.fn}")) %>%
mutate(across(where(is.numeric), ~ round(., 2)))

write.csv(class_summary, file = file.path(output_dir, "dataset1_GMM_3D_class_proportions.csv"), row.names = FALSE)
write.csv(qc_ranges, file = file.path(output_dir, "dataset1_GMM_3D_QC_ranges_by_class.csv"), row.names = FALSE)

# console print
options(width = 1000)
cat("\nGMM Class Proportions Per Sample and Overall\n")
print(as.data.frame(class_summary), row.names = FALSE)

cat("\nQC Metric Ranges per GMM Class\n")
print(as.data.frame(qc_ranges), row.names = FALSE)

cat("\nGMM Class 3 Proportions Per Sample and Overall\n")
print(as.data.frame(class_3_summary), row.names = FALSE)
```
```r
GMM Class 3 Proportions Per Sample and Overall
  orig.ident gmm_distance_class     n percentage
 SRR12539462                  3  2663      36.33
 SRR12539463                  3  4733      29.75
 SRR14615558                  3  4744      33.28
  SRR9897622                  3  2250      22.62
  SRR9897623                  3  4486      19.08
  SRR9897624                  3  3205      10.33
  SRR9897625                  3  2350      25.12
     OVERALL                  3 24431      20.67
```
```r
QC Metric Ranges per GMM Class
 gmm_distance_class nFeature_RNA_min nFeature_RNA_max nCount_RNA_min nCount_RNA_max percent.mt_min percent.mt_max percent.ribo_min percent.ribo_max
                  1              344             7339            503          36103           0.06          10.54             1.59            48.06
                  2              254             7981            500          39986           0.00          19.98             1.64            51.59
                  3              200             7725            500          39975           0.00          19.99             0.91            59.55
```

```r
GMM Class Proportions Per Sample and Overall
  orig.ident gmm_distance_class     n percentage
 SRR12539462                  1  2407      32.83
 SRR12539462                  2  2261      30.84
 SRR12539462                  3  2663      36.33
 SRR12539463                  1  5461      34.33
 SRR12539463                  2  5713      35.92
 SRR12539463                  3  4733      29.75
 SRR14615558                  1  5007      35.13
 SRR14615558                  2  4502      31.59
 SRR14615558                  3  4744      33.28
  SRR9897621                  1  4906      71.26
  SRR9897621                  2  1979      28.74
  SRR9897622                  1  2435      24.48
  SRR9897622                  2  5261      52.90
  SRR9897622                  3  2250      22.62
  SRR9897623                  1  6278      26.71
  SRR9897623                  2 12744      54.21
  SRR9897623                  3  4486      19.08
  SRR9897624                  1  7762      25.02
  SRR9897624                  2 20053      64.65
  SRR9897624                  3  3205      10.33
  SRR9897625                  1  2586      27.65
  SRR9897625                  2  4418      47.23
  SRR9897625                  3  2350      25.12
     OVERALL                  1 36842      31.17
     OVERALL                  2 56931      48.16
     OVERALL                  3 24431      20.67
```

