# Dataset 1 metrics and plots

## violin plots per sample and global (25 april)
```r
library(Seurat)
library(ggplot2)
library(patchwork)

# make output directory
output_dir <- "/global/scratch/username/R_output_dataset1"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# prep merged object
merged_qc <- merge(dataset1_list[[1]], 
                   y = dataset1_list[2:length(dataset1_list)], 
                   add.cell.ids = names(dataset1_list))

# single group ID for global plot 
merged_qc$all_cells <- "Total Dataset"

# plot 1: individual SRRs
png(filename = file.path(output_dir, "dataset1_QC_by_SRR.png"), 
    width = 2400, height = 1200, res = 300) 

VlnPlot(merged_qc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 4, 
        group.by = "orig.ident", 
        pt.size = 0,
        raster = FALSE) + 
  plot_annotation(title = "Dataset 1 SRP217277 (Lai et al., 2021) - Individual Samples") & 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.x = element_blank())

dev.off()

# plot 2: full dataset 
png(filename = file.path(output_dir, "dataset1_QC_Merged_Global.png"), 
    width = 1800, height = 1200, res = 300)

VlnPlot(merged_qc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 4, 
        group.by = "all_cells", # Forces all 8 SRRs into ONE violin per feature
        pt.size = 0,
        raster = FALSE) + 
  NoLegend() +
  plot_annotation(title = "Dataset 1 SRP217277 (Lai et al., 2021) - Global QC") & 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank())

dev.off()
```

## ranks and thresholds plot (26 april)
```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(scales)

# setup paths
input_path <- "/global/scratch/hpc6140/R_output_dataset1/dataset1_allQCmetrics.rds"
output_dir <- "/global/scratch/hpc6140/R_output_dataset1"

# 1. load and aggregate
dataset1_list <- readRDS(input_path)
all_meta <- bind_rows(lapply(dataset1_list, function(x) x@meta.data))

# save master metadata for downstream analysis
saveRDS(all_meta, file = file.path(output_dir, "dataset1_aggregated_QCmetadata_raw.rds"))

# 2. generate summary metrics with full-spectrum resolution
summary_stats <- all_meta %>%
select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo) %>%
pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
group_by(metric) %>%
summarise(
min = min(value),
q01 = quantile(value, 0.01),
q05 = quantile(value, 0.05),
q10 = quantile(value, 0.10),
q25 = quantile(value, 0.25),
median = median(value),
mean = mean(value),
q75 = quantile(value, 0.75),
q80 = quantile(value, 0.80),
q90 = quantile(value, 0.90),
q95 = quantile(value, 0.95),
q99 = quantile(value, 0.99),
max = max(value),
sd = sd(value)
) %>%
mutate(across(where(is.numeric), ~ round(., 2))) %>%
ungroup()

# print summary to console
cat("\nGlobal QC Quantiles Summary (2 Decimal Places)\n")
print(as.data.frame(summary_stats), row.names = FALSE)
cat("\n")

# save summary stats
write.csv(summary_stats, file = file.path(output_dir, "dataset1_qc_summary_quantiles.csv"), row.names = FALSE)

# 3. generate ranked knee data
get_ranked_data <- function(df, metric) {
vals <- sort(df[[metric]], decreasing = TRUE)
return(data.frame(rank = seq_along(vals), value = vals, metric = metric))
}

knee_data_list <- list(
counts = get_ranked_data(all_meta, "nCount_RNA"),
features = get_ranked_data(all_meta, "nFeature_RNA"),
mt = get_ranked_data(all_meta, "percent.mt"),
ribo = get_ranked_data(all_meta, "percent.ribo")
)

# save ranked data for plot reproducibility
saveRDS(knee_data_list, file = file.path(output_dir, "dataset1_knee_plot_data_ranked.rds"))

# 4. hybrid plotting function with multiple quantile overlays
plot_knee_hybrid <- function(df, title, h_line = NULL, log_scale = TRUE) {
metric_name <- unique(df$metric)

p <- ggplot(df, aes(x = rank, y = value)) +
geom_line(color = "firebrick", linewidth = 1) +
theme_minimal() +
labs(title = title, x = "Rank", y = metric_name)

if(log_scale) {
p <- p + 
scale_x_log10(labels = label_comma(), breaks = breaks_log(n = 6)) + 
scale_y_log10(labels = label_comma(), breaks = breaks_log(n = 6))
} else {
p <- p + scale_x_continuous(labels = label_comma())
}

if(!is.null(h_line)) {
for (val in h_line) {
p <- p + 
geom_hline(yintercept = val, linetype = 1, color = "blue", linewidth = 0.2) +
annotate("text", x = 1, y = val, label = paste0(val), 
vjust = -0.5, hjust = 0, size = 4, color = "blue", fontface = "bold")
}
}

return(p)
}

# 5. individual plot calls with floors and ceilings
p_count <- plot_knee_hybrid(knee_data_list$counts, "nCount_RNA Inflection", h_line = c(500, 40000), log_scale = TRUE)
p_feat <- plot_knee_hybrid(knee_data_list$features, "nFeature_RNA Inflection", h_line = c(200, 8000), log_scale = TRUE)
p_mt <- plot_knee_hybrid(knee_data_list$mt, "percent.mt Distribution", h_line = 20, log_scale = FALSE)
p_ribo <- plot_knee_hybrid(knee_data_list$ribo, "percent.ribo Distribution", h_line = NULL, log_scale = FALSE)

# 6. finalize and save composite
png(filename = file.path(output_dir, "dataset1_Global_Knee_Plots_Sandwich_.png"), 
width = 2400, height = 2400, res = 300)

(p_count | p_feat) / (p_mt | p_ribo) + 
plot_annotation(title = "Dataset 1 Global Ranked QC Metrics (Sandwich Thresholds)") &
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
axis.title = element_text(size = 10),
axis.text = element_text(size = 8))

dev.off()

cat("Diagnostic sandwich plots saved to:", output_dir, "\n")
```

## plot pairwise relationships for QC metrics
```r
library(Seurat)
library(ggplot2)
library(patchwork)

# set output path
output_dir <- "/global/scratch/hpc6140/R_output_dataset1"

# setup plots
plot1 <- FeatureScatter(merged_qc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)
plot2 <- FeatureScatter(merged_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)
plot3 <- FeatureScatter(merged_qc, feature1 = "percent.ribo", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)
plot4 <- FeatureScatter(merged_qc, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)

# define tick intervals
count_breaks <- seq(0, 150000, by = 25000)
pct_breaks <- seq(0, 100, by = 10)

# apply breaks to x-axes
plot1 <- plot1 + scale_x_continuous(breaks = count_breaks)
plot2 <- plot2 + scale_x_continuous(breaks = count_breaks)
plot3 <- plot3 + scale_x_continuous(breaks = pct_breaks)
plot4 <- plot4 + scale_x_continuous(breaks = count_breaks)

# save 2x2 grid with angled labels
png(filename = file.path(output_dir, "dataset1_QC_Comprehensive_Scatters.png"), 
    width = 2400, height = 2400, res = 300)

((plot1 | plot2) / (plot3 | plot4)) + 
  plot_annotation(title = "Dataset 1 SRP217277 (Lai et al., 2021) - Comprehensive QC Scatters") & 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

dev.off()
```
## per sample QC metrics
```r
library(dplyr)
library(purrr)

output_dir <- "/global/scratch/hpc6140/R_output_dataset1"

qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
srr_qc_list <- merged_qc@meta.data %>%
split(.$orig.ident) %>%
map(~ summary(.x[, qc_vars]))

saveRDS(srr_qc_list, file = file.path(output_dir, "dataset1_per_sample_qc_summary_list.rds"))

cat("QC Metric Report: All 8 SRRs\n")
iwalk(srr_qc_list, function(summ, srr_id) {
cat("\nSample ID:", srr_id, "\n")
print(summ)
})
```

```
Sample ID: SRR12539462
   nCount_RNA      nFeature_RNA    percent.mt        percent.ribo
 Min.   :   500   Min.   : 202   Min.   : 0.07325   Min.   : 3.296
 1st Qu.:  7988   1st Qu.:2350   1st Qu.: 2.47945   1st Qu.:34.610
 Median : 12554   Median :3115   Median : 3.31381   Median :39.353
 Mean   : 12507   Mean   :2918   Mean   : 4.05823   Mean   :36.791
 3rd Qu.: 15887   3rd Qu.:3561   3rd Qu.: 4.03932   3rd Qu.:42.333
 Max.   :123970   Max.   :9624   Max.   :71.46995   Max.   :59.554

Sample ID: SRR12539463
   nCount_RNA     nFeature_RNA    percent.mt       percent.ribo
 Min.   :  500   Min.   : 217   Min.   : 0.3596   Min.   : 1.243
 1st Qu.: 1297   1st Qu.: 685   1st Qu.: 4.5013   1st Qu.:21.147
 Median : 2669   Median :1217   Median : 6.1788   Median :24.437
 Mean   : 5011   Mean   :1611   Mean   : 7.7290   Mean   :23.466
 3rd Qu.: 6511   3rd Qu.:2177   3rd Qu.: 8.6192   3rd Qu.:27.412
 Max.   :83166   Max.   :9129   Max.   :93.7974   Max.   :50.696

Sample ID: SRR14615558
   nCount_RNA     nFeature_RNA    percent.mt       percent.ribo
 Min.   :  500   Min.   : 200   Min.   : 0.2926   Min.   : 2.205
 1st Qu.:  961   1st Qu.: 530   1st Qu.: 2.2361   1st Qu.:20.357
 Median : 2166   Median :1005   Median : 3.7874   Median :25.293
 Mean   : 3848   Mean   :1357   Mean   : 5.3620   Mean   :25.315
 3rd Qu.: 4871   3rd Qu.:1830   3rd Qu.: 6.2868   3rd Qu.:30.522
 Max.   :49945   Max.   :7249   Max.   :78.3385   Max.   :49.518

Sample ID: SRR9897621
   nCount_RNA      nFeature_RNA     percent.mt       percent.ribo
 Min.   :   500   Min.   :  205   Min.   : 0.0000   Min.   : 1.595
 1st Qu.:  2356   1st Qu.: 1185   1st Qu.: 0.7987   1st Qu.:12.005
 Median :  7165   Median : 2598   Median : 1.8003   Median :16.452
 Mean   : 11833   Mean   : 2918   Mean   : 3.1639   Mean   :16.328
 3rd Qu.: 15013   3rd Qu.: 4120   3rd Qu.: 3.7469   3rd Qu.:20.452
 Max.   :147367   Max.   :11127   Max.   :82.9004   Max.   :51.592

Sample ID: SRR9897622
   nCount_RNA     nFeature_RNA     percent.mt      percent.ribo
 Min.   :  500   Min.   :  217   Min.   : 0.182   Min.   : 1.504
 1st Qu.: 2598   1st Qu.: 1356   1st Qu.: 3.296   1st Qu.:16.412
 Median : 5679   Median : 2506   Median : 4.530   Median :20.630
 Mean   : 7734   Mean   : 2687   Mean   : 5.604   Mean   :20.567
 3rd Qu.: 9751   3rd Qu.: 3591   3rd Qu.: 6.284   3rd Qu.:24.595
 Max.   :98026   Max.   :10219   Max.   :84.854   Max.   :54.421

Sample ID: SRR9897623
   nCount_RNA     nFeature_RNA    percent.mt      percent.ribo
 Min.   :  498   Min.   : 205   Min.   : 0.000   Min.   : 2.131
 1st Qu.:  799   1st Qu.: 518   1st Qu.: 1.393   1st Qu.:21.403
 Median : 1268   Median : 747   Median : 2.587   Median :24.480
 Mean   : 2482   Mean   :1059   Mean   : 3.778   Mean   :23.877
 3rd Qu.: 2368   3rd Qu.:1206   3rd Qu.: 4.727   3rd Qu.:27.273
 Max.   :58323   Max.   :8289   Max.   :70.560   Max.   :48.506

Sample ID: SRR9897624
   nCount_RNA     nFeature_RNA      percent.mt       percent.ribo
 Min.   :  498   Min.   : 258.0   Min.   : 0.4815   Min.   : 2.182
 1st Qu.:  846   1st Qu.: 537.0   1st Qu.: 4.3743   1st Qu.:18.987
 Median : 1330   Median : 765.0   Median : 5.9524   Median :21.382
 Mean   : 2004   Mean   : 964.8   Mean   : 6.6480   Mean   :21.323
 3rd Qu.: 2259   3rd Qu.:1157.0   3rd Qu.: 8.1566   3rd Qu.:23.739
 Max.   :52961   Max.   :6910.0   Max.   :62.1840   Max.   :59.257

Sample ID: SRR9897625
   nCount_RNA     nFeature_RNA    percent.mt        percent.ribo
 Min.   :  501   Min.   : 226   Min.   : 0.09149   Min.   : 0.8385
 1st Qu.: 5633   1st Qu.:2080   1st Qu.: 2.37985   1st Qu.:13.4331
 Median : 8658   Median :2730   Median : 3.39528   Median :17.5577
 Mean   : 8896   Mean   :2664   Mean   : 4.50149   Mean   :16.9031
 3rd Qu.:11306   3rd Qu.:3257   3rd Qu.: 4.87614   3rd Qu.:21.1085
 Max.   :47087   Max.   :6949   Max.   :66.61952   Max.   :50.7623
```
## global summary
 ```r
output_dir <- "/global/scratch/hpc6140/R_output_dataset1"
qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")

# 1. generate global summary
global_qc_summary <- summary(merged_qc@meta.data[, qc_vars])

# 2. save to scratch
saveRDS(global_qc_summary, file = file.path(output_dir, "dataset1_global_qc_summary.rds"))

# 3. print to console
cat("Global QC Metric Report:\n")
print(global_qc_summary)
```
```
   nCount_RNA      nFeature_RNA     percent.mt      percent.ribo
 Min.   :   498   Min.   :  200   Min.   : 0.000   Min.   : 0.8385
 1st Qu.:  1026   1st Qu.:  609   1st Qu.: 2.660   1st Qu.:18.7797
 Median :  2113   Median : 1070   Median : 4.412   Median :22.6376
 Mean   :  4896   Mean   : 1614   Mean   : 5.477   Mean   :22.7968
 3rd Qu.:  6279   3rd Qu.: 2320   3rd Qu.: 6.811   3rd Qu.:26.5163
 Max.   :147367   Max.   :11127   Max.   :93.797   Max.   :59.5536
```
## high resolution quantiles to justify standard floors/ceilings 
coded together with the ranked plot 
```r
# 2. generate summary metrics with full-spectrum resolution
summary_stats <- all_meta %>%
select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo) %>%
pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
group_by(metric) %>%
summarise(
min = min(value),
q01 = quantile(value, 0.01),
q05 = quantile(value, 0.05),
q10 = quantile(value, 0.10),
q25 = quantile(value, 0.25),
median = median(value),
mean = mean(value),
q75 = quantile(value, 0.75),
q80 = quantile(value, 0.80),
q90 = quantile(value, 0.90),
q95 = quantile(value, 0.95),
q99 = quantile(value, 0.99),
max = max(value),
sd = sd(value)
) %>%
mutate(across(where(is.numeric), ~ round(., 2))) %>%
ungroup()

# print summary to console
cat("\nGlobal QC Quantiles Summary (2 Decimal Places)\n")
print(as.data.frame(summary_stats), row.names = FALSE)
cat("\n")

# save summary stats
write.csv(summary_stats, file = file.path(output_dir, "dataset1_qc_summary_quantiles.csv"), row.names = FALSE)
```
```
Global QC Quantiles Summary (2 Decimal Places)
       metric    min    q01    q05    q10     q25  median    mean     q75
   nCount_RNA 498.00 519.00 598.00 697.00 1026.00 2113.00 4895.58 6279.00
 nFeature_RNA 200.00 350.00 399.00 448.00  609.00 1070.00 1614.11 2320.00
   percent.mt   0.00   0.32   0.89   1.42    2.66    4.41    5.48    6.81
 percent.ribo   0.84   4.51  10.39  14.18   18.78   22.64   22.80   26.52
     q80      q90      q95      q99       max      sd
 7896.00 12542.00 16989.00 31968.42 147367.00 6780.58
 2688.00  3465.00  4275.00  6240.71  11127.00 1347.52
    7.54    10.04    13.14    25.12     93.80    4.87
   27.60    31.21    36.42    43.60     59.55    7.35
```
### calculate attrition from standard floors and ceilings
```r
library(dplyr)

# define ecdf functions
f_count <- ecdf(all_meta$nCount_RNA)
f_feat <- ecdf(all_meta$nFeature_RNA)
f_mt <- ecdf(all_meta$percent.mt)

# thresholds
count_floor <- 500
count_ceil <- 40000
feat_floor <- 200
feat_ceil <- 8000
mt_ceil <- 20

# cumulative probabilities
p_count_floor <- f_count(count_floor)
p_feat_floor <- f_feat(feat_floor)
p_count_ceil <- f_count(count_ceil)
p_feat_ceil <- f_feat(feat_ceil)
p_mt_ceil <- f_mt(mt_ceil)

# output results
cat("\nfloors\n")
cat(paste0("count < 500: ", round(p_count_floor * 100, 4), "%\n"))
cat(paste0("feat < 200:  ", round(p_feat_floor * 100, 4), "%\n"))

cat("\nceilings\n")
cat(paste0("count > 40k: ", round((1 - p_count_ceil) * 100, 4), "%\n"))
cat(paste0("feat > 8k:  ", round((1 - p_feat_ceil) * 100, 4), "%\n"))
cat(paste0("mt > 20%:    ", round((1 - p_mt_ceil) * 100, 4), "%\n"))

cat("\ndataset attrition estimate\n")
total_loss_est <- (p_count_floor + (1 - p_count_ceil) + (1 - p_mt_ceil)) * 100
cat(paste0("total estimated loss: ~", round(total_loss_est, 2), "%\n"))
```
```
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
library(Seurat)

input_path <- "/global/scratch/hpc6140/R_output_dataset1/dataset1_allQCmetrics.rds"
output_dir <- "/global/scratch/hpc6140/R_output_dataset1"

dataset1_list <- readRDS(input_path)

count_floor <- 500
feat_floor <- 200
mt_ceil <- 20

dataset1_pre_dblfilter_list <- lapply(dataset1_list, function(obj) {
subset(obj, subset = nCount_RNA >= count_floor & 
nFeature_RNA >= feat_floor & 
percent.mt < mt_ceil)
})

raw_counts <- sapply(dataset1_list, ncol)
filtered_counts <- sapply(dataset1_pre_dblfilter_list, ncol)

cat("\npre dblfinder filter summary\n")
print(data.frame(
sample = names(dataset1_list),
raw_cells = raw_counts,
pre_dblfilter_cells = filtered_counts,
removed = raw_counts - filtered_counts,
removed_pct = round(((raw_counts - filtered_counts) / raw_counts) * 100, 2)
), row.names = FALSE)

saveRDS(dataset1_pre_dblfilter_list, file = file.path(output_dir, "dataset1_pre_dblfilter_list.rds"))

cat("\npre-dblfilter objects saved successfully.\n")
```
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
### verify range of pre_dblfilter objects
```r
# 1. extract and combine metadata from the filtered list
pre_dbl_meta <- bind_rows(lapply(dataset1_pre_dblfilter_list, function(x) x@meta.data))

# 2. generate tabular summary by orig.ident
qc_verification_table <- pre_dbl_meta %>%
group_by(orig.ident) %>%
summarise(
n_cells = n(),
min_count = min(nCount_RNA),
median_count = median(nCount_RNA),
min_feat = min(nFeature_RNA),
median_feat = median(nFeature_RNA),
max_mt = max(percent.mt),
mean_mt = mean(percent.mt),
mean_ribo = mean(percent.ribo)
) %>%
mutate(across(where(is.numeric), ~ round(., 2)))

cat("\nper-sample qc verification (pre-dblfilter)\n")
print(as.data.frame(qc_verification_table), row.names = FALSE)
```
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

## violin plot for pre_dblfilter (28 april)

```r
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# 1. extract and combine metadata from the filtered list
pre_dbl_meta <- bind_rows(lapply(dataset1_pre_dblfilter_list, function(x) x@meta.data))

# 2. generate tabular summary by orig.ident
qc_verification_table <- pre_dbl_meta %>%
group_by(orig.ident) %>%
summarise(
n_cells = n(),
min_count = min(nCount_RNA),
median_count = median(nCount_RNA),
min_feat = min(nFeature_RNA),
median_feat = median(nFeature_RNA),
max_mt = max(percent.mt),
mean_mt = mean(percent.mt),
mean_ribo = mean(percent.ribo)
) %>%
mutate(across(where(is.numeric), ~ round(., 2)))

cat("\nper-sample qc verification (pre-dblfilter)\n")
print(as.data.frame(qc_verification_table), row.names = FALSE)

# 3. generate violin plots for visual confirmation
v_count <- ggplot(pre_dbl_meta, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
geom_violin(scale = "width", trim = FALSE) +
scale_y_log10(labels = label_comma()) +
theme_minimal() +
theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "nCount_RNA (Floor: 500)", x = NULL)

v_feat <- ggplot(pre_dbl_meta, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
geom_violin(scale = "width", trim = FALSE) +
scale_y_log10(labels = label_comma()) +
theme_minimal() +
theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "nFeature_RNA (Floor: 200)", x = NULL)

v_mt <- ggplot(pre_dbl_meta, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
geom_violin(scale = "width", trim = FALSE) +
theme_minimal() +
theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "percent.mt (Ceiling: < 20)", x = NULL)

v_ribo <- ggplot(pre_dbl_meta, aes(x = orig.ident, y = percent.ribo, fill = orig.ident)) +
geom_violin(scale = "width", trim = FALSE) +
theme_minimal() +
theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "percent.ribo", x = NULL)

# 4. save composite plot
png(filename = file.path(output_dir, "dataset1_pre_dblfilter_QC_violins.png"), width = 2400, height = 1600, res = 300)
(v_count | v_feat) / (v_mt | v_ribo) + 
plot_annotation(title = "Pre-Doublet Filter QC Distributions by Sample") &
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

cat("\npre-dblfilter qc verification plots saved.\n")
```

## run scDblFinder per sample 

```r
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)
library(BiocParallel)

# setup directory and load data
output_dir <- "/global/scratch/hpc6140/R_output_dataset1"
raw_list_path <- file.path(output_dir, "dataset1_allQCmetrics.rds")
pre_list_path <- file.path(output_dir, "dataset1_pre_dblfilter_list.rds")

dataset1_list <- readRDS(raw_list_path)
dataset1_pre_dblfilter_list <- readRDS(pre_list_path)

sample_names <- names(dataset1_pre_dblfilter_list)

# print samples to be processed
cat("\nprocessing the following samples:\n")
print(sample_names)

# utilize multicore processing with a reproducible seed
bp <- MulticoreParam(8, RNGseed=123)

# run doublet detection per sample using documented random approach
dataset1_dblrandom_annotated_list <- bplapply(sample_names, function(s_name) {
seu <- dataset1_pre_dblfilter_list[[s_name]]

# determine doublet rate based on original 10x loading counts
original_cell_count <- ncol(dataset1_list[[s_name]])
expected_dbr <- (original_cell_count / 1000) * 0.008

# convert to singlecellexperiment for compatibility
sce <- as.SingleCellExperiment(seu)

# print sample id before scdblfinder output starts
cat(paste0("\nstarting scdblfinder for sample ", s_name, "\n"))

# run scdblfinder using default normalization and internal feature selection
sce <- scDblFinder(sce, dbr = expected_dbr)

# save diagnostic histogram of doublet scores per sample
png(file.path(output_dir, paste0(s_name, "_scDblFinder_scores.png")))
hist(sce$scDblFinder.score, main=paste("Doublet Scores:", s_name), xlab="Score")
dev.off()

# transfer results back to the seurat object metadata
seu$scDblFinder_class <- sce$scDblFinder.class
seu$scDblFinder_score <- sce$scDblFinder.score
return(seu)
}, BPPARAM = bp)

names(dataset1_dblrandom_annotated_list) <- sample_names

# save the annotated list as rds
saveRDS(dataset1_dblrandom_annotated_list, file = file.path(output_dir, "dataset1_dblrandom_annotated_list.rds"))

# apply biological thresholds and singlet filtering
dataset1_dblrandom_final_list <- lapply(dataset1_dblrandom_annotated_list, function(seu) {
subset(seu, subset = scDblFinder_class == "singlet" &
nCount_RNA < 40000 &
nFeature_RNA < 8000)
})

# save final clean singlets list as rds
saveRDS(dataset1_dblrandom_final_list, file = file.path(output_dir, "dataset1_dblrandom_final_list.rds"))

# generate filtering summary statistics
raw_counts <- sapply(dataset1_list, ncol)
pre_dbl_counts <- sapply(dataset1_pre_dblfilter_list, ncol)
final_counts <- sapply(dataset1_dblrandom_final_list, ncol)

filter_summary <- data.frame(
sample = names(dataset1_list),
raw_cells = raw_counts,
pre_dblfilter_cells = pre_dbl_counts,
final_singlets = final_counts,
removed_stage2 = pre_dbl_counts - final_counts,
removed_total_pct = round(((raw_counts - final_counts) / raw_counts) * 100, 2)
)

# output filtering summary with sample headers
cat("\npost dblfinder filter summary\n")
print(filter_summary, row.names = FALSE)

# aggregate final qc metrics
final_meta <- bind_rows(lapply(dataset1_dblrandom_final_list, function(x) x@meta.data))

final_qc_table <- final_meta %>%
group_by(orig.ident) %>%
summarise(
n_cells = n(),
max_count = max(nCount_RNA),
median_count = median(nCount_RNA),
max_feat = max(nFeature_RNA),
median_feat = median(nFeature_RNA),
max_mt = max(percent.mt),
mean_mt = mean(percent.mt),
mean_ribo = mean(percent.ribo)
) %>%
mutate(across(where(is.numeric), ~ round(., 2)))

# output final qc metrics with sample identifiers
cat("\npost dblfinder qc metrics by sample (orig.ident):\n")
print(as.data.frame(final_qc_table), row.names = FALSE)
```
sample console output
```
Creating ~6158 artificial doublets...
Dimensional reduction
Evaluating kNN...
Training model...
iter=0, 756 cells excluded from training.
iter=1, 768 cells excluded from training.
iter=2, 647 cells excluded from training.
Threshold found:0.449
301 (3.9%) doublets called
```
### scDBlfinder output
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
```
### additional calculations from the output table 
theoretical 0.8% doublets per 1000 cells per 10x genomics v2
```
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
### QC summary post scDblFinder
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
***
## MASS + GMM on scDblfinder filtered dataset
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

# Updated directory for lustre07 scratch
output_dir <- "/lustre07/scratch/delphine/R_output_dataset1"

# Load scDblFinder filtered dataset as the starting point for GMM QC
# Note: Path updated to reflect the existing file confirmed in previous step
hardfiltered <- readRDS(file.path(output_dir, "dataset1_dblrandom_final_list.rds"))

# Merge list if necessary or handle metadata extraction
if (inherits(hardfiltered, "list")) {
hardfiltered <- merge(hardfiltered[[1]], y = hardfiltered[-1], add.cell.ids = names(hardfiltered))
}

meta_df <- hardfiltered@meta.data %>% rownames_to_column("cell_id")
meta_split <- split(meta_df, meta_df$orig.ident)

# Parallel processing with 8 cores
meta_processed <- mclapply(meta_split, function(.x) {

# Log transform metrics for normality in covariance modeling
.x <- .x %>% mutate(
log_nCount = log1p(nCount_RNA),
log_nFeature = log1p(nFeature_RNA),
log_mt = log1p(percent.mt),
log_ribo = log1p(percent.ribo)
)

# Matrix construction for mahalanobis calculation (4D QC metrics)
mat <- as.matrix(.x[, c("log_nCount", "log_nFeature", "log_mt", "log_ribo")])

# MCD to find the robust core of high-quality cells
mcd <- cov.rob(mat, method = "mcd")

# Calculate multivariate distance
.x$mahal_dist <- mahalanobis(mat, center = mcd$center, cov = mcd$cov)
.x$log_mahal_dist <- log1p(.x$mahal_dist)

# GMM fitting (searching for 2-3 clusters of stress/health populations)
gmm_fit <- Mclust(.x$log_mahal_dist, G = 2:3)

.x$gmm_distance_class <- as.character(gmm_fit$classification)
.x$gmm_uncertainty <- gmm_fit$uncertainty

return(.x)

}, mc.cores = 8)

new_meta <- bind_rows(meta_processed) %>% column_to_rownames("cell_id")
hardfiltered <- AddMetaData(hardfiltered, metadata = new_meta)

# Save the GMM annotated object
saveRDS(hardfiltered, file = file.path(output_dir, "dataset1_scDbl_GMM_annotated.rds"))

# Visual diagnostics: GMM classification and Uncertainty
plot_gmm_density <- ggplot(hardfiltered@meta.data, aes(x = log_mahal_dist, y = orig.ident, fill = gmm_distance_class)) +
geom_density_ridges(alpha = 0.7, scale = 0.95) +
theme_minimal() +
labs(title = "GMM classification", x = "log(mahalanobis distance)", y = "sample")

plot_uncertainty <- FeatureScatter(hardfiltered, 
feature1 = "log_mahal_dist", 
feature2 = "gmm_uncertainty", 
group.by = "gmm_distance_class", 
pt.size = 0.5, 
raster = FALSE) +
labs(title = "GMM uncertainty", 
x = "log(mahalanobis distance)", 
y = "uncertainty score")

# Save diagnostic composite
png(filename = file.path(output_dir, "dataset1_GMM_QC_diagnostics.png"), 
width = 2400, height = 1200, res = 300)

(plot_gmm_density | plot_uncertainty) +
plot_annotation(title = "Multivariate QC: MASS (MCD) + GMM Classification") &
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
legend.text = element_text(size = 6),
axis.title = element_text(size = 8),
axis.text = element_text(size = 7))

dev.off()
```

### MASS GMM class metrics
this is the main run with all 4 QC metrics used
```r
# 1. Calculate class percentages per sample
class_stats_sample <- hardfiltered@meta.data %>%
group_by(orig.ident, gmm_distance_class) %>%
summarise(n = n(), .groups = "drop_last") %>%
mutate(percentage = round((n / sum(n)) * 100, 2))

# 2. Calculate class percentages for the overall dataset
class_stats_overall <- hardfiltered@meta.data %>%
group_by(gmm_distance_class) %>%
summarise(n = n(), .groups = "drop") %>%
mutate(
orig.ident = "OVERALL",
percentage = round((n / sum(n)) * 100, 2)
)

# Combine sample and overall stats
class_summary <- bind_rows(class_stats_sample, class_stats_overall)

# 3. Calculate QC ranges (Min/Max) per class for the 4 metrics
qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")

qc_ranges <- hardfiltered@meta.data %>%
group_by(gmm_distance_class) %>%
summarise(across(all_of(qc_metrics), list(
min = ~ min(.),
max = ~ max(.)
), .names = "{.col}_{.fn}")) %>%
mutate(across(where(is.numeric), ~ round(., 2)))

# Save summary tables to CSV
write.csv(class_summary, file = file.path(output_dir, "dataset1_GMM_class_proportions.csv"), row.names = FALSE)
write.csv(qc_ranges, file = file.path(output_dir, "dataset1_GMM_QC_ranges_by_class.csv"), row.names = FALSE)

# Output to console for immediate review
options(width = 1000)
cat("\nGMM Class Proportions Per Sample and Overall\n")
print(as.data.frame(class_summary), row.names = FALSE)

cat("\nQC Metric Ranges per GMM Class\n")
print(as.data.frame(qc_ranges), row.names = FALSE)

cat("\nGMM Class 3 Proportions Per Sample and Overall\n")
print(as.data.frame(class_3_summary), row.names = FALSE)
```
```r
GMM Class Proportions Per Sample and Overall
  orig.ident gmm_distance_class     n percentage
 SRR12539462                  1  3012      41.09
 SRR12539462                  2  2033      27.73
 SRR12539462                  3  2286      31.18
 SRR12539463                  1  9477      59.58
 SRR12539463                  2  5454      34.29
 SRR12539463                  3   976       6.14
 SRR14615558                  1  6690      46.94
 SRR14615558                  2  7460      52.34
 SRR14615558                  3   103       0.72
  SRR9897621                  1  6227      90.44
  SRR9897621                  2   658       9.56
  SRR9897622                  1  3474      34.93
  SRR9897622                  2  3024      30.40
  SRR9897622                  3  3448      34.67
  SRR9897623                  1  9468      40.28
  SRR9897623                  2  9675      41.16
  SRR9897623                  3  4365      18.57
  SRR9897624                  1  9481      30.56
  SRR9897624                  2 17563      56.62
  SRR9897624                  3  3976      12.82
  SRR9897625                  1  2903      31.03
  SRR9897625                  2  3686      39.41
  SRR9897625                  3  2765      29.56
     OVERALL                  1 50732      42.92
     OVERALL                  2 49553      41.92
     OVERALL                  3 17919      15.16
```
```r
QC Metric Ranges per GMM Class
 gmm_distance_class nFeature_RNA_min nFeature_RNA_max nCount_RNA_min nCount_RNA_max percent.mt_min percent.mt_max percent.ribo_min percent.ribo_max
                  1              313             7730            500          39986              0          16.31             4.81            47.88
                  2              254             7499            500          39885              0          19.98             1.59            53.36
                  3              200             7981            500          39975              0          19.99             0.91            59.55
```
```r
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
## Median 3MAD filter
```r
library(Seurat)
library(dplyr)
library(tibble)

# Set directory
output_dir <- "/lustre07/scratch/delphine/R_output_dataset1"

# Load the post-scDblFinder dataset
hardfiltered <- readRDS(file.path(output_dir, "dataset1_dblrandom_final_list.rds"))

# Merge list if necessary to process metadata collectively
if (inherits(hardfiltered, "list")) {
hardfiltered <- merge(hardfiltered[[1]], y = hardfiltered[-1], add.cell.ids = names(hardfiltered))
}

# Extract metadata for calculation
meta_df <- hardfiltered@meta.data %>% rownames_to_column("cell_id")

# Define target QC metrics
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")

# Identify outliers using 3MAD per sample
# This handles the depth heterogeneity observed in samples like SRR...462
outlier_analysis <- meta_df %>%
group_by(orig.ident) %>%
group_modify(~ {
df <- .x
for (metric in qc_metrics) {
m_val <- median(df[[metric]], na.rm = TRUE)
mad_val <- mad(df[[metric]], na.rm = TRUE)
lower_bound <- m_val - (3 * mad_val)
upper_bound <- m_val + (3 * mad_val)
# Identify outliers for specific metric
df[[paste0(metric, "_outlier")]] <- df[[metric]] < lower_bound | df[[metric]] > upper_bound
}
# Cell is an outlier if it fails any one of the 4 metrics
df$is_outlier <- rowSums(df[, grep("_outlier", colnames(df))]) > 0
return(df)
})

# 1. Calculate proportions per sample
sample_outlier_stats <- outlier_analysis %>%
group_by(orig.ident) %>%
summarise(
total_cells = n(),
n_outliers = sum(is_outlier),
percent_outliers = round((n_outliers / total_cells) * 100, 2)
)

# 2. Calculate proportions overall
overall_outlier_stats <- data.frame(
orig.ident = "OVERALL",
total_cells = nrow(outlier_analysis),
n_outliers = sum(outlier_analysis$is_outlier),
percent_outliers = round((sum(outlier_analysis$is_outlier) / nrow(outlier_analysis)) * 100, 2)
)

full_proportion_summary <- bind_rows(sample_outlier_stats, overall_outlier_stats)

# 3. Calculate QC metrics of the outlier population
outlier_qc_summary <- outlier_analysis %>%
filter(is_outlier == TRUE) %>%
group_by(orig.ident) %>%
summarise(across(all_of(qc_metrics), list(
min = ~ min(.),
max = ~ max(.),
median = ~ median(.)
), .names = "{.col}_{.fn}")) %>%
mutate(across(where(is.numeric), ~ round(., 2)))

# Save results as CSV
write.csv(full_proportion_summary, file = file.path(output_dir, "dataset1_3MAD_outlier_proportions.csv"), row.names = FALSE)
write.csv(outlier_qc_summary, file = file.path(output_dir, "dataset1_3MAD_outlier_QC_metrics.csv"), row.names = FALSE)

cat("\n3MAD Outlier Proportions Per Sample and Overall\n")
print(as.data.frame(full_proportion_summary), row.names = FALSE)

cat("\nQC Metrics of 3MAD Outliers (Min/Max/Median)\n")
print(as.data.frame(outlier_qc_summary), row.names = FALSE)
```
```r
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
