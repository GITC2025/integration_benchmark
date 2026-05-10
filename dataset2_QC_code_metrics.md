- we run a similar pipeline for dataset2, but with the addition of running decontX, so we need raw matrix as well
- to run decontX we need to process both filtered and raw matrices
- the raw matrices is a learning model for decontX to learn ambient RNA profile
- decontX then corrects the filtered matrix
- creating a snakemake pipeline can automate the process across datasets

## prep celescope output for seurat
- this is celescope output, we have to adjust accordingly.
- first we need to zip the raw and filtered outputs for seurat to be able to process them
  
```r
system("gzip /lustre07/scratch/delphine/dataset2_R_input/*/*.tsv")
system("gzip /lustre07/scratch/delphine/dataset2_R_input/*/*.mtx")
```

```r
base_path <- "/lustre07/scratch/delphine/dataset2_R_input"
all_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)

for (d in all_dirs) {
cat("DIRECTORY:", basename(d), "\n")
files_to_check <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
for (f in files_to_check) {
file_path <- file.path(d, f)
if (file.exists(file_path)) {
cat("File:", f, "\n")
con <- gzfile(file_path, "rt")
print(readLines(con, n = 5))
close(con)
} else {
cat("File:", f, "- [NOT FOUND]\n")
}
}
cat("\n")
}
```
- we can see this is the GEXSCOPE-V1 pattern with 3 x C8 segments, linkers are not visible
- https://github.com/singleron-RD/CeleScope/blob/master/celescope/chemistry_dict.py
```
DIRECTORY: SRR17259465filtered
File: barcodes.tsv.gz
[1] "AATGTTGC_AAACATCG_AAACATCG" "GTCTGTCA_AACAACCA_AAACATCG"
[3] "CGCTGATC_AACGCTTA_AAACATCG" "AGTGGTCA_AACTCACC_AAACATCG"
[5] "TAGGATGA_AACTCACC_AAACATCG"
File: features.tsv.gz
[1] "ENSG00000142611\tPRDM16\tGene Expression"
[2] "ENSG00000284616\tENSG00000284616\tGene Expression"
[3] "ENSG00000260972\tENSG00000260972\tGene Expression"
[4] "ENSG00000232596\tLINC01646\tGene Expression"
[5] "ENSG00000231510\tLINC02782\tGene Expression"
File: matrix.mtx.gz
[1] "%%MatrixMarket matrix coordinate integer general"
[2] "%"
[3] "55516 10908 17444617"
[4] "43 1 1"
[5] "291 1 1"
```
## process raw matrix
- read10x can still be used - but adjust accordingly to format, e.g. gene is in the 2nd column
- you can also use the generic ReadMtx() but requires manual mapping
```r
library(parallel)
library(Seurat)

base_path <- "/lustre07/scratch/delphine/dataset2_R_input"
raw_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
raw_dirs <- raw_dirs[grep("raw$", raw_dirs)]
sample_names <- basename(raw_dirs)

# read10x but adjusted to celescope output
dataset2_list <- mclapply(seq_along(raw_dirs), function(i) {
counts <- Read10X(data.dir = raw_dirs[i], gene.column = 2)
obj <- CreateSeuratObject(counts = counts, project = sample_names[i])
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent.ribo")
return(obj)
}, mc.cores = 4)

names(dataset2_list) <- sample_names

save_path <- "/lustre07/scratch/delphine/dataset2_R_output/dataset2_raw_QCmetrics.rds"
saveRDS(dataset2_list, file = save_path)

cat("Successfully processed", length(dataset2_list), "samples.\n")
cat("Output saved to:", save_path, "\n")
```
```r
head(dataset2_list[[1]]@meta.data, 5)
```
## raw matrix
- this sparsity/low quality is typical of raw matrices
```
# raw matrices
> head(dataset2_list[[1]]@meta.data, 5)
                           orig.ident nCount_RNA nFeature_RNA percent.mt
AAACATCG_AAACATCG_AAACATCG   AAACATCG        143          112   16.08392
AACAACCA_AAACATCG_AAACATCG   AACAACCA          9            9    0.00000
AACCGAGA_AAACATCG_AAACATCG   AACCGAGA          6            6   33.33333
AACGCTTA_AAACATCG_AAACATCG   AACGCTTA          4            4   25.00000
AACGTGAT_AAACATCG_AAACATCG   AACGTGAT          2            2    0.00000
                           percent.ribo
AAACATCG_AAACATCG_AAACATCG     11.88811
AACAACCA_AAACATCG_AAACATCG      0.00000
AACCGAGA_AAACATCG_AAACATCG      0.00000
AACGCTTA_AAACATCG_AAACATCG     25.00000
AACGTGAT_AAACATCG_AAACATCG      0.00000
```
## process filtered matrix
```r
library(parallel)
library(Seurat)

base_path <- "/lustre07/scratch/delphine/dataset2_R_input"
filtered_dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
filtered_dirs <- filtered_dirs[grep("filtered$", filtered_dirs)]
sample_names <- basename(filtered_dirs)

# Process filtered matrices using gene.column = 2 for Celescope compatibility
dataset2_list_filtered <- mclapply(seq_along(filtered_dirs), function(i) {
counts <- Read10X(data.dir = filtered_dirs[i], gene.column = 2)
obj <- CreateSeuratObject(counts = counts, project = sample_names[i])
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent.ribo")
return(obj)
}, mc.cores = 4)

names(dataset2_list_filtered) <- sample_names

# Save to a dedicated filtered results file
save_path_filtered <- "/lustre07/scratch/delphine/dataset2_R_output/dataset2_filtered_QCmetrics.rds"
saveRDS(dataset2_list_filtered, file = save_path_filtered)

cat("Successfully processed", length(dataset2_list_filtered), "filtered samples.\n")
```

```r
head(dataset2_list_filtered[[1]]@meta.data, 5)
```
in the filtered matrix - the dropouts have been removed through cellranger emptydrops
```
                           orig.ident nCount_RNA nFeature_RNA percent.mt
CGACTGGA_AACAACCA_AAACATCG   CGACTGGA      14400         3941  14.444444
TCCGTCTA_AACAACCA_AAACATCG   TCCGTCTA      12207         3016  25.018432
AAGACGGA_AACGCTTA_AAACATCG   AAGACGGA      16471         3294  34.248072
GTCTGTCA_AACGCTTA_AAACATCG   GTCTGTCA        735          489   8.027211
TCCGTCTA_AACGCTTA_AAACATCG   TCCGTCTA       5530         1471  48.607595
                           percent.ribo
CGACTGGA_AACAACCA_AAACATCG    16.423611
TCCGTCTA_AACAACCA_AAACATCG    21.069878
AAGACGGA_AACGCTTA_AAACATCG    16.574586
GTCTGTCA_AACGCTTA_AAACATCG    19.047619
TCCGTCTA_AACGCTTA_AAACATCG     1.826401
```
## violin plots on filtered matrix
```r
library(Seurat)
library(ggplot2)
library(patchwork)

# fix ID metadata
# This replaces the barcode strings in orig.ident with the SRR IDs stored in the list names
for (i in seq_along(dataset2_list_filtered)) {
dataset2_list_filtered[[i]]$orig.ident <- names(dataset2_list_filtered)[i]
}

# 2. Re-merge the objects with corrected identities
merged_qc_filtered <- merge(dataset2_list_filtered[[1]], 
y = dataset2_list_filtered[2:length(dataset2_list_filtered)], 
add.cell.ids = names(dataset2_list_filtered))

merged_qc_filtered$all_cells <- "Total Dataset"

# 1. Strip the "filtered" suffix from the identity labels
merged_qc_filtered$orig.ident <- gsub("filtered", "", merged_qc_filtered$orig.ident)

# 2. Generate Plot 1: Sample Breakdown (Clean Labels)
png(filename = file.path(output_dir, "dataset2_filtered_QC_by_SRR_clean.png"), 
width = 2400, height = 1200, res = 300) 

VlnPlot(merged_qc_filtered, 
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
ncol = 4, 
group.by = "orig.ident", 
pt.size = 0,
raster = FALSE) + 
plot_annotation(title = "Dataset 2 - Filtered Matrices QC by Sample") & 
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
axis.title.x = element_blank())

dev.off()

# 3. Generate Plot 2: Global Summary
png(filename = file.path(output_dir, "dataset2_filtered_QC_Merged_Global_clean.png"), 
width = 1800, height = 1200, res = 300)

VlnPlot(merged_qc_filtered, 
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
ncol = 4, 
group.by = "all_cells", 
pt.size = 0,
raster = FALSE) + 
NoLegend() +
plot_annotation(title = "Dataset 2 - Global QC Summary") & 
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
axis.text.x = element_text(size = 12),
axis.title.x = element_blank())

dev.off()
```
## rank and quantiles and metrics on filtered
```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(scales)

# Setup paths for Dataset 2
output_dir <- "/lustre07/scratch/delphine/dataset2_R_output"

all_meta <- merged_qc_filtered@meta.data

saveRDS(all_meta, file = file.path(output_dir, "dataset2_aggregated_QCmetadata_filtered.rds"))

# 2. Generate summary metrics with full-spectrum resolution
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

# Print summary to console
options(width = 1000)
cat("\nDataset 2 Global QC Quantiles Summary\n")
print(as.data.frame(summary_stats), row.names = FALSE)

# Save summary stats
write.csv(summary_stats, file = file.path(output_dir, "dataset2_qc_summary_quantiles.csv"), row.names = FALSE)

# 3. Generate ranked knee data
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

# Save ranked data
saveRDS(knee_data_list, file = file.path(output_dir, "dataset2_knee_plot_data_ranked.rds"))

# 4. Hybrid plotting function
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

# 5. Individual plot calls (Thresholds can be adjusted based on Dataset 2 summary_stats)
p_count <- plot_knee_hybrid(knee_data_list$counts, "nCount_RNA Inflection", h_line = c(500, 40000), log_scale = TRUE)
p_feat <- plot_knee_hybrid(knee_data_list$features, "nFeature_RNA Inflection", h_line = c(200, 8000), log_scale = TRUE)
p_mt <- plot_knee_hybrid(knee_data_list$mt, "percent.mt Distribution", h_line = 20, log_scale = FALSE)
p_ribo <- plot_knee_hybrid(knee_data_list$ribo, "percent.ribo Distribution", h_line = NULL, log_scale = FALSE)

# 6. Finalize and save composite
png(filename = file.path(output_dir, "dataset2_Global_Knee_Plots_Sandwich.png"), 
width = 2400, height = 2400, res = 300)

(p_count | p_feat) / (p_mt | p_ribo) + 
plot_annotation(title = "Dataset 2 Global Ranked QC Metrics (Sandwich Thresholds)") &
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
axis.title = element_text(size = 10),
axis.text = element_text(size = 8))

dev.off()

cat("Diagnostic sandwich plots saved to:", output_dir, "\n")
```
## global quantiles summary
```
Dataset 2 Global QC Quantiles Summary
       metric    min    q01    q05     q10     q25  median    mean     q75     q80      q90      q95      q99      max      sd
   nCount_RNA 500.00 535.00 715.40 1039.00 2234.00 4491.00 5894.27 8186.00 9275.00 12385.00 15663.00 23212.80 49721.00 4965.22
 nFeature_RNA  61.00 337.00 430.00  548.00  890.00 1562.00 1750.13 2459.00 2674.00  3212.00  3649.60  4515.00  7175.00 1031.65
   percent.mt   0.20   3.38   7.19    9.44   14.26   20.70   24.69   33.00   36.19    45.58    53.70    66.05    95.65   14.40
 percent.ribo   0.26   1.06   1.93    2.84    5.97   10.24   10.33   13.93   14.84    17.55    20.04    24.46    46.55    5.57
 ```
## pairwise QC metrics
```r
library(Seurat)
library(ggplot2)
library(patchwork)

output_dir <- "/lustre07/scratch/delphine/dataset2_R_output"

# Setup plots
plot1 <- FeatureScatter(merged_qc_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)
plot2 <- FeatureScatter(merged_qc_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)
plot3 <- FeatureScatter(merged_qc_filtered, feature1 = "percent.ribo", feature2 = "percent.mt", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)
plot4 <- FeatureScatter(merged_qc_filtered, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by = "orig.ident", pt.size = 0.5, raster = FALSE)

# Intervals
count_breaks <- seq(0, 50000, by = 10000)
pct_breaks <- seq(0, 100, by = 10)

# Apply breaks
plot1 <- plot1 + scale_x_continuous(breaks = count_breaks)
plot2 <- plot2 + scale_x_continuous(breaks = count_breaks)
plot3 <- plot3 + scale_x_continuous(breaks = pct_breaks)
plot4 <- plot4 + scale_x_continuous(breaks = count_breaks)

# Save 2x2 grid with enlarged legends
png(filename = file.path(output_dir, "dataset2_QC_Comprehensive_Scatters.png"), 
width = 2400, height = 2400, res = 300)

((plot1 | plot2) / (plot3 | plot4)) + 
plot_annotation(title = "Dataset 2 - Filtered Comprehensive QC Scatters") & 
theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
legend.title = element_text(size = 14),
legend.text = element_text(size = 12),
axis.title = element_text(size = 8),
axis.text = element_text(size = 7),
axis.text.x = element_text(angle = 45, hjust = 1)) &
guides(color = guide_legend(override.aes = list(size = 5)))

dev.off()
```
## run deContX at contam <0.2 (strict)
```r
## version 2 w doc input - contam 0.2 
```r
library(Seurat)
library(celda)
library(dplyr)
library(tibble)

input_dir <- "/lustre07/scratch/delphine/dataset2_R_input"
output_dir <- "/lustre07/scratch/delphine/dataset2_R_output"
samples <- c("SRR17259462", "SRR17259463", "SRR17259464", "SRR17259465")

seurat_original_list <- list()
seurat_decont_list <- list()

for (s in samples) {
filt_dir <- file.path(input_dir, paste0(s, "filtered"))
raw_dir <- file.path(input_dir, paste0(s, "raw"))

counts_filt <- Read10X(data.dir = filt_dir)
counts_raw <- Read10X(data.dir = raw_dir)

sobj_orig <- CreateSeuratObject(counts = counts_filt, project = s)
sobj_orig[["percent.mt"]] <- PercentageFeatureSet(sobj_orig, pattern = "^MT-")
sobj_orig[["percent.ribo"]] <- PercentageFeatureSet(sobj_orig, pattern = "^RP[SL]")
seurat_original_list[[s]] <- sobj_orig

decont_res <- decontX(x = counts_filt, background = counts_raw)

sobj_decont <- CreateSeuratObject(counts = round(decont_res$decontXcounts), project = s)
sobj_decont <- AddMetaData(sobj_decont, metadata = decont_res$contamination, col.name = "decontX_contamination")
sobj_decont[["percent.mt"]] <- PercentageFeatureSet(sobj_decont, pattern = "^MT-")
sobj_decont[["percent.ribo"]] <- PercentageFeatureSet(sobj_decont, pattern = "^RP[SL]")
seurat_decont_list[[s]] <- sobj_decont
}

merged_orig <- merge(seurat_original_list[[1]], y = seurat_original_list[-1], add.cell.ids = samples)
merged_decont <- merge(seurat_decont_list[[1]], y = seurat_decont_list[-1], add.cell.ids = samples)

saveRDS(merged_decont, file = file.path(output_dir, "dataset2_decontX_merged_unfiltered.rds"))

contam_threshold <- 0.20
cells_to_keep <- rownames(merged_decont@meta.data[merged_decont@meta.data$decontX_contamination < contam_threshold, ])
merged_decont_filtered <- subset(merged_decont, cells = cells_to_keep)

total_orig <- ncol(merged_orig)
total_retained <- ncol(merged_decont_filtered)
attrition_rate <- ((total_orig - total_retained) / total_orig) * 100

cat("\nAttrition Summary\n")
cat("Total Original Cells:", total_orig, "\n")
cat("Total Retained Cells (Contamination <", contam_threshold, "):", total_retained, "\n")
cat("Attrition Rate:", round(attrition_rate, 2), "%\n\n")

get_stats <- function(df, stage) {
stats <- df %>% 
summarise(
Stage = stage,
nCount_min = min(nCount_RNA), nCount_max = max(nCount_RNA), nCount_median = median(nCount_RNA), nCount_mean = mean(nCount_RNA),
nFeature_min = min(nFeature_RNA), nFeature_max = max(nFeature_RNA), nFeature_median = median(nFeature_RNA), nFeature_mean = mean(nFeature_RNA),
mt_min = min(percent.mt), mt_max = max(percent.mt), mt_median = median(percent.mt), mt_mean = mean(percent.mt),
ribo_min = min(percent.ribo), ribo_max = max(percent.ribo), ribo_median = median(percent.ribo), ribo_mean = mean(percent.ribo)
)
return(stats)
}

stats_orig <- get_stats(merged_orig@meta.data, "Before_decontX_Filter")
stats_decont <- get_stats(merged_decont_filtered@meta.data, "After_decontX_Filter")
summary_stats <- bind_rows(stats_orig, stats_decont)

write.csv(summary_stats, file = file.path(output_dir, "dataset2_decontX_attrition_stats.csv"), row.names = FALSE)

cat("Dataset 2 decontX Dry Run QC Stats\n")
print(as.data.frame(summary_stats), row.names = FALSE)
```
contam <0.2 is too aggressive for this dataset
- over half of captured droplets contain at least 20% ambient RNA
- it doesn't mean we have to filter out more than half our dataset
- all cells will be cleaned up by decontX, improving gene profiles
```
Attrition Summary
Total Original Cells: 33909
Total Retained Cells (Contamination < 0.2 ): 15562
Attrition Rate: 54.11 %

Dataset 2 decontX Dry Run QC Stats

Stage                  nCount_min nCount_max nCount_median nCount_mean nFeature_min nFeature_max nFeature_median nFeature_mean mt_min    mt_max   mt_median mt_mean  ribo_min  ribo_max  ribo_median ribo_mean
Before_decontX_Filter  500        49721      4491          5894.274    61           7175         1562            1750.134      0.1984127 95.65406 20.69648  24.68912 0.2563164 46.54896  10.23680    10.33034
After_decontX_Filter   448        42366      4032          5651.807    60           7099         1489            1769.102      0.1996008 94.45430 16.49805  21.67215 0.2487562 44.23696  10.89439    10.18093
```
## run deContX at contam <0.6 (lenient)
```r
library(Seurat)
library(celda)
library(dplyr)
library(tibble)

input_dir <- "/lustre07/scratch/delphine/dataset2_R_input"
output_dir <- "/lustre07/scratch/delphine/dataset2_R_output"
samples <- c("SRR17259462", "SRR17259463", "SRR17259464", "SRR17259465")

seurat_original_list <- list()
seurat_decont_list <- list()

for (s in samples) {
filt_dir <- file.path(input_dir, paste0(s, "filtered"))
raw_dir <- file.path(input_dir, paste0(s, "raw"))

counts_filt <- Read10X(data.dir = filt_dir)
counts_raw <- Read10X(data.dir = raw_dir)

sobj_orig <- CreateSeuratObject(counts = counts_filt, project = s)
sobj_orig[["percent.mt"]] <- PercentageFeatureSet(sobj_orig, pattern = "^MT-")
sobj_orig[["percent.ribo"]] <- PercentageFeatureSet(sobj_orig, pattern = "^RP[SL]")
seurat_original_list[[s]] <- sobj_orig

decont_res <- decontX(x = counts_filt, background = counts_raw)

sobj_decont <- CreateSeuratObject(counts = round(decont_res$decontXcounts), project = s)
sobj_decont <- AddMetaData(sobj_decont, metadata = decont_res$contamination, col.name = "decontX_contamination")
sobj_decont[["percent.mt"]] <- PercentageFeatureSet(sobj_decont, pattern = "^MT-")
sobj_decont[["percent.ribo"]] <- PercentageFeatureSet(sobj_decont, pattern = "^RP[SL]")
seurat_decont_list[[s]] <- sobj_decont
}

merged_orig <- merge(seurat_original_list[[1]], y = seurat_original_list[-1], add.cell.ids = samples)
merged_decont <- merge(seurat_decont_list[[1]], y = seurat_decont_list[-1], add.cell.ids = samples)

saveRDS(merged_decont, file = file.path(output_dir, "dataset2_decontX_merged_unfiltered.rds"))

contam_threshold <- 0.60
cells_to_keep <- rownames(merged_decont@meta.data[merged_decont@meta.data$decontX_contamination < contam_threshold, ])
merged_decont_filtered <- subset(merged_decont, cells = cells_to_keep)

total_orig <- ncol(merged_orig)
total_retained <- ncol(merged_decont_filtered)
attrition_rate <- ((total_orig - total_retained) / total_orig) * 100

cat("\nAttrition Summary\n")
cat("Total Original Cells:", total_orig, "\n")
cat("Total Retained Cells (Contamination <", contam_threshold, "):", total_retained, "\n")
cat("Attrition Rate:", round(attrition_rate, 2), "%\n\n")

get_stats <- function(df, stage) {
stats <- df %>% 
summarise(
Stage = stage,
nCount_min = min(nCount_RNA), nCount_max = max(nCount_RNA), nCount_median = median(nCount_RNA), nCount_mean = mean(nCount_RNA),
nFeature_min = min(nFeature_RNA), nFeature_max = max(nFeature_RNA), nFeature_median = median(nFeature_RNA), nFeature_mean = mean(nFeature_RNA),
mt_min = min(percent.mt), mt_max = max(percent.mt), mt_median = median(percent.mt), mt_mean = mean(percent.mt),
ribo_min = min(percent.ribo), ribo_max = max(percent.ribo), ribo_median = median(percent.ribo), ribo_mean = mean(percent.ribo)
)
return(stats)
}

stats_orig <- get_stats(merged_orig@meta.data, "Before_decontX_Filter")
stats_decont <- get_stats(merged_decont_filtered@meta.data, "After_decontX_Filter")
summary_stats <- bind_rows(stats_orig, stats_decont)

write.csv(summary_stats, file = file.path(output_dir, "dataset2_decontX_attrition_stats_60_threshold.csv"), row.names = FALSE)

options(width = 1000)
cat("Dataset 2 decontX Dry Run QC Stats\n")
print(as.data.frame(summary_stats), row.names = FALSE)
```
contam <0.6 is reasonable for this dataset
- 26% of cells have ambient RNA contam >0.6
- these cells are usually removed due to high contamination, although all cells are cleaned in deContX
```
Attrition Summary
Total Original Cells: 33909
Total Retained Cells (Contamination < 0.6 ): 24964
Attrition Rate: 26.38 %

Dataset 2 decontX Dry Run QC Stats
                 Stage nCount_min nCount_max nCount_median nCount_mean nFeature_min nFeature_max nFeature_median nFeature_mean    mt_min   mt_max mt_median  mt_mean  ribo_min ribo_max ribo_median ribo_mean
 Before_decontX_Filter        500      49721        4491.0    5894.274           61         7175            1562      1750.134 0.1984127 95.65406  20.69648 24.68912 0.2563164 46.54896    10.23680  10.33034
  After_decontX_Filter        167      42366        3605.5    4994.588           60         7099            1399      1654.708 0.1996008 94.45430  17.07233 23.17881 0.0000000 44.23696    10.70319  10.07724
```
## violin plot post deContX (<0.6)
```r
library(Seurat)
library(ggplot2)
library(patchwork)

# fix metadata
fix_meta <- function(obj) {
obj$orig.ident <- sapply(strsplit(colnames(obj), "_"), `[`, 1)
obj$orig.ident <- gsub("filtered", "", obj$orig.ident)
obj$orig.ident <- factor(obj$orig.ident, levels = c("SRR17259462", "SRR17259463", "SRR17259464", "SRR17259465"))
return(obj)
}

merged_decont <- fix_meta(merged_decont)
merged_decont_filtered <- fix_meta(merged_decont_filtered)

# sync y axis
metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
y_lims <- list()
for (m in metrics) {
y_lims[[m]] <- c(0, max(merged_decont@meta.data[[m]], na.rm = TRUE) * 1.05)
}

# plot row helper
make_qc_row <- function(obj, row_label) {
p_list <- VlnPlot(
obj, 
features = metrics, 
ncol = 4, 
group.by = "orig.ident", 
pt.size = 0, 
raster = FALSE
) & 
theme(
plot.title = element_text(size = 10, face = "bold"),
plot.subtitle = element_text(size = 9, face = "italic"),
axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
axis.title.x = element_blank(),
legend.position = "none"
)

# apply scales to match original violin plot
for (i in 1:4) {
p_list[[i]] <- p_list[[i]] + coord_cartesian(ylim = y_lims[[metrics[i]]])
if (i == 1) {
p_list[[i]] <- p_list[[i]] + labs(subtitle = row_label)
}
}
return(p_list)
}

png(filename = file.path(output_dir, "dataset2_decontX_comparison_60_threshold_SRR_fixed.png"), 
width = 3200, height = 2400, res = 300)

top_row <- make_qc_row(merged_decont, "top: all cells (post-decontx)")
bottom_row <- make_qc_row(merged_decont_filtered, "bottom: retained cells (contam < 0.6)")

combined_plot <- top_row / bottom_row + 
plot_annotation(
title = "dataset 2: decontamination filtering by sample id",
theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
)

print(combined_plot)
dev.off()
```
