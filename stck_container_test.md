# custom vs official sctk docker comparison 

* test run on PBMC 1K dataset on the HPC
* debugging underlying issues in the official sctk_runQC script and official sctk 2.14.0 container
* refining our custom container sctk 2.18.0 built on a docker base of R 4.5.3
* many big bugs solved, but small bugs remain as of 23 april 2026 - debugging continues
* apptainer = docker for HPC (but without root access, only fakeroot)
* ensure your working node/job node matches CPU specs according to documentation
* https://camplab.net/sctk/current/articles/cmd_qc.html#running-sctk-qc-with-docker 

## bugs found in official container/script sctk 2.14.0 (and official script for sctk 2.18.0)
* designed for a permissive local environment - does not work on restrictive HPCs
  * attempts to download dependencies during runtime near the end
  * in our custom container we create a sufficient python venv and lockdown against downloading - solves this issue
  * container is entirely self-sufficient and locks versions
* unable to detect mitochondrial genes even with correct settings - cellrangerV3, human symbol for PBMC 1k cellranger v3 dataset
  * same bug in the official script for sctk 2.18.0 - added shims for fallback mito detection and calculation from raw counts
  * will not show up in overall metrics even with shim
* Scrublet does not work
* Seurat 5 bridging issue - we patched this
* our custom container with 5 shims in the script solves most of these issues - but requires more refinement

## official campbio sctk 2.14.0 docker
```bash
apptainer pull /global/scratch/$USER/sctk_v2.14.0.sif docker://campbio/sctk_qc:latest
```
```bash
# check versions in the latest container
apptainer exec sctk_qc_latest.sif R -e "packageVersion('singleCellTK')"

R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
...
> packageVersion('singleCellTK')
[1] ‘2.14.0’
```
```bash
apptainer exec --cleanenv /global/scratch/$USER/sctk_v.2.14.0.sif Rscript -e '
pkgs <- as.data.frame(installed.packages())
pkg_versions <- paste(pkgs$Package, pkgs$Version, sep=" ")
write.table(sort(pkg_versions), file="sctk_v2.14.0_deps.txt", row.names=F, col.names=F, quote=F)
'
```
[sctk_v2.14.0_deps.txt](https://github.com/user-attachments/files/27032839/sctk_v2.14.0_deps.txt)

### test run PBMC 1K output from official sctk_v2.14.0.sif
[pbmc_sctkv2.14.0.txt](https://github.com/user-attachments/files/27032910/pbmc_sctkv2.14.0.txt)
* execution halted when R attempted to download py dependencies but SSL certs were not valid
* a container should ideally contain all its own dependencies 

## our custom container sctk v2.18.0 built on R 4.5.3
[building instructions](https://github.com/GITC2025/integration_benchmark/blob/main/custom_sctk_container.md)

[sctk_v2.18.0_deps.txt](https://github.com/user-attachments/files/27032840/sctk_v2.18.0_deps.txt)

[sctk_v2.18.0_tree_23april.txt](https://github.com/user-attachments/files/27032884/sctk_v2.18.0_tree_23april.txt)

### shimmied R script for custom container, built on official sctk script
[sctk_mito3_newnet_23april.R](https://github.com/GITC2025/integration_benchmark/blob/main/sctk_mito3_newnet_23april.R)

[sctk v2.18.0 original script](https://raw.githubusercontent.com/compbiomed/singleCellTK/v2.18.0/exec/SCTK_runQC.R)

### test run PBMC 1K output from custom container and custom script

full output log here

```bash
apptainer exec \
--bind /global/scratch/$USER:/global/scratch/$USER \
/global/scratch/$USER/sctk_v2.18.0_R_v4.5.3.sif \
Rscript /global/scratch/$USER/sctk_mito3_newnet_23april.R \
-P CellRangerV3 \
-R /usr/local/bin/sctk_mito3_newnet_23april.R \
-C /global/scratch/$USER/test_sctk/filtered_feature_bc_matrix \
-s pbmc_1k \
-o /global/scratch/$USER/test_sctk/pbmc_qc_mito3 \
-d Both \
-S TRUE \
-F SCE,Seurat,AnnData,FlatFile \
-M TRUE \
-E human-symbol \
-n 8 \
2>&1 | tee pmbc_1k_mito3.txt
```

```bash
cat SCTK_pbmc_1k_cellQC_summary.csv
```
* good match with online data https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_web_summary.html 
* from here you can see Scrublet did not work, mito gene check did not work
* similar issues with official container, though we debugged some issues
* this will affect downstream normalization

```bash
# from our custom container

"","All Samples"
"Number of Cells",1222
"Mean counts",7557.4
"Median counts",6628
"Mean features detected",2044.4
"Median features detected",1919
"Scrublet - Number of doublets",0
"Scrublet - Percentage of doublets",0
"scDblFinder - Number of doublets",41
"scDblFinder - Percentage of doublets",3.36
"DoubletFinder - Number of doublets, Resolution 1.5",92
"DoubletFinder - Percentage of doublets, Resolution 1.5",7.53
"CXDS - Number of doublets",178
"CXDS - Percentage of doublets",14.6
"BCDS - Number of doublets",96
"BCDS - Percentage of doublets",7.86
"SCDS Hybrid - Number of doublets",104
"SCDS Hybrid - Percentage of doublets",8.51
"DecontX - Mean contamination",0.0399
"DecontX - Median contamination",0.00977
```
* Note that the Python folders are unfortunately empty, as the script - whether original or modified with the safetynet still runs into issues converting SCE to h5ad. while it's possible to write a conversion patch into the script, it may be best to convert manually as needed, and check for data integrity. a common issue across all tools.

```bash
├── FlatFile
│   ├── Cells
│   │   ├── assays
│   │   │   ├── pbmc_1k_counts.mtx.gz
│   │   │   ├── pbmc_1k_decontXcounts_bg.mtx.gz
│   │   │   ├── pbmc_1k_decontXcounts.mtx.gz
│   │   │   ├── pbmc_1k_SoupX_bg.mtx.gz
│   │   │   └── pbmc_1k_SoupX.mtx.gz
│   │   ├── metadata
│   │   │   └── pbmc_1k_metadata.rds
│   │   ├── pbmc_1k_cellData.txt.gz
│   │   ├── pbmc_1k_featureData.txt.gz
│   │   └── reducedDims
│   │       ├── pbmc_1k_decontX_UMAP_bg.txt.gz
│   │       ├── pbmc_1k_decontX_UMAP.txt.gz
│   │       ├── pbmc_1k_doubletFinder_UMAP.txt.gz
│   │       ├── pbmc_1k_scrublet_TSNE.txt.gz
│   │       ├── pbmc_1k_scrublet_UMAP.txt.gz
│   │       ├── pbmc_1k_SoupX_bg_UMAP_all_cells.txt.gz
│   │       └── pbmc_1k_SoupX_UMAP_all_cells.txt.gz
│   └── Droplets
│       ├── assays
│       │   └── pbmc_1k_counts.mtx.gz
│       ├── metadata
│       │   └── pbmc_1k_metadata.rds
│       ├── pbmc_1k_cellData.txt.gz
│       └── pbmc_1k_featureData.txt.gz
├── pbmc_1k_QCParameters.yaml
├── Python
│   ├── Cells
│   └── Droplets
├── R
│   ├── pbmc_1k_Cells.rds
│   └── pbmc_1k_Droplets.rds
└── SCTK_pbmc_1k_cellQC_summary.csv
```
### partial results from custom sctk container on dataset 1 (4 out of 8)
* for comparison later with our multi step QC
* [merged_sctk_dataset1.txt](https://github.com/user-attachments/files/27086942/merged_sctk_dataset1.txt)
***
## CPU architecture specs
* https://camplab.net/sctk/current/articles/cmd_qc.html#running-sctk-qc-with-docker 
* in official documentation: `cpu_arch=broadwell|haswell|skylake|cascadelake`
* names may differ according to your HPC
* check via `sinfo -o "%N %f"`
* sample job script with CPU specs `#SBATCH --constraint="avx512"`
```bash
#!/bin/bash
#SBATCH --job-name=sctk_deftest
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --output=/global/scratch/%u/sctk_deftest_%j.out
#SBATCH --error=/global/scratch/%u/sctk_deftest_%j.err
#SBATCH --constraint="avx512"

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

module load apptainer/1.4.5

# specify tmpdir to slurm to prevent I/O bottlenecks, but cache can go into scratch
export APPTAINER_TMPDIR="$SLURM_TMPDIR"
export APPTAINER_CACHEDIR="/global/scratch/$USER/apptainer_cache"
export APPTAINERENV_TMPDIR="$SLURM_TMPDIR"

cd /global/scratch/$USER

apptainer exec \
--bind /global/scratch/$USER:/global/scratch/$USER \
/global/scratch/$USER/sctk_testdef_23april.sif \
Rscript /usr/local/bin/sctk_mito3_newnet_23april.R \
-P CellRangerV3 \
-R /global/scratch/$USER/test_sctk/raw_feature_bc_matrix \
-C /global/scratch/$USER/test_sctk/filtered_feature_bc_matrix \
-s pbmc_1k \
-o /global/scratch/$USER/test_sctk/pbmc_testdef \
-d Both \
-S TRUE \
-F SCE,Seurat,AnnData,FlatFile \
-M TRUE \
-E human-symbol \
-n "$SLURM_CPUS_PER_TASK"

echo -e "\n[ALERT] Run complete."
```



