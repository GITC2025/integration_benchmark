# Integration benchmark scRNA seq analysis

_figuring out the Frontenac HPC, babysitting the downloads, and benchmarking integration • HSCI591 W2026 BHSC, Queen's Uni_

* a first foray into HPC scRNA analysis 
* start date: mid March 2026 after lit review complete (Jan to March 2026) • work pace: intermittent ☕︎
* IDE: Positron • Environment: Frontenac Cluster HPC
* we're going thru each step in detail for learning purposes - the full scale integration will be containerised pipelines
* [**all shortcuts**](https://github.com/GITC2025/integration_benchmark/blob/main/shortcuts.md)
* [**sequencing protocols**](https://github.com/GITC2025/integration_benchmark/wiki/Protocols)

# Workflow
<img align="center" width="450" height="1300" alt="smallhsci591workflow" src="https://github.com/user-attachments/assets/1a8dbc5f-8b51-41f4-b46b-d6128e6b290b"/>
<br>when building the containerised pipelines: the pink ovals are logic gates to continue or pause the pipeline as appropriate.

## Methodology
```r
Part One: Data acquisition, alignment and QC
├──1. download fastq.gz from ENA database + md5 checksum
├──2. FastQC
├──3. Alignment: zcat check chemistry, Cellranger dataset 1 + Celescope dataset 2 
	├── run an alternative kallisto pseudoalignment to compare downstream analysis (optional)
├──4. QC with sctk vs Seurat QC)

Bridge
├──5. Normalization SCTransform (predictive modeling method)
├──6. Feature selection: ID HVGs (Seurat FindVariableFeatures etc.)

Part Two: Core Analysis
├──7. Dimension reduction
	├── PCA: scree plot to check PCs cover 70-90% of variance
	├── choose PCs for correction (cover 90-95% variance)

├──8. Integration benchmarking
	├── control vs scANVI vs Harmony etc. benchmark across diff categories of tools
	├── evaluate latent space with scIB

├──9. Clustering 
	├── Leiden clustering: adjust resolution 
	├── UMAP as qualitative benchmark
		├── set seeds, try diff k values for stability
	├── overlay Leiden cluster labels on UMAP
	├── ID outliers for final QC

├──10. Findmarkers Seurat to ID celltypes via gene enrichment (unsupervised)
	├── reference with database
```
# Dataset 1
Lai, H., Cheng, X., Liu, Q., Luo, W., Liu, M., Zhang, M., Miao, J., Ji, Z., Lin, G. N., Song, W., Zhang, L., Bo, J., Yang, G., Wang, J., & Gao, W. Q. (2021). Single-cell RNA sequencing reveals the epithelial cell heterogeneity and invasive subpopulation in human bladder cancer. International journal of cancer, 149(12), 2099–2115. https://doi.org/10.1002/ijc.33794

* https://www.ebi.ac.uk/ena/browser/view/SRP217277 
* 7 primary tumour samples + 1 normal sample = 8 SRA files
* for this analysis download fastq.gz as standard practice

# Dataset 2
Luo, Y., Tao, T., Tao, R., Huang, G., & Wu, S. (2022). Single-Cell Transcriptome Comparison of Bladder Cancer Reveals Its Ecosystem. Frontiers in oncology, 12, 818147. https://doi.org/10.3389/fonc.2022.818147
* https://www.ebi.ac.uk/ena/browser/view/SRP351272
* 4 samples total: 2x pri BCa, 1 recurrent BCa and one cystitis glandularis

## [download, md5 checksum scripts](https://github.com/GITC2025/integration_benchmark/blob/main/download_script.md) 

## Run fastqc on fastq.gz
* Reports per-base quality scores across read length (boxplots of Phred scores).
* Checks per-sequence quality (distribution of mean quality per read).
* Assesses GC content distribution vs expected
* Detects overrepresented sequences (adapters, primers, contaminants)
* Flags sequence duplication levels (PCR/optical duplicates)
* Examines per-base N content and other composition biases
* Summarizes results in an HTML report with pass/warn/fail flags for each module
* Recommend if trimming/adapter removal is needed before alignment.

**24 March fastqc**
* no falco on HPC
* unable to install falco on HPC - figure this out later

# fastqc
## dataset 1
```bash
#!/bin/bash
#SBATCH --job-name=fastqc_all
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err

set -euo pipefail

# fastQC html output goes here
cd /global/scratch/$USER/12marchENA_GSE
mkdir -p logs fastqc_out

module load fastqc

# run 8 fastQC parallel
find . -type f -name "*.fastq.gz" -print0 \
  | xargs -0 -n 1 -P "${SLURM_CPUS_PER_TASK:-8}" fastqc -t 1 -o fastqc_out

realpath fastqc_out
``` 
## dataset 2
```bash
#!/bin/bash
#SBATCH --job-name=fastqc_4samples
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/SRPfastqc_%j.out
#SBATCH --error=logs/SRPfastqc_%j.err

set -euo pipefail

cd /global/scratch/$USER/4samplesSRP351272
mkdir -p logs SRPfastqc_out

module load fastqc

# parallel 2 fastq, threads 2 per fastq
find . -type f -name "*.fastq.gz" -print0 \
  | xargs -0 -n 1 -P 2 fastqc -t 2 -o SRPfastqc_out

realpath SRPfastqc_out
```
## download FastQC html reports - from local console
```bash
scp -r YOUR_USERNAME@YOUR_SERVER_ADDRESS:/PATH/TO/REMOTE/FOLDER ./LOCAL_DESTINATION_FOLDER
```
## Fastqc_data.txt parsing Dataset 1

the txt are in the zip folders • [shortcuts here](https://github.com/GITC2025/integration_benchmark/blob/main/shortcuts.md#unzipfastq-shortcut)

```bash
# unzipfastq shortcut and move txt to current folder shortcut
unzipfastq

# shortcut checks
fastqcfail
fastqcwarn
```
* [fastqc parser results](https://github.com/GITC2025/integration_benchmark/blob/main/fastQCparser.md)
* nothing of concern here since no warn/fail for Per base sequence quality - most important module. 

## Reading FastQC output
<img width="650" height="85" alt="image" src="https://github.com/user-attachments/assets/cc260527-b231-4c35-a77b-4affb4f2dbc9" />

Illumina 1.9 = Phred +33

* Most important module
* There are 2 main types of Phred quality scores
* Phred +33 (Sanger/Illumina 1.8+) modern standard - most common
* Phred +64 (Illumina 1.3-1.7) older format
* FastQC will evaluate the Phred system used 
* Occasionally fastQC will guess wrongly - check the headers vs the fastQC data
* In general, quality may drop towards the end for fastq2 (actual sequence), this is normal
* fastq 1 (barcode/UMI) usually shows up as a flat line

<img width="450" height="370" alt="image" src="https://github.com/user-attachments/assets/0cb8210a-7631-42c9-bf2e-52f832f9d1dc" />

* 10x Genomics data often fails this module - not an issue
* fastq1 - fixed barcodes at beginning have specific sequences, creating a bias
	* UMI at the end of fastq1 have specific patterns 
* fastq2 - TSO at the begininning, creates bias
* random hexamers is biased towards specific sequences in the beginning

<img width="650" height="90" alt="image" src="https://github.com/user-attachments/assets/f9bd114b-a3c2-48f0-8dab-5bb0e6965220" />

* scRNA data often fails this module - not an issue
* fastQC was originally designed for genomic DNA, vs low-diversity libraries in scRNA-seq
* limited transcripts per cell create natural duplicates 
* deep sequencing - high coverage - duplicates
* Next step Cellranger will handle dedup

<img width="750" height="180" alt="image" src="https://github.com/user-attachments/assets/ac4afb9f-416d-4f57-be8a-a17dc5d06fd2" />

* these are usually TSOs - not an issue
* template switching oligos enable full length cDNA seq

<img width="500" height="380" alt="image" src="https://github.com/user-attachments/assets/bc6cb758-e4a2-41d3-8114-498e5bf72e9b" />

* frequent warns/fails here for read 2
* 3' read through - adaptor contam. increases towards 3'-end of the reads - sequencing continues beyond cDNA
* Cellranger will account for this via soft-clipping during STAR alignment
* additional trimming tools not advisable as may damage barcode/UMI

***
# Chemistry Detection with zcat 

* zcat verify R1 length manually 
* found R1 = 26bp = 16 bp barcode + 10 bp UMI ► v2 chemistry
* do not set auto to Cellranger - may be inaccurate

**v2 chemistry** (R1 = 26)

<img width="600" height="180" alt="image" src="https://github.com/user-attachments/assets/d5974ee4-e3ad-4062-af76-e0f05889e0ad" />

***

**v3,4 chemistry** (R1 = 28)

<img width="750" height="180" alt="image" src="https://github.com/user-attachments/assets/a730e6b5-8da7-4089-b6e1-18f5dd1cfab4" />

***

## Chemistry in brief
* v2 dominant prior to 2018 (16 + 10)
* post 2018 - v3 and 3.1 (16 + 12)
* 2023 onwards - v4 (16 + 12)
* v4 supports feature barcodes - run scRNA parallel w/ other assays
* multiplex for multi-omic output
* improved cell recovery 

***
## zcat Dataset 1 

[zcat shortcuts](https://github.com/GITC2025/integration_benchmark/blob/main/shortcuts.md#zcat-check-chemistry)

```bash
# zcat shortcuts
zcatR1

# output
SRR12539462_1.fastq.gz:
26
SRR12539463_1.fastq.gz:
26
SRR14615558_1.fastq.gz:
26
SRR9897621_1.fastq.gz:
26
SRR9897622_1.fastq.gz:
26
SRR9897623_1.fastq.gz:
26
SRR9897624_1.fastq.gz:
26
SRR9897625_1.fastq.gz:
26
```
## zcat Dataset 2

```bash
# use shortcut
zcatR1
zcatR2

# output
SRR17259462_1.fastq.gz:
150
SRR17259463_1.fastq.gz:
150
SRR17259464_1.fastq.gz:
150
SRR17259465_1.fastq.gz:
150

SRR17259462_2.fastq.gz:
150
SRR17259463_2.fastq.gz:
150
SRR17259464_2.fastq.gz:
150
SRR17259465_2.fastq.gz:
150
```

* Library construction: Singleron GEXSCOPER protocol (150x150 paired end chemistry)
* Cellranger cannot be used 
* use celescope 
* good stress test for cross-platform integration benchmarking later

# Alignment for dataset 1: Cellranger
* Cell Ranger 10.0.0 (standalone binaries)
* map to 10x GRCh38 2024-A reference genome
* actually cellranger v4 + don't need lane numbers but best practice include it
* since we don't have lane number info we use placeholder L001
* modern workflows merge lanes by default (Illumina `bcl2fastq` with the `--no-lane-splitting` option)

## dataset 1: rename fastq.gz files for cellranger via symlinks

```bash
# create new subfolder w/ symlink names formatted for cellranger

mkdir GSEfastqsymlinks

for f in SRR*_1.fastq.gz; do
  base="${f%_1.fastq.gz}"
  ln -s -- "$PWD/$f" "GSEfastqsymlinks/${base}_S1_L001_R1_001.fastq.gz"
done

for f in SRR*_2.fastq.gz; do
  base="${f%_2.fastq.gz}"
  ln -s -- "$PWD/$f" "GSEfastqsymlinks/${base}_S1_L001_R2_001.fastq.gz"
done
```

dataset 1 (first 4 files)
<img width="900" height="130" alt="image" src="https://github.com/user-attachments/assets/c7864647-af49-4375-ac44-af3fe0e2487e" />

## prep for Cellranger

1. get the latest human reference genome GRCh38
2. generate a samples.txt file as input for cellranger
3. run cellranger as an array: 8 arrays for 8 distinct biological samples (8 pairs of fastq)
* Since we're only running cellranger on one dataset of 8 paired fastq.gz files, we'll use the scratch space
* larger datasets require more careful space management to ensure you hv sufficient space during cellranger count
* Each 25GB fastq R2 can generate 500-700GB intermediate files which are deleted, and 200-280GB final output files
* 8 arrays: Each distinct biological sample requires its own cellranger count job
* In our dataset, each fastq pair = one patient = one distinct biological sample = one well
* There was no lane number in the original fastq - cellranger auto detects lanes/ merged lanes (modern standard) (Lanes are from PCR seq)
* If this is the older style of split lanes - you pass all lanes from one sample/well to cellranger for one cellranger count job

Note: In other datasets, one sample could be one GEM well of multiple patient samples (multiplexed), or one patient sample could be split into multiple GEM wells. Check the context.

```bash
# in our scratch space make a working directory
mkdir run_cellranger_count_GSE

# Download the latest precompiled human reference genome and unpack it - just takes a min on HPC
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

# create samples.txt in our fastq.gz folder then move samples.txt to cellranger working directory
ls | grep -o 'SRR[0-9]\+' | sort -u > samples.txt

# check output
cat samples.txt

SRR12539462
SRR12539463
SRR14615558
SRR9897621
SRR9897622
SRR9897623
SRR9897624
SRR9897625
```
```bash
# check reference genome GRCh38
[hpc6140@frnt150 run_cellranger_count_GSE]$ ls -lh refdata-gex-GRCh38-2024-A
total 2.0K
drwxr-x---. 2 hpc6140 hpc6140 4.0K Mar  8  2024 fasta
drwxr-x---. 2 hpc6140 hpc6140 4.0K Mar  8  2024 genes
-rw-r-----. 1 hpc6140 hpc6140  451 Mar  8  2024 reference.json
drwxr-x---. 2 hpc6140 hpc6140 4.0K Mar  8  2024 star
```
* fasta contains the genome sequence
* gene contains annotations
* star has pre computed STAR index files - the core engine [spliced transcripts alignment to reference]

## run the cellranger array
```bash
#!/bin/bash
#SBATCH --job-name=GSE_array
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --output=/global/scratch/hpc6140/run_cellranger_count_GSE/logs2/job_%A_%a.out
#SBATCH --error=/global/scratch/hpc6140/run_cellranger_count_GSE/logs2/error_%A_%a.err

# set our soft limits to system hard limits to optimise
ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

# verify limits in output log
echo "Environment Setup"
echo "Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Max open files (-n): $(ulimit -n)"
echo "Max user processes (-u): $(ulimit -u)"

mkdir -p /global/scratch/hpc6140/run_cellranger_count_GSE/logs2

export PATH=/global/project/hpcg6049/software/cellranger-10.0.0/bin:$PATH

# our working directory containing sample.txt and where output goes
cd /global/scratch/hpc6140/run_cellranger_count_GSE

# the fastqsymlinks are formatted for cellranger
FASTQ_DIR="/global/scratch/hpc6140/12marchENA_GSE/GSEfastqsymlinks"
SAMPLE_LIST="samples.txt"

# each unique sample ID (fastq pair) gets one array for 8 arrays total
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

cellranger count --id=${SAMPLE_ID}_Aligned \
--create-bam=true \
--transcriptome=/global/scratch/hpc6140/run_cellranger_count_GSE/refdata-gex-GRCh38-2024-A \
--fastqs=${FASTQ_DIR} \
--sample=${SAMPLE_ID} \
--include-introns=true \
--localcores=$SLURM_CPUS_PER_TASK \
--localmem=90
```
* last line localmem must be set lower than requested RAM (about 4-8GB lower) to prevent crashes
* bash/linux HPC systems do not requrie setting localvmem - may lead to crashes due to conflict
* cellranger v8+ requires create-bam=true/false (will fail immed if you don't specify)
* Takes a while of tweaking the specs and timing to ensure cellranger doesn't crash or run out of time, yet balance against reasonable job queueing time.
* In case cellranger times out, you can resume the job by submitting another one. (kallisto does not have this feature)
* This takes a day or two of waiting for resources and running the jobs. 
* we run a [kallisto alignment](https://github.com/GITC2025/integration_benchmark/blob/main/kallisto_analysis.md) as an alternative when cellranger is done

```bash
# head output in job logs - array 1 has started after about 24h of waiting
# our dynamic ulimits set up nicely, beyond cellranger's minimum needs
Environment Setup
Job ID: 7793386
Task ID: 1
Max open files (-n): 131072
Max user processes (-u): 2059486

# several arrays have completed, other arrays are ongoing - completion output for one array 
# start time
2026-04-12 09:27:39 [runtime] (ready)           ID.SRR12539462_Aligned.SC_RNA_COUNTER_CS.SC_MULTI_CORE.SANITIZE_MAP_CALLS

# end time 
2026-04-12 19:37:02 [runtime] (join_complete)   ID.SRR12539462_Aligned.SC_RNA_COUNTER_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.COUNT_GEM_WELL_PROCESSOR._BASIC_SC_RNA_COUNTER._POST_MATRIX_COMPUTATION.WRITE_POS_BAM

Outputs:
- Run summary HTML:                         /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/web_summary.html     
- Run summary CSV:                          /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/metrics_summary.csv  
- BAM:                                      /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/possorted_genome_bam.bam
- BAM BAI index:                            /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/possorted_genome_bam.bam.bai
- BAM CSI index:                            null
- Filtered feature-barcode matrices MEX:    /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/analysis
- Per-molecule read information:            /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/molecule_info.h5     
- CRISPR-specific analysis:                 null
- Antibody aggregate barcodes:              null
- Loupe Browser file:                       /global/scratch/hpc6140/run_cellranger_count_GSE/SRR12539462_Aligned/outs/cloupe.cloupe        
- Feature Reference:                        null
- Target Panel File:                        null
- Probe Set File:                           null
- cell_types:                               null
- Symlink web_summary_cell_types.html:      null

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

# elapsed time for all arrays to optimise future job time settings
[hpc6140@login1 logs2]$ sacct -j 7793386 -X --format=JobID,Elapsed
JobID           Elapsed
------------ ----------
7793386_1      10:10:11
7793386_2      04:04:51
7793386_3      09:50:44
7793386_4      03:52:25
7793386_5      05:21:57
7793386_6      22:48:10
7793386_7      09:26:23
7793386_8      07:00:41
```
## assess cellranger metrics quality
[SRP217277_cellranger_reports.zip](https://github.com/user-attachments/files/26834408/SRP217277_cellranger_reports.zip)

if you only have a few files, compile all web_summary.html into a staging folder with the SRR prefixes
```bash
cd /global/scratch/hpc6140/run_cellranger_count_GSE
mkdir -p staging_htmls
find . -mindepth 2 -type f -name 'web_summary.html' | while read file; do
prefix=$(echo "$file" | cut -d'/' -f2)
cp "$file" "staging_htmls/${prefix}_web_summary.html"
echo "Staged: ${prefix}_web_summary.html"
done
```
from your local console pull the staging folder in the HPC into your local folder
```bash
scp hpc6140@frontenac.cac.queensu.ca:/global/scratch/hpc6140/run_cellranger_count_GSE/staging_htmls/* "C:\Users\Finnt\Downloads\"
```
<img width="1669" height="1065" alt="CRreport" src="https://github.com/user-attachments/assets/05f4ba4c-e035-45e8-8048-a60976d2d46a" />
<br>
* the junction at the knee plot separates genuine reads from background via UMI count vs barcode ranking (transcript abundance vs cells)

## explanations 
**1. Sequencing Quality**
* **Valid Barcodes:** > 85%, lower indicates chemistry failure or incorrect barcode whitelist.
* **Q30 Bases in Barcode / RNA Read / UMI:** > 85%, Measures base calling accuracy; lower values indicate poor sequencing run quality 

**2. Mapping Efficiency**
* **Reads Mapped Confidently to Genome:** > 85%, lower values suggest contamination or poor sample quality → run background remover
* **Reads Mapped Confidently to Transcriptome:** scRNA-seq: > 60%
    * snRNA-seq: > 30-50% - Nuclei have high % unspliced intronic RNA - maps to genome but not confidently to mature transcriptome, unless reference include pre-mRNA (e.g. in cellranger, use --include-introns flag, similar to using the human nac index for kallisto)

**3. Cell Capture & Background**
* **Estimated Number of Cells:** Should align with your targeted recovery if you know (e.g., within `± 20%` of loaded cell count).
* **Fraction Reads in Cells:** scRNA-seq: > 70%		• 	snRNA-seq: 40% - 60%  • lower values need background removal
	* _if you see a sharp drop in the knee plot, and have a fraction reads in cells > 70% (meaning 30% is background), you generally don't need to run a background remover. However it may be a good idea to remove background if you have low abundance cell types. when your fraction reads is very high >90%, bg removal may lead to over correction, so it's not needed._ 

**4. Library Complexity & Depth**
* **Mean Reads per Cell:** > 20,000 
* **Median UMI Counts per Cell:** tissue-dependent. Generally > 1000 (scRNA) or > 500 (snRNA).
* **Median Genes per Cell:** tissue-dependent. 500 - 4,000+`. (Drops significantly in metabolically inactive cells or snRNA assays).
* **Sequencing Saturation:** > 50%

## parsing for key metrics

if you have a lot of files you can parse the metrics_summary.csv directly
* shortcut logic gate [`CRmetrics()`](https://github.com/GITC2025/integration_benchmark/blob/main/shortcuts.md#cellranger-qc-parser-for-scrna-data) customise the QC parameters to what you want to know, adjust QC thresholds for snRNA data (under construction)
* [MultiQC](https://github.com/MultiQC/MultiQC) if you need a polished aggregated html report

## logic gate: check all cellranger arrays completed
* shortcut parser [`pipesuccess`](https://github.com/GITC2025/integration_benchmark/blob/main/shortcuts.md#cellranger-pipestance-completion-parser)
* parses tail -n 5 of all .out logs for Pipestance completed successfully!
* prints log numbers without success statement
* if time out error- resend the job for that SRR and cellranger would resume where it left off

***
# Align Dataset 2 with celescope 

| Feature | 10x Genomics (Chromium) | Singleron (GEXSCOPE) |
| :--- | :--- | :--- |
| **Technology** | Droplet-based | Microwell-based |
| **Throughput** | High (Up to >2M cells) | Medium-High |
| **Workflow** | Automated (Instrument) | Automated/Manual |
| **Best For** | Massive, complex tissue atlases | Clinical samples, small biopsies |
| **Sensitivity** | Very high (low dropout) | High |
| **Fastq format** | R1 26 (v2) or 28 (v3,4) | R1, R2 both 150 each |

## set up celescope for Dataset 2
use the apptainer for the standalone unit with everything included
* bare metal celescope is a distributed bundle unlike cellranger standalone binaries
* w/o apptainer: python interpreter is system dependent, binary dependencies need to be pathed
```bash
apptainer pull celescope_v2.11.4.sif docker://quay.io/singleron-rd/celescope:v2.11.4
```
remake existing reference genome to be compatible for celescope - make another copy first just for celescope
```bash
# first unzip your genes.gtf.gz if its gz - if downloaded from 10X precompiled (will throw error if zipped)
cd /global/scratch/hpc6140/refdata-gex-GRCh38-2024-A/genes/
gunzip genes.gtf.gz
```
```bash
#!/bin/bash
#SBATCH --job-name=celescope_mkref
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --output=/global/scratch/hpc6140/run_celescope_mkref_%j.out
#SBATCH --error=/global/scratch/hpc6140/run_celescope_mkref_%j.err

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

module load apptainer/1.4.5

mkdir -p /global/scratch/hpc6140/celescope_GRCh38_ref_output
cd /global/scratch/hpc6140/celescope_GRCh38_ref_output

apptainer exec --bind /global/scratch/hpc6140:/global/scratch/hpc6140 /global/scratch/hpc6140/celescope_v2.11.4.sif celescope rna mkref \
--genome_name celescope_GRCh38_ref \
--fasta /global/scratch/hpc6140/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
--gtf /global/scratch/hpc6140/refdata-gex-GRCh38-2024-A/genes/genes.gtf \
--thread $SLURM_CPUS_PER_TASK
```
```bash
# output from logs
tail -1 run_celescope_mkref_8018948.out
Apr 13 23:58:56 ..... finished successfully

tail -1 run_celescope_mkref_8018948.err
2026-04-13 23:58:56,589 - celescope.rna.mkref.run - INFO - done. time used: 0:17:45.353121
```
this is how celescope rna mkref works

| Argument 				| Description |
| :--- | :--- |
| `genome_name` | **Required.** The name of the output directory where the reference will be built. |
| `fasta` | **Required.** Path to the genome FASTA file. |
| `gtf` | **Required.** Path to the genome GTF annotation file. |
| `thread` | Default is 6; 16 is recommended for GRCh38. soft code as $SLURM_CPUS_PER_TASK |
| `mt_gene_list` 	| **Optional.** Path to a file containing a list of mito genes (QC metrics in report). we're doing mito content QC downstream so we have skipped mt_gene_list here.  |

this is what you should get when the ref genome is made
```bash
[hpc6140@frnt128 celescope_GRCh38_ref]$ tree
.
├── celescope_genome.config
├── chrLength.txt
├── chrNameLength.txt
├── chrName.txt
├── chrStart.txt
├── exonGeTrInfo.tab
├── exonInfo.tab
├── geneInfo.tab
├── Genome
├── genomeParameters.txt
├── Log.out
├── SA
├── SAindex
├── sjdbInfo.txt
├── sjdbList.fromGTF.out.tab
├── sjdbList.out.tab
└── transcriptInfo.tab
```
**Most important components**

| File Name | Description |
| :--- | :--- |
| `celescope_genome.config` | config file for celescope |
| `Genome` | ref genome in fast-access binary format. |
| `SA` | suffix array - the search engine |

## prep a sample file for celescope

```bash
# print all unique SRRs to mapfile.tsv
ls *.fastq* | cut -d'_' -f1 | sort -u | while read id; do
printf "${id}\t${PWD}\t${id}\n"
done > mapfile.tsv

# output
cat mapfile.tsv

SRR17259462     /global/scratch/hpc6140/4samplesSRP351272       SRR17259462
SRR17259463     /global/scratch/hpc6140/4samplesSRP351272       SRR17259463
SRR17259464     /global/scratch/hpc6140/4samplesSRP351272       SRR17259464
SRR17259465     /global/scratch/hpc6140/4samplesSRP351272       SRR17259465
```
after 2 days, our cellranger last array is wrapping up. 

## run celescope array
```bash
#!/bin/bash
#SBATCH --job-name=celescope_align
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --array=1-4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --output=/global/scratch/hpc6140/logs/celescope_%A_%a.out
#SBATCH --error=/global/scratch/hpc6140/logs/celescope_%A_%a.err

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)
module load apptainer/1.4.5

SIF="/global/scratch/hpc6140/celescope_v2.11.4.sif"
REF="/global/scratch/hpc6140/celescope_GRCh38_ref"
FQ_DIR="/global/scratch/hpc6140/4samplesSRP351272"
BASE_OUTDIR="/global/scratch/hpc6140/celescope_output"

SRRS=($(ls ${FQ_DIR}/*_1.fastq.gz | xargs -n1 basename | sed 's/_1.fastq.gz//'))
SRR_ID=${SRRS[$((SLURM_ARRAY_TASK_ID - 1))]}

SAMPLE_NAME=$SRR_ID
SAMPLE_OUTDIR="${BASE_OUTDIR}/${SAMPLE_NAME}"

mkdir -p "$SAMPLE_OUTDIR"
cd "$SAMPLE_OUTDIR"

apptainer exec --bind /global/scratch/hpc6140:/global/scratch/hpc6140 "$SIF" \
celescope rna starsolo \
--fq1 "${FQ_DIR}/${SRR_ID}_1.fastq.gz" \
--fq2 "${FQ_DIR}/${SRR_ID}_2.fastq.gz" \
--genomeDir "$REF" \
--sample "$SAMPLE_NAME" \
--thread "$SLURM_CPUS_PER_TASK" \
--chemistry auto \
--outdir "$SAMPLE_OUTDIR"
```
## celescope complete in 1.5h per array
```bash
        STAR version: 2.7.11a   compiled: 2023-09-15T02:58:53+0000 :/opt/conda/conda-bld/star_1694746407721/work/source
Apr 14 21:58:20 ..... started STAR run
Apr 14 21:58:20 ..... loading genome
Apr 14 21:58:47 ..... started mapping
Apr 14 22:57:37 ..... finished mapping
Apr 14 22:57:39 ..... started Solo counting
Apr 14 23:08:23 ..... finished Solo counting
Apr 14 23:08:23 ..... started sorting BAM
Apr 14 23:20:30 ..... finished successfully
```
## Singleron barcodes
* combinatorial barcode synthesis on the bead
* $T = N^3$ : exponential increase in combinations with each round 	
* (compare w/ Parse Biosciences that combi barcode inside the cell)
* n rounds of barcoding = n barcodes in a row sep by linkers

```bash
# from official docs
GEXSCOPE Microbead (single barcode): C12U8
2018 GEXSCOPE V1: C8L16C8L16C8L1U12
2020 GEXSCOPE V2: C9L16C9L16C9L1U12 (Linker 16bp)
2024 GEXSCOPE V3: C9L6C9L6C9L1U12 (Linker 6bp)
```
* barcode structures more complex than cellranger
* many GEXSCOPE R1 are 150bp, so you can't just check the # bp in R1 to detect chemistry quickly
* the dataset is from 2022, so it could be microbead, V1 or V2 (not specified in paper)
* for manual verification: check barcode structures here
* official doc: [chemistry](https://github.com/singleron-RD/CeleScope/blob/master/doc/chemistry.md) • [chemistry dictionary](https://github.com/singleron-RD/CeleScope/blob/master/celescope/chemistry_dict.py)

## verify chemistry auto-detection
[view celescope QC reports](https://github.com/GITC2025/integration_benchmark/blob/main/SRP351272_celescope_reports.pdf)

* relatively high valid reads and low corrected barcodes means chem detection was correct
* the out logs also show the results of the barcode subsampling - there is majority vote for GEXSCOPE V1
* the minority votes for other chemistry are seq. errors detected as other chem, not of concern

```bash
[hpc6140@frnt153 logs]$ tail -n 3 *.out
==> celescope_8024838_1.out <==
Apr 14 23:20:30 ..... finished successfully
[('GEXSCOPE-V1', 7904), ('flv_rna', 98), ('GEXSCOPE-MicroBead', 5)]
/global/scratch/hpc6140/4samplesSRP351272/SRR17259462_1.fastq.gz: GEXSCOPE-V1

==> celescope_8024838_2.out <==
Apr 14 23:48:36 ..... finished successfully
[('GEXSCOPE-V1', 8176), ('flv_rna', 67), ('GEXSCOPE-MicroBead', 4)]
/global/scratch/hpc6140/4samplesSRP351272/SRR17259463_1.fastq.gz: GEXSCOPE-V1

==> celescope_8024838_3.out <==
Apr 14 23:50:13 ..... finished successfully
[('GEXSCOPE-V1', 8291), ('flv_rna', 56), ('GEXSCOPE-MicroBead', 1)]
/global/scratch/hpc6140/4samplesSRP351272/SRR17259464_1.fastq.gz: GEXSCOPE-V1

==> celescope_8024838_4.out <==
Apr 14 23:35:27 ..... finished successfully
[('GEXSCOPE-V1', 8314), ('flv_rna', 77), ('GEXSCOPE-MicroBead', 1)]
/global/scratch/hpc6140/4samplesSRP351272/SRR17259465_1.fastq.gz: GEXSCOPE-V1
```
## read QC reports
* the UMI elbow plots look fine for differentiating background vs real UMI counts
* of slight concern is the borderline scores for reads mapped to unique loci, and relatively high reads mapped to multiple loci
* we should investigate this further 
***
# PCA in a nutshell

<img width="800" height="440" alt="PCAsteps" src="https://github.com/user-attachments/assets/adcc1405-3f5c-4d00-b761-0593a56f8734" />

***
