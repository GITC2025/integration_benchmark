# Alignment for dataset 1: Cellranger
* [Dataset1 Cellranger html reports](https://github.com/GITC2025/integration_benchmark/blob/main/SRP217277_cellranger_reports.zip)
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
[$USER@frnt150 run_cellranger_count_GSE]$ ls -lh refdata-gex-GRCh38-2024-A
total 2.0K
drwxr-x---. 2 $USER $USER 4.0K Mar  8  2024 fasta
drwxr-x---. 2 $USER $USER 4.0K Mar  8  2024 genes
-rw-r-----. 1 $USER $USER  451 Mar  8  2024 reference.json
drwxr-x---. 2 $USER $USER 4.0K Mar  8  2024 star
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
#SBATCH --mail-user=<EMAIL>
#SBATCH --output=/global/scratch/%u/run_cellranger_count_GSE/logs2/job_%A_%a.out
#SBATCH --error=/global/scratch/%u/run_cellranger_count_GSE/logs2/error_%A_%a.err

# set our soft limits to system hard limits to optimise
ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

# verify limits in output log
echo "Environment Setup"
echo "Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Max open files (-n): $(ulimit -n)"
echo "Max $USER processes (-u): $(ulimit -u)"

mkdir -p /global/scratch/$USER/run_cellranger_count_GSE/logs2

export PATH=/global/project/hpcg6049/software/cellranger-10.0.0/bin:$PATH

# our working directory containing sample.txt and where output goes
cd /global/scratch/$USER/run_cellranger_count_GSE

# the fastqsymlinks are formatted for cellranger
FASTQ_DIR="/global/scratch/$USER/12marchENA_GSE/GSEfastqsymlinks"
SAMPLE_LIST="samples.txt"

# each unique sample ID (fastq pair) gets one array for 8 arrays total
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

cellranger count --id=${SAMPLE_ID}_Aligned \
--create-bam=true \
--transcriptome=/global/scratch/$USER/run_cellranger_count_GSE/refdata-gex-GRCh38-2024-A \
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
Max $USER processes (-u): 2059486

# several arrays have completed, other arrays are ongoing - completion output for one array 
# start time
2026-04-12 09:27:39 [runtime] (ready)           ID.SRR12539462_Aligned.SC_RNA_COUNTER_CS.SC_MULTI_CORE.SANITIZE_MAP_CALLS

# end time 
2026-04-12 19:37:02 [runtime] (join_complete)   ID.SRR12539462_Aligned.SC_RNA_COUNTER_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.COUNT_GEM_WELL_PROCESSOR._BASIC_SC_RNA_COUNTER._POST_MATRIX_COMPUTATION.WRITE_POS_BAM

Outputs:
- Run summary HTML:                         /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/web_summary.html     
- Run summary CSV:                          /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/metrics_summary.csv  
- BAM:                                      /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/possorted_genome_bam.bam
- BAM BAI index:                            /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/possorted_genome_bam.bam.bai
- BAM CSI index:                            null
- Filtered feature-barcode matrices MEX:    /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/analysis
- Per-molecule read information:            /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/molecule_info.h5     
- CRISPR-specific analysis:                 null
- Antibody aggregate barcodes:              null
- Loupe Browser file:                       /global/scratch/$USER/run_cellranger_count_GSE/SRR12539462_Aligned/outs/cloupe.cloupe        
- Feature Reference:                        null
- Target Panel File:                        null
- Probe Set File:                           null
- cell_types:                               null
- Symlink web_summary_cell_types.html:      null

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

# elapsed time for all arrays to optimise future job time settings
[$USER@login1 logs2]$ sacct -j 7793386 -X --format=JobID,Elapsed
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

if you only have a few files, compile all web_summary.html into a staging folder with the SRR prefixes
```bash
cd /global/scratch/$USER/run_cellranger_count_GSE
mkdir -p staging_htmls
find . -mindepth 2 -type f -name 'web_summary.html' | while read file; do
prefix=$(echo "$file" | cut -d'/' -f2)
cp "$file" "staging_htmls/${prefix}_web_summary.html"
echo "Staged: ${prefix}_web_summary.html"
done
```
from your local console pull the staging folder in the HPC into your local folder
```bash
scp $USER@frontenac.cac.queensu.ca:/global/scratch/$USER/run_cellranger_count_GSE/staging_htmls/* "C:\$USERs\Finnt\Downloads\"
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
