
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
scp -r YOUR_$USERNAME@YOUR_SERVER_ADDRESS:/PATH/TO/REMOTE/FOLDER ./LOCAL_DESTINATION_FOLDER
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
