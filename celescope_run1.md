# celescope align first run 

**resulted in low mapping rates - rerun with pri assembly GRCH38**

[celescope run1 reports](https://github.com/GITC2025/integration_benchmark/blob/main/SRP351272_run1_celescope_reports.pdf)

remake existing reference genome to be compatible for celescope - make another copy first just for celescope
```bash
# first unzip your genes.gtf.gz if its gz - if downloaded from 10X precompiled (will throw error if zipped)
cd /global/scratch/$USER/refdata-gex-GRCh38-2024-A/genes/
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
#SBATCH --output=/global/scratch/$USER/run_celescope_mkref_%j.out
#SBATCH --error=/global/scratch/$USER/run_celescope_mkref_%j.err

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

module load apptainer/1.4.5

mkdir -p /global/scratch/$USER/celescope_GRCh38_ref_output
cd /global/scratch/$USER/celescope_GRCh38_ref_output

apptainer exec --bind /global/scratch/$USER:/global/scratch/$USER /global/scratch/$USER/celescope_v2.11.4.sif celescope rna mkref \
--genome_name celescope_GRCh38_ref \
--fasta /global/scratch/$USER/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
--gtf /global/scratch/$USER/refdata-gex-GRCh38-2024-A/genes/genes.gtf \
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
[$USER@frnt128 celescope_GRCh38_ref]$ tree
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

SRR17259462     /global/scratch/$USER/4samplesSRP351272       SRR17259462
SRR17259463     /global/scratch/$USER/4samplesSRP351272       SRR17259463
SRR17259464     /global/scratch/$USER/4samplesSRP351272       SRR17259464
SRR17259465     /global/scratch/$USER/4samplesSRP351272       SRR17259465
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
#SBATCH --output=/global/scratch/$USER/logs/celescope_%A_%a.out
#SBATCH --error=/global/scratch/$USER/logs/celescope_%A_%a.err

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)
module load apptainer/1.4.5

SIF="/global/scratch/$USER/celescope_v2.11.4.sif"
REF="/global/scratch/$USER/celescope_GRCh38_ref"
FQ_DIR="/global/scratch/$USER/4samplesSRP351272"
BASE_OUTDIR="/global/scratch/$USER/celescope_output"

SRRS=($(ls ${FQ_DIR}/*_1.fastq.gz | xargs -n1 basename | sed 's/_1.fastq.gz//'))
SRR_ID=${SRRS[$((SLURM_ARRAY_TASK_ID - 1))]}

SAMPLE_NAME=$SRR_ID
SAMPLE_OUTDIR="${BASE_OUTDIR}/${SAMPLE_NAME}"

mkdir -p "$SAMPLE_OUTDIR"
cd "$SAMPLE_OUTDIR"

apptainer exec --bind /global/scratch/$USER:/global/scratch/$USER "$SIF" \
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
[$USER@frnt153 logs]$ tail -n 3 *.out
==> celescope_8024838_1.out <==
Apr 14 23:20:30 ..... finished successfully
[('GEXSCOPE-V1', 7904), ('flv_rna', 98), ('GEXSCOPE-MicroBead', 5)]
/global/scratch/$USER/4samplesSRP351272/SRR17259462_1.fastq.gz: GEXSCOPE-V1

==> celescope_8024838_2.out <==
Apr 14 23:48:36 ..... finished successfully
[('GEXSCOPE-V1', 8176), ('flv_rna', 67), ('GEXSCOPE-MicroBead', 4)]
/global/scratch/$USER/4samplesSRP351272/SRR17259463_1.fastq.gz: GEXSCOPE-V1

==> celescope_8024838_3.out <==
Apr 14 23:50:13 ..... finished successfully
[('GEXSCOPE-V1', 8291), ('flv_rna', 56), ('GEXSCOPE-MicroBead', 1)]
/global/scratch/$USER/4samplesSRP351272/SRR17259464_1.fastq.gz: GEXSCOPE-V1

==> celescope_8024838_4.out <==
Apr 14 23:35:27 ..... finished successfully
[('GEXSCOPE-V1', 8314), ('flv_rna', 77), ('GEXSCOPE-MicroBead', 1)]
/global/scratch/$USER/4samplesSRP351272/SRR17259465_1.fastq.gz: GEXSCOPE-V1
```
## read QC reports
* the UMI elbow plots look fine for differentiating background vs real UMI counts
* of slight concern is the borderline scores for reads mapped to unique loci, and relatively high reads mapped to multiple loci
* we should investigate this further 
***
## dataset 2 celescope QC metrics analysis

**1. Sequencing Quality**
* **Valid Barcodes:** > 85%
  * borderline passes - range 84.6% to 89%
* **Q30 Bases in Barcode / RNA Read / UMI:** > 85%
  * clear pass in all criteria, ranges from 86.6% to 96%

**2. Mapping Efficiency**
* **Reads Mapped to Unique Loci = ( cellranger Reads Mapped Confidently to Genome):** > 85%
  * 43.3% to 60.9% - troubleshoot this
* **Reads Mapped Uniquely/Confidently to Transcriptome:** scRNA-seq: > 60%
  * 36.3% to 56.6% - very low
* generally much higher genomic mapping vs transcriptomic mapping - could indicate annotation mismatch  
* seq quality module pass indicates that low quality mapping is not due to incorrect chemistry

**3. Cell Capture & Background**
* **Fraction Reads in Cells:** scRNA-seq: > 70%
  * 66.7 % to 75% - borderline (mainly chemistry and sample quality issue)

**4. Library Complexity & Depth**
* **Mean Reads per Cell:** > 20,000 
  * clear pass - but many reads discarded, so used reads were low
* **Median UMI Counts per Cell:** tissue-dependent. Generally > 1000 (scRNA) or > 500 (snRNA).
  * clear pass
* **Median Genes per Cell:** tissue-dependent. 500 - 4,000+`
  * ranges from 1341 to 2770
* **Sequencing Saturation:** > 50%
  * borderline pass

to diagnose the reasons - in the celescope output folder find the logs

```bash
[$USER@frnt153 celescope_output]$ tree

├── SRR17259462
│   ├── SRR17259462_Log.final.out
│   ├── SRR17259462_Log.out
│   ├── SRR17259462_Log.progress.out
│   ├── SRR17259462_SJ.out.tab

tail -n 40 SRR17259462_Log.final.out

# partial stats
                                    UNIQUE READS:
                   Uniquely mapped reads number |       220754281
                        Uniquely mapped reads % |       54.68%
                          Average mapped length |       138.18
                       Number of splices: Total |       69618509
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       82003592
             % of reads mapped to multiple loci |       20.31%
        Number of reads mapped to too many loci |       486668
             % of reads mapped to too many loci |       0.12%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       67650560
                 % of reads unmapped: too short |       16.76%
                Number of reads unmapped: other |       32806376
                     % of reads unmapped: other |       8.13%
```

**next steps - rerun alignment with some adjustments to see if things improve**

* much higher genomic mapping vs transcriptomic mapping could indicate annotation mismatch
* we previously built from 10x precompiled GRch38 human ref genome
   * may cause attribute name mismatch, loss of non standard genes, alt haplotypes multi mapping 
* rebuild ref genome using raw primary assembly 
  * could improve mapping, full dictionary and compatibility
* add include-introns flag to account for pre-mRNA in scRNA data
  * reduces discarding of pre-mRNA reads
