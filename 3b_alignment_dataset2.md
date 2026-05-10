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
* [celescope mkref and align run 1](https://github.com/GITC2025/integration_benchmark/blob/main/celescope_run1.md)
* [celescope run1 reports](https://github.com/GITC2025/integration_benchmark/blob/main/3b1_SRP351272_run1_celescope_reports.pdf)
* low mappping rates, rerun with original primary assembly GRch38

## celescope run 2
### make reference genome for celescope
**download the Enseml release 115**

* use primary assembly paired with pri. assembly GTF 
  * for standard aligners like STAR, HISAT2, Bowtie2
  * prevents multi mapping issues
* p.a. = single, non redundant representation
* contains all chromosomes: 1-22, X, Y and mito + scaffolds/patches
* contrast with Toplevel (full) assembly - contains redundancies
  * contains chroms, scaffold, haplotypes and patches
  * more complete view of human genetic diversity
  * for high resolution variant calling
  * requires alt-aware mappers like BWA MEM or minimap2
  * may results in multi-mapping for typical usage

```bash
cd /global/scratch/$USER
wget http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
```
* we have to process the GTF gene transfer format for our usage - use `celescope utils mkgtf`
  * otherwise using the full GTF results in high multi mapping rates 
* raw GTF contains: non standard chromosomes such as scaffolds and patches 
  * pseudogenes, various types of ncRNA, special categories like TEC (to be experimentally confirmed), artifacts
* use recommended GTF filtering settings: protein coding, lncRNA and immune genes

```bash
# to see the gene biotypes in your original gtf
cut -f9 Homo_sapiens.GRCh38.115.gtf | grep -o 'gene_biotype "[^"]*"' | sort | uniq
```
[full GTF categories here](https://github.com/GITC2025/integration_benchmark/blob/main/GTF_gene_biotypes_115.md)

Filter the GTF for scRNA data
```bash
module load apptainer/1.4.5

SIF="/global/scratch/$USER/celescope_v2.11.4.sif"

apptainer exec --bind /global/scratch/$USER:/global/scratch/$USER "$SIF" \
celescope utils mkgtf \
--attributes "gene_biotype=protein_coding,lncRNA,IG_V_gene,IG_J_gene,TR_V_gene,TR_D_gene,TR_J_gene,TR_C_gene;" \
Homo_sapiens.GRCh38.115.gtf \
Homo_sapiens.GRCh38.115.filtered.gtf
```
check the filtered GTF
```bash
cut -f9 Homo_sapiens.GRCh38.115.filtered.gtf | grep -o 'gene_biotype "[^"]*"' | sort | uniq
```
```bash
gene_biotype "IG_J_gene"
gene_biotype "IG_V_gene"
gene_biotype "lncRNA"
gene_biotype "protein_coding"
gene_biotype "TR_C_gene"
gene_biotype "TR_D_gene"
gene_biotype "TR_J_gene"
gene_biotype "TR_V_gene"
```
make the reference for celescope now, using our filtered GTF and primary assembly 115
```bash
#!/bin/bash
#SBATCH --job-name=celescope_mkref
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<EMAIL>
#SBATCH --output=/global/scratch/%u/19april_celescope_mkref_%j.out
#SBATCH --error=/global/scratch/%u/19april_celescope_mkref_%j.err

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

module load apptainer/1.4.5

mkdir -p /global/scratch/$USER/celescope_GRCh38_ref_output_19april
cd /global/scratch/$USER/celescope_GRCh38_ref_output_19april

apptainer exec --bind /global/scratch/$USER:/global/scratch/$USER /global/scratch/$USER/celescope_v2.11.4.sif celescope rna mkref \
--genome_name celescope_GRCh38_ref_19april \
--fasta /global/scratch/$USER/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /global/scratch/$USER/Homo_sapiens.GRCh38.115.filtered.gtf \
--thread $SLURM_CPUS_PER_TASK
```
rerun celescope align to our new reference
```bash
#!/bin/bash
#SBATCH --job-name=celescope_align
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --array=1-4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<EMAIL>
#SBATCH --output=/global/scratch/%u/logs/19aprilcelescope_%A_%a.out
#SBATCH --error=/global/scratch/%u/logs/19aprilcelescope_%A_%a.err

set -euo pipefail

ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)
module load apptainer/1.4.5

SIF="/global/scratch/$USER/celescope_v2.11.4.sif"
REF="/global/scratch/$USER/celescope_GRCh38_ref_output_19april"
FQ_DIR="/global/scratch/$USER/4samplesSRP351272"
BASE_OUTDIR="/global/scratch/$USER/celescope_output_19april"

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
* [celescope run2 reports](https://github.com/GITC2025/integration_benchmark/blob/main/3b2_SRP351272_run2_celescope_reports.pdf)
* unfortunately not much difference with run1, just took longer
* run background remover for this dataset
