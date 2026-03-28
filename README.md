# Integration benchmark scRNA seq analysis

_figuring out the HPC, babysitting the downloads, and benchmarking integration • HSCI591 W2026 BHSC, Queen's Uni_

**Using Positron to have an all-in-one environment for SSH connection to HPC, run R in Apptainer, use Python natively etc.**

* [**all shortcuts**](#shortcuts)
* [**sequencing protocols**](https://github.com/GITC2025/integration_benchmark/wiki/Protocols)

# Dataset 1

Lai, H., Cheng, X., Liu, Q., Luo, W., Liu, M., Zhang, M., Miao, J., Ji, Z., Lin, G. N., Song, W., Zhang, L., Bo, J., Yang, G., Wang, J., & Gao, W. Q. (2021). Single-cell RNA sequencing reveals the epithelial cell heterogeneity and invasive subpopulation in human bladder cancer. International journal of cancer, 149(12), 2099–2115. https://doi.org/10.1002/ijc.33794

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135337 (for SRA files)
https://www.ebi.ac.uk/ena/browser/view/SRP217277 (for fastq.gz)

* 7 primary tumour samples + 1 normal sample = 8 SRA files
* SRA contains: metadata, metadata objects, and raw sequence data, encompassing SRP, SRS, SRX, and SRR
* 10X Genomics seq - we expect paired ends (fastq 1 is barcode/UMI, fast 2 is actual genome)
* for this analysis download fastq.gz as standard practice
* SRA download and conversion to fastq.gz pipeline [here](https://github.com/GITC2025/SRA-convert-to-fastq)
  
```bash
#!/bin/bash
#SBATCH --job-name=batch1_xarg4HTTP
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=batch1xarg4HTTP_sra_out_%j.log
#SBATCH --error=batch1xarg4_sra_err_%j.log
#SBATCH --mail-user=username@email.com
#SBATCH --mail-type=END,FAIL

cd /global/scratch/$USER/12marchENA_GSE

# EOF is master command
# for larger downloads use a transfer node - FTP/aspera etc.

xargs -n 1 -P 4 wget -nc -q --show-progress << 'EOF'
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/002/SRR9897622/SRR9897622_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/063/SRR12539463/SRR12539463_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/003/SRR9897623/SRR9897623_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/062/SRR12539462/SRR12539462_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/001/SRR9897621/SRR9897621_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/058/SRR14615558/SRR14615558_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/005/SRR9897625/SRR9897625_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/062/SRR12539462/SRR12539462_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/004/SRR9897624/SRR9897624_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/063/SRR12539463/SRR12539463_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/058/SRR14615558/SRR14615558_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/005/SRR9897625/SRR9897625_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/001/SRR9897621/SRR9897621_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/002/SRR9897622/SRR9897622_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/004/SRR9897624/SRR9897624_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/003/SRR9897623/SRR9897623_2.fastq.gz
EOF

echo "Generating MD5 checksums..."
find . -maxdepth 1 -type f -name '*.fastq.gz' -print0 \
  | xargs -0 -n1 -P4 md5sum > ENA_GSE_checksums.md5.tmp \
  && mv ENA_GSE_checksums.md5.tmp ENA_GSE_checksums.md5
echo "Checksum generation complete."
```
# checksum match

```bash
# download ENA original checksum
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP217277&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv&limit=0" \
| awk -F'\t' '
NR==1 {
  for (i=1;i<=NF;i++) { if ($i=="fastq_ftp") f=i; if ($i=="fastq_md5") m=i }
  next
}
{
  n = split($(f), urls, ";")
  nm = split($(m), mds, ";")
  for (i=1;i<=n;i++) {
    last = split(urls[i], parts, "/")
    md = (i<=nm ? mds[i] : "")
    if (md != "" && md != "-") print md "  " parts[last]
  }
}
' | LC_ALL=C sort -k2,2 > downloaded_ENA_GSE_checksums.md5

# ENA online checksums - sort by SRR number as needed and overwrite the md5 file
cat downloaded_ENA_GSE_checksums.md5

933f127475b63b99b7efadae016e3341  SRR12539462_1.fastq.gz
5571a9c0d43fea7e81f18effe74f7b41  SRR12539462_2.fastq.gz
764d58d6cc38dfdd3291dddcef6d6ecd  SRR12539463_1.fastq.gz
0d03914dc8ae8f63a2ea2ed1848ff17c  SRR12539463_2.fastq.gz
6b96bac8b3aa8054123727eb83175ce6  SRR14615558_1.fastq.gz
ce7739585bec2a444309f18dbd215de4  SRR14615558_2.fastq.gz
a0a34e2ac001de6d853d5da7c21ac578  SRR9897621_1.fastq.gz
8cc84af901016a68412b0b84f90272a5  SRR9897621_2.fastq.gz
c41c7f08df4f7c72993a45ca78a4c955  SRR9897622_1.fastq.gz
de4055f3a3405e5dcb1dfe601a62ccb3  SRR9897622_2.fastq.gz
41134f5fb9e687ea8c11c310b1be9183  SRR9897623_1.fastq.gz
e947c3f8f950a0800e00f366b64c088c  SRR9897623_2.fastq.gz
d132d15355b72345fdd8f3b27169a87d  SRR9897624_1.fastq.gz
35a57e8bf023b74da3a7eef14f2d680c  SRR9897624_2.fastq.gz
859854cd12660fefed73e4e94ce17354  SRR9897625_1.fastq.gz
9267d8354cb2c141d52f8ce4d5c8cebd  SRR9897625_2.fastq.gz

# locally generated checksums
cat ENA_GSE_checksums.md5

933f127475b63b99b7efadae016e3341  SRR12539462_1.fastq.gz
5571a9c0d43fea7e81f18effe74f7b41  SRR12539462_2.fastq.gz
764d58d6cc38dfdd3291dddcef6d6ecd  SRR12539463_1.fastq.gz
0d03914dc8ae8f63a2ea2ed1848ff17c  SRR12539463_2.fastq.gz
6b96bac8b3aa8054123727eb83175ce6  SRR14615558_1.fastq.gz
ce7739585bec2a444309f18dbd215de4  SRR14615558_2.fastq.gz
a0a34e2ac001de6d853d5da7c21ac578  SRR9897621_1.fastq.gz
8cc84af901016a68412b0b84f90272a5  SRR9897621_2.fastq.gz
c41c7f08df4f7c72993a45ca78a4c955  SRR9897622_1.fastq.gz
de4055f3a3405e5dcb1dfe601a62ccb3  SRR9897622_2.fastq.gz
41134f5fb9e687ea8c11c310b1be9183  SRR9897623_1.fastq.gz
e947c3f8f950a0800e00f366b64c088c  SRR9897623_2.fastq.gz
d132d15355b72345fdd8f3b27169a87d  SRR9897624_1.fastq.gz
35a57e8bf023b74da3a7eef14f2d680c  SRR9897624_2.fastq.gz
859854cd12660fefed73e4e94ce17354  SRR9897625_1.fastq.gz
9267d8354cb2c141d52f8ce4d5c8cebd  SRR9897625_2.fastq.gz
```
# dataset 1 md5 checksum 

**[shell fxn shortcut below](#md5check-shortcut)**
```bash
# match check
awk '
NR==FNR {
  ena[$2] = $1
  next
}
{
  loc = $1
  $1 = ""
  sub(/^ +/, "")
  f = $0
  print f, (ena[f] == loc ? "MATCH" : "MISMATCH")
}
' ENA_GSE_checksums.md5 downloaded_ENA_GSE_checksums.md5

SRR12539462_1.fastq.gz MATCH
SRR12539462_2.fastq.gz MATCH
SRR12539463_1.fastq.gz MATCH
SRR12539463_2.fastq.gz MATCH
SRR14615558_1.fastq.gz MATCH
SRR14615558_2.fastq.gz MATCH
SRR9897621_1.fastq.gz MATCH
SRR9897621_2.fastq.gz MATCH
SRR9897622_1.fastq.gz MATCH
SRR9897622_2.fastq.gz MATCH
SRR9897623_1.fastq.gz MATCH
SRR9897623_2.fastq.gz MATCH
SRR9897624_1.fastq.gz MATCH
SRR9897624_2.fastq.gz MATCH
SRR9897625_1.fastq.gz MATCH
SRR9897625_2.fastq.gz MATCH
```
# seqkit check

* all num_seqs for same SRR should be exactly the same - paired end reads
* later files being _2 and _3 : it happened that _1.fastq was an index read (library tag) and was discarded
* as long as we have 2 paired reads per SRR with identical num_seqs
* check num_seqs against vdp dump metadata SEQ with code for longer lists

```yaml
# dataset 1 seqkit check
seqkit stat -j 8 ./*.fastq.gz

# dataset 1 GSE135337 seqkit check output
[mii] loading StdEnv/2023 seqkit/2.5.1 ...
processed files:  16 / 16 [======================================] ETA: 0s. done
file                    format  type     num_seqs         sum_len  min_len  avg_len  max_len
SRR9897621_1.fastq.gz   FASTQ   DNA   351,909,721   9,149,652,746       26       26       26
SRR9897621_2.fastq.gz   FASTQ   DNA   351,909,721  53,138,367,871      151      151      151
SRR9897622_1.fastq.gz   FASTQ   DNA   373,103,899   9,700,701,374       26       26       26
SRR9897622_2.fastq.gz   FASTQ   DNA   373,103,899  56,338,688,749      151      151      151
SRR9897623_1.fastq.gz   FASTQ   DNA   168,058,967   4,369,533,142       26       26       26
SRR9897623_2.fastq.gz   FASTQ   DNA   168,058,967  25,376,904,017      151      151      151
SRR9897624_1.fastq.gz   FASTQ   DNA   374,515,865   9,737,412,490       26       26       26
SRR9897624_2.fastq.gz   FASTQ   DNA   374,515,865  56,551,895,615      151      151      151
SRR9897625_1.fastq.gz   FASTQ   DNA   307,211,122   7,987,489,172       26       26       26
SRR9897625_2.fastq.gz   FASTQ   DNA   307,211,122  46,388,879,422      151      151      151
SRR12539462_1.fastq.gz  FASTQ   DNA   359,075,221   9,335,955,746       26       26       26
SRR12539462_2.fastq.gz  FASTQ   DNA   359,075,221  54,220,358,371      151      151      151
SRR12539463_1.fastq.gz  FASTQ   DNA   367,846,395   9,564,006,270       26       26       26
SRR12539463_2.fastq.gz  FASTQ   DNA   367,846,395  55,544,805,645      151      151      151
SRR14615558_1.fastq.gz  FASTQ   DNA   374,150,612   9,727,915,912       26       26       26
SRR14615558_2.fastq.gz  FASTQ   DNA   374,150,612  56,496,742,412      151      151      151

```
# check 8 lines
```bash
for f in *.fastq.gz; do
  echo "$f"
  zcat "$f" | awk 'NR<=8 && (NR%4==2 || NR%4==0)'
  echo
done
```
```bash
# output
SRR17259462_1.fastq.gz
NTCAATGAATCCACGTGCTTGAGAACTATGCATCAGCATGCGGCTACGAACGCGTACTCATCAAGCTATTTTTTTTTTTTTTTTTTTTGTATTATTTCCGAAGCAAAATATACTCATCAAATTCAAATAAAAATGCCATAGTGTATTTGA
#FFFFFF:FFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFF:,:FFFF:FFFFFFF:FFFFFFFFFFFFFFF:FFFF,FFFFFFFF:FFFF,FFFFFFF:F,:,F,:F,FFFFFFFFFFFFFF,:F,:F,:FFFF:F,,
NGGCTTCAATCCACGTGCTGTCTCTTATACACATCTCCGAGCCCACGAGACTCTGCTGTATCTCGTATGCCGTCTTCTGCTTGAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFF:FFFFFF:FFFF:,,,::F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
* fastq 1 contains barcodes and UMIs
* each line is followed by quality score line (Phred 33 in this case - auto detected by fastQC tool later)

```

| File                     | Content                                                         | Length (10X Genomics)|
| ------------------------ | --------------------------------------------------------------- | ---------------------|
| SRR12539462_1.fastq (R1) | Cell barcode (16 bp) + UMI (10-12 bp) + poly-T/adapter (~10 bp) | ~26-91 bp |
| SRR12539462_2.fastq (R2) | Actual cDNA transcript (gene expression data)                   | 91-150 bp |


* _1.fastq.gz has barcode/UMI
* _2.fastq.gz actual genomic sequence and quality score (Phred)
* peek at Phred score to know if it's phred 33 or phred 64
* this one has phred 33
* https://people.duke.edu/~ccc14/duke-hts-2018/bioinformatics/quality_scores.html

***
# Dataset 2

**4 samples total

2x pri BCa, 1 recurrent BCa and one cystitis glandularis**

Luo, Y., Tao, T., Tao, R., Huang, G., & Wu, S. (2022). Single-Cell Transcriptome Comparison of Bladder Cancer Reveals Its Ecosystem. Frontiers in oncology, 12, 818147. https://doi.org/10.3389/fonc.2022.818147

https://www.ebi.ac.uk/ena/browser/view/SRP351272

We download the fastq.gz directly from ENA this time, xarg for 4 parallel downloads.

```bash
#!/bin/bash
#SBATCH --job-name=4samplesHTTPxarg4
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=4samplesHTTP_sra_out_%j.log
#SBATCH --error=4samplesHTTP_sra_err_%j.log
#SBATCH --mail-user=username@email.com
#SBATCH --mail-type=END,FAIL

# didn't create or pt to new directory


xargs -n 1 -P 4 wget -nc -q --show-progress << 'EOF'
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/062/SRR17259462/SRR17259462_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/063/SRR17259463/SRR17259463_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/065/SRR17259465/SRR17259465_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/064/SRR17259464/SRR17259464_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/062/SRR17259462/SRR17259462_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/063/SRR17259463/SRR17259463_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/064/SRR17259464/SRR17259464_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/065/SRR17259465/SRR17259465_2.fastq.gz
EOF

# compute local checksum for checking later
echo "Generating MD5 checksums..."
md5sum *.fastq.gz > downloaded_checksums.md5
echo "Checksum generation complete."

```
# dataset 2 md5 checksum 

([shell fxn shortcut below](#md5check-shortcut))

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP351272&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv&limit=0" \
| awk -F'\t' '
NR==1 {
  for (i=1;i<=NF;i++) { if ($i=="fastq_ftp") f=i; if ($i=="fastq_md5") m=i }
  next
}
{
  n = split($(f), urls, ";")
  nm = split($(m), mds, ";")
  for (i=1;i<=n;i++) {
    last = split(urls[i], parts, "/")
    md = (i<=nm ? mds[i] : "")
    if (md != "" && md != "-") print md "  " parts[last]
  }
}
' | LC_ALL=C sort -k2,2 > ena_checksums.md5
```
```
# md5checksums
cat downloaded_checksums.md5

3fe6b6c7a818a249df2d2d3a06f8410c  SRR17259462_1.fastq.gz
f32e69227a30d2f89adcae45d0bd2730  SRR17259462_2.fastq.gz
3320d2287ef1086260cbbc1b5dcf034f  SRR17259463_1.fastq.gz
7f7b117baf3f96a58929fec2054224b9  SRR17259463_2.fastq.gz
9c02719c112f38c003cb4ae1477ba4ed  SRR17259464_1.fastq.gz
68e6cb9488f450a4c2b46021af0b1837  SRR17259464_2.fastq.gz
74532db5fc0b65effb4321719bfddab4  SRR17259465_1.fastq.gz
fd91c016fd1cd15814044f375b00d3af  SRR17259465_2.fastq.gz

cat ena_checksums.md5

3fe6b6c7a818a249df2d2d3a06f8410c  SRR17259462_1.fastq.gz
f32e69227a30d2f89adcae45d0bd2730  SRR17259462_2.fastq.gz
3320d2287ef1086260cbbc1b5dcf034f  SRR17259463_1.fastq.gz
7f7b117baf3f96a58929fec2054224b9  SRR17259463_2.fastq.gz
68e6cb9488f450a4c2b46021af0b1837  SRR17259464_2.fastq.gz
9c02719c112f38c003cb4ae1477ba4ed  SRR17259464_1.fastq.gz
74532db5fc0b65effb4321719bfddab4  SRR17259465_1.fastq.gz
fd91c016fd1cd15814044f375b00d3af  SRR17259465_2.fastq.gz
```
```bash
# check match
awk '
NR==FNR {
  ena[$2] = $1
  next
}
{
  loc = $1
  $1 = ""
  sub(/^ +/, "")
  f = $0
  print f, (ena[f] == loc ? "MATCH" : "MISMATCH")
}
' ena_checksums.md5 downloaded_checksums.md5

# output
SRR17259462_1.fastq.gz MATCH
SRR17259462_2.fastq.gz MATCH
SRR17259463_1.fastq.gz MATCH
SRR17259463_2.fastq.gz MATCH
SRR17259464_1.fastq.gz MATCH
SRR17259464_2.fastq.gz MATCH
SRR17259465_1.fastq.gz MATCH
SRR17259465_2.fastq.gz MATCH
```

# seqkit check for pairs num seq match

```
module load StdEnv/2023 seqkit/2.5.1

for file in *.fastq.gz; do
  (
    echo "=== Stats for $file ==="
    seqkit stats -j 2 "$file"
    echo ""
  ) &
done
wait

# output
file                    format  type     num_seqs         sum_len  min_len  avg_len  max_len
SRR17259462_1.fastq.gz  FASTQ   DNA   403,701,477   60,555,221,550      150      150      150
SRR17259462_2.fastq.gz  FASTQ   DNA   403,701,477   60,555,221,550      150      150      150
SRR17259463_1.fastq.gz  FASTQ   DNA   379,705,781   56,955,867,150      150      150      150
SRR17259463_2.fastq.gz  FASTQ   DNA   379,705,781   56,955,867,150      150      150      150
SRR17259464_1.fastq.gz  FASTQ   DNA   506,432,186   75,964,827,900      150      150      150
SRR17259464_2.fastq.gz  FASTQ   DNA   506,432,186   75,964,827,900      150      150      150
SRR17259465_1.fastq.gz  FASTQ   DNA   464,865,330   69,729,799,500      150      150      150
SRR17259465_2.fastq.gz  FASTQ   DNA   464,865,330   69,729,799,500      150      150      150
```

check against metadata if you wish

# Run fastqc on fastq.gz
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
* unable to install falco on HPC

# dataset 1 fastqc
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
# dataset 2 fastqc
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
# download FastQC html reports - from local console
```bash
scp -r YOUR_USERNAME@YOUR_SERVER_ADDRESS:/PATH/TO/REMOTE/FOLDER ./LOCAL_DESTINATION_FOLDER
```
<img width="900" height="400" alt="image" src="https://github.com/user-attachments/assets/4a9a927b-5bbe-4c7f-bccf-239590ac993b" />

<img width="900" height="400" alt="image" src="https://github.com/user-attachments/assets/5df125c0-8775-46ef-84f7-c79a1c5c739c" />

# Fastqc_data.txt parsing Dataset 1

the txt are in the zip folders

```bash
# unzip all 
for z in *.zip; do [ -f "$z" ] || continue; unzip -o "$z"; done

# move to current folder
# move each fastqc_data.txt from subdirs into cwd, prefixing with parent-dir to avoid overwrites
find . -mindepth 2 -type f -name 'fastqc_data.txt' -print0 \
  | while IFS= read -r -d '' f; do
      d=$(basename "$(dirname "$f")")
      mv -n "$f" "./${d}_fastqc_data.txt"
    done

# parse for FAIL
awk 'BEGIN{OFS="\t"; print "FILE","MODULE","STATUS"}
/^FAIL[[:space:]]/{
  status=$1; $1=""; sub(/^[[:space:]]+/,""); print FILENAME, $0, status
}
/^>>/{
  line=$0; sub(/^>>/,"",line); split(line,a,"\t");
  if(toupper(a[2])=="FAIL") print FILENAME, a[1], a[2]
}' *.txt

# parse for WARN
awk 'BEGIN{OFS="\t"; print "FILE","MODULE","STATUS"}
/^WARN[[:space:]]/{
  status=$1; $1=""; sub(/^[[:space:]]+/,""); print FILENAME, $0, status
}
/^>>/{
  line=$0; sub(/^>>/,"",line); split(line,a,"\t");
  if(toupper(a[2])=="WARN") print FILENAME, a[1], a[2]
}' *.txt
```
[**fastqc txt parser shortcut below**](#FastQC-txt-parser-shortcut)

```yaml
FILE    MODULE  STATUS
SRR12539462_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR12539462_2_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR12539462_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR12539463_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR12539463_2_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR12539463_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR14615558_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR14615558_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR9897621_1_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897621_2_fastqc_fastqc_data.txt     Per base sequence content       fail
SRR9897621_2_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897621_2_fastqc_fastqc_data.txt     Adapter Content fail
SRR9897622_1_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897622_2_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897623_2_fastqc_fastqc_data.txt     Per base sequence content       fail
SRR9897623_2_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897623_2_fastqc_fastqc_data.txt     Adapter Content fail
SRR9897624_1_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897624_2_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897625_1_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail
SRR9897625_2_fastqc_fastqc_data.txt     Per base sequence content       fail
SRR9897625_2_fastqc_fastqc_data.txt     Sequence Duplication Levels     fail

FILE    MODULE  STATUS
SRR12539462_1_fastqc_fastqc_data.txt    Per base sequence content       warn
SRR12539462_2_fastqc_fastqc_data.txt    Adapter Content warn
SRR12539463_1_fastqc_fastqc_data.txt    Per base sequence content       warn
SRR12539463_2_fastqc_fastqc_data.txt    Overrepresented sequences       warn
SRR12539463_2_fastqc_fastqc_data.txt    Adapter Content warn
SRR14615558_1_fastqc_fastqc_data.txt    Per base sequence content       warn
SRR14615558_2_fastqc_fastqc_data.txt    Per base sequence content       warn
SRR14615558_2_fastqc_fastqc_data.txt    Adapter Content warn
SRR9897621_2_fastqc_fastqc_data.txt     Overrepresented sequences       warn
SRR9897622_1_fastqc_fastqc_data.txt     Per base sequence content       warn
SRR9897622_2_fastqc_fastqc_data.txt     Per base sequence content       warn
SRR9897622_2_fastqc_fastqc_data.txt     Overrepresented sequences       warn
SRR9897622_2_fastqc_fastqc_data.txt     Adapter Content warn
SRR9897623_1_fastqc_fastqc_data.txt     Per base sequence content       warn
SRR9897623_1_fastqc_fastqc_data.txt     Sequence Duplication Levels     warn
SRR9897623_2_fastqc_fastqc_data.txt     Per sequence GC content warn
SRR9897623_2_fastqc_fastqc_data.txt     Overrepresented sequences       warn
SRR9897624_1_fastqc_fastqc_data.txt     Per base sequence content       warn
SRR9897624_2_fastqc_fastqc_data.txt     Per base sequence content       warn
SRR9897624_2_fastqc_fastqc_data.txt     Overrepresented sequences       warn
SRR9897624_2_fastqc_fastqc_data.txt     Adapter Content warn
SRR9897625_1_fastqc_fastqc_data.txt     Per base sequence content       warn
SRR9897625_1_fastqc_fastqc_data.txt     Overrepresented sequences       warn
SRR9897625_2_fastqc_fastqc_data.txt     Per sequence GC content warn
SRR9897625_2_fastqc_fastqc_data.txt     Overrepresented sequences       warn
SRR9897625_2_fastqc_fastqc_data.txt     Adapter Content warn
```

Nothing of concern here since no warn/fail for Per base sequence quality - most important module. 

# dataset 2 fastQC

**[unzipfastq shortcut](#unzipfastq-shortcut)**
```bash
# unzipfastq and move txt to current folder shortcut
unzipfastq

# checks
fastqcfail
fastqcwarn
```
```yaml
# output
FILE    MODULE  STATUS
SRR17259462_1_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259462_1_fastqc_fastqc_data.txt    Per sequence GC content fail
SRR17259462_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259462_1_fastqc_fastqc_data.txt    Adapter Content fail
SRR17259462_2_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259462_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259462_2_fastqc_fastqc_data.txt    Adapter Content fail
SRR17259463_1_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259463_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259463_1_fastqc_fastqc_data.txt    Adapter Content fail
SRR17259463_2_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259463_2_fastqc_fastqc_data.txt    Per sequence GC content fail
SRR17259463_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259463_2_fastqc_fastqc_data.txt    Overrepresented sequences       fail
SRR17259463_2_fastqc_fastqc_data.txt    Adapter Content fail
SRR17259464_1_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259464_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259464_1_fastqc_fastqc_data.txt    Adapter Content fail
SRR17259464_2_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259464_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259464_2_fastqc_fastqc_data.txt    Overrepresented sequences       fail
SRR17259464_2_fastqc_fastqc_data.txt    Adapter Content fail
SRR17259465_1_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259465_1_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259465_2_fastqc_fastqc_data.txt    Per base sequence content       fail
SRR17259465_2_fastqc_fastqc_data.txt    Sequence Duplication Levels     fail
SRR17259465_2_fastqc_fastqc_data.txt    Adapter Content fail
FILE    MODULE  STATUS
SRR17259462_2_fastqc_fastqc_data.txt    Overrepresented sequences       warn
SRR17259463_1_fastqc_fastqc_data.txt    Per sequence GC content warn
SRR17259464_1_fastqc_fastqc_data.txt    Per sequence GC content warn
SRR17259464_2_fastqc_fastqc_data.txt    Per sequence GC content warn
SRR17259465_1_fastqc_fastqc_data.txt    Per sequence GC content warn
SRR17259465_1_fastqc_fastqc_data.txt    Adapter Content warn
SRR17259465_2_fastqc_fastqc_data.txt    Per sequence GC content warn
SRR17259465_2_fastqc_fastqc_data.txt    Overrepresented sequences       warn
```

# Reading FastQC output
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

# Chemistry in brief
* v2 dominant prior to 2018 (16 + 10)
* post 2018 - v3 and 3.1 (16 + 12)
* 2023 onwards - v4 (16 + 12)
* v4 supports feature barcodes - run scRNA parallel w/ other assays
* multiplex for multi-omic output
* improved cell recovery 

Sequencing protocol https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html

v2 doc https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000108_AssayConfiguration_SC3v2.pdf

v3 doc https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf

v4 doc https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000731_ChromiumGEM-X_SingleCell3_ReagentKits_v4_UserGuide_RevA.pdf

https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/10x_LIT000220_Product_Sheet_GEM-X-Single-Cell-Gene-Expression_Letter_Digital.pdf

***
# zcat Dataset 1 

* [zcat shortcuts](#zcat-check-chemistry)

```bash
# zcat extract all R1 lengths
for f in *1.fastq.gz; do
echo "$f:"
zcat "$f" | head -n 10 | awk 'NR%4==2 {print length($0)}' | sort -u
done

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
# zcat Dataset 2

```bash
# use shortcut
zcatR1

# output
SRR17259462_1.fastq.gz:
150
SRR17259463_1.fastq.gz:
150
SRR17259464_1.fastq.gz:
150
SRR17259465_1.fastq.gz:
150

# this is not 10x genomics! 

# use shortcut
zcatR2

# output 
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
* Use this: https://github.com/singleron-RD/scrna
* good stress test for cross-platform integration benchmarking later

# Alignment
* Cell Ranger 10.0.0 
* map to 10x GRCh38 2024-A reference genome

  
# Next steps
1. raw QC - FastQC - complete
2. zcat check chemistry - complete 
   Align and quantify - CellRanger
3. QC with sctk
4. Normalization and feature selection - SCTrasnform, ID HVGs (2000 highly variable genes)
5. 6. Findmarkers Seurat for biological grounding
6. Dimension reduction - PCA, scree plot to check to find where PC levels off, stop when PC# covers 70-90% of variance
7. Integration benchmarking here - control vs scArch vs Harmony etc. benchmark across diff categories of tools
        *evaluate against ASQ, ARI, LISI etc.
8. UMAP as qualitative benchmark - set seed, try out diff seeds, try out diff k's to find stability and convergence
	Findmarkers Seurat for biological grounding


# PCA in a nutshell

<img width="800" height="440" alt="PCAsteps" src="https://github.com/user-attachments/assets/adcc1405-3f5c-4d00-b761-0593a56f8734" />

***
# shortcuts

```bash
alias myjobs='queue -u "$USER'
alias smallnode='salloc --time=04:00:00 --cpus-per-task=8 --mem=16G'
alias bignode='salloc --time=04:00:00 --cpus-per-task=8 --mem=128G'
# dataset 1
alias GSE='cd /global/scratch/$USER/12marchENA_GSE'

# dataset 2
alias SRP='cd /global/scratch/$USER/4samplesSRP351272'

# bind the alias to persist across sessions
echo "alias myalias='command here'" >> ~/.bashrc
source ~/.bashrc
```

# bind shell fxns
```bash
cat >> ~/.bashrc <<'EOF'

myfunc() {
  # function body
  echo "doing something"
}
EOF

source ~/.bashrc
```

# md5check shortcut
```bash
md5check() {
  local local_file="${1:-localchecksumfile.md5}"
  local target_file="${2:-targetchecksumfile.md5}"

  awk '
  NR==FNR {
    ena[$2] = $1
    next
  }
  {
    loc = $1
    $1 = ""
    sub(/^ +/, "")
    f = $0
    print f, (ena[f] == loc ? "MATCH" : "MISMATCH")
  }
  ' "$target_file" "$local_file"
}

# use shortcut
md5check localchecksumfile.md5 targetchecksumfile.md5
```

# unzipfastq shortcut
```bash
# unzip fastq.zip and move fastqc txt into current folder
cat >> ~/.bashrc <<'EOF'

unzipfastq() {
  # unzip all zip files in current directory
  for z in *.zip; do
    [ -f "$z" ] || continue
    unzip -o "$z"
  done

  # move each fastqc_data.txt from subdirectories into current directory
  # prefix with parent directory name to avoid collisions
  find . -mindepth 2 -type f -name 'fastqc_data.txt' -print0 \
    | while IFS= read -r -d '' f; do
        d=$(basename "$(dirname "$f")")
        mv -n "$f" "./${d}_fastqc_data.txt"
      done
}
EOF

source ~/.bashrc
```


# FastQC txt parser shortcut
```bash
fastqcfail() {
  shopt -s nullglob
  local files=( *.txt )
  [ ${#files[@]} -gt 0 ] || { echo "no *.txt files here"; return 1; }

  awk 'BEGIN{OFS="\t"; print "FILE","MODULE","STATUS"}
  /^FAIL[[:space:]]/{
    status=$1; $1=""; sub(/^[[:space:]]+/,""); print FILENAME, $0, status
  }
  /^>>/{
    line=$0; sub(/^>>/,"",line); split(line,a,"\t");
    if(toupper(a[2])=="FAIL") print FILENAME, a[1], a[2]
  }' "${files[@]}"
}
```
```bash
fastqcwarn() {
  shopt -s nullglob
  local files=( *.txt )
  [ ${#files[@]} -gt 0 ] || { echo "no *.txt files here"; return 1; }

  awk 'BEGIN{OFS="\t"; print "FILE","MODULE","STATUS"}
  /^WARN[[:space:]]/{
    status=$1; $1=""; sub(/^[[:space:]]+/,""); print FILENAME, $0, status
  }
  /^>>/{
    line=$0; sub(/^>>/,"",line); split(line,a,"\t");
    if(toupper(a[2])=="WARN") print FILENAME, a[1], a[2]
  }' "${files[@]}"
}
```
```bash
# use shortcuts
fastqcfail
fastqcwarn
```
# zcat check chemistry
```bash
# check zcatR1
cat >> ~/.bashrc <<'EOF'
zcatR1() {
  local files=("$@")
  if [ ${#files[@]} -eq 0 ]; then files=( *1.fastq.gz ); fi
  for f in "${files[@]}"; do
    echo "$f:"
    zcat "$f" | head -n 10 | awk 'NR%4==2 {print length($0)}' | sort -u
  done
}
EOF
source ~/.bashrc

# check zcatR2
cat >> ~/.bashrc <<'EOF'
zcatR2() {
  local files=("$@")
  if [ ${#files[@]} -eq 0 ]; then files=( *2.fastq.gz ); fi
  for f in "${files[@]}"; do
    echo "$f:"
    zcat "$f" | head -n 10 | awk 'NR%4==2 {print length($0)}' | sort -u
  done
}
EOF
source ~/.bashrc

# shortcuts
zcatR1
zcatR2
```
