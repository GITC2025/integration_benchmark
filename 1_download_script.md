# download Dataset 1
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
## checksum match

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
## dataset 1 md5 checksum 
[shortcuts](https://github.com/GITC2025/integration_benchmark/blob/main/shortcuts.md)

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
## seqkit check

* all num_seqs for same SRR should be exactly the same - paired end reads
* if you do SRA convert to fastq.gz you may sometimes get _2 and _3 fastq.gz as _1 fastq was an index read and discarded

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
## check 8 lines
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

SRR17259462_2.fastq.gz
CAGTAATTATTAAGAACACTTTCCCTTTAAAGGAAACAAAAGTGAATACCATATTGTTTTTACTGTCATAGTGTTGCTTTCTTGCCTGTCCTGCTTAGTTTTTACTTGCTGGATGATACCATAATGTATCAAGGAGCGTCCATGGATACA
FFFFFFFFFFFFF:FFFFFFFFFFF,FFFFFFFFFFFFFF:FFFFFFFFFFFFFF,:FFFFF:FF:FFFFFFFFF:FFFFFFFFF:FFF,F:F,FFFFFFFFFFFFF:F:FFFFFFFF:FFFFFFFFFFFFFFFFFFF,::FF,FFFFFF
CACGTGGATTGAAGCCTAGATCGAAGAGCGTCGTGTAGGGAAAGAGTGGAGTTGTCGGGGGGGGGCGGGGGGTTTGGATTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
FFFF,F:F:FFFFFFFF:FFFFFFFFFFFF,FF:FFF:FF,FFF,FFF,,F,,,,,FF,F:,:FF,:F,,,,,,,,,,,,:FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFF
```
* fastq 1 contains barcodes and UMIs
* each line is followed by quality score line (Phred 33 in this case - auto detected by fastQC tool later)


| File                     | Content                                                         | Length (10X Genomics)|
| ------------------------ | --------------------------------------------------------------- | ---------------------|
| SRR12539462_1.fastq (R1) | Cell barcode (16 bp) + UMI (10-12 bp) + poly-T/adapter (~10 bp) | ~26-91 bp |
| SRR12539462_2.fastq (R2) | Actual cDNA transcript (gene expression data)                   | 91-150 bp |


* _1.fastq.gz has barcode/UMI
* _2.fastq.gz actual genomic sequence and quality score (Phred)
* peek at Phred score to know if it's phred 33 or phred 64
* this one has phred 33
* https://people.duke.edu/~ccc14/duke-hts-2018/bioinformatics/quality_scores.html

# Download dataset 2
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

# didn't create or pt to new directory - add line as needed

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
## dataset 2 md5 checksum 

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
# md5check match shortcut
md5check ena_checksums.md5 downloaded_checksums.md5

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

## seqkit check for pairs num seq match

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
