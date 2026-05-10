## gene transfer format categories Ensembl 115 Grch38 primary assembly
```bash
wget http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
```
```bash
cut -f9 Homo_sapiens.GRCh38.115.gtf | grep -o 'gene_biotype "[^"]*"' | sort | uniq

# if compressed, use zcat
zcat Homo_sapiens.GRCh38.115.gtf.gz | cut -f9 | cut -d' ' -f2,3 | sort | uniq

gene_biotype "IG_J_gene"       
gene_biotype "IG_J_pseudogene"
gene_biotype "IG_pseudogene"
gene_biotype "IG_V_gene"
gene_biotype "IG_V_pseudogene"
gene_biotype "lncRNA"
gene_biotype "miRNA"
gene_biotype "misc_RNA"
gene_biotype "Mt_rRNA"
gene_biotype "Mt_tRNA"
gene_biotype "processed_pseudogene"
gene_biotype "protein_coding" 
gene_biotype "pseudogene"
gene_biotype "ribozyme"
gene_biotype "rRNA"
gene_biotype "rRNA_pseudogene"
gene_biotype "scaRNA"
gene_biotype "snoRNA"
gene_biotype "snRNA"
gene_biotype "sRNA"
gene_biotype "TEC"
gene_biotype "transcribed_processed_pseudogene"
gene_biotype "transcribed_unitary_pseudogene"
gene_biotype "transcribed_unprocessed_pseudogene"
gene_biotype "translated_processed_pseudogene"
gene_biotype "TR_C_gene"
gene_biotype "TR_D_gene"
gene_biotype "TR_J_gene"
gene_biotype "TR_J_pseudogene"
gene_biotype "TR_V_gene"
gene_biotype "TR_V_pseudogene"
gene_biotype "unitary_pseudogene"
gene_biotype "unprocessed_pseudogene"
gene_biotype "vault_RNA"
```
## explanations
* **key categories to keep: protein_coding, lncRNA and immune related genes**
* IG : immunoglobulin
* Mt: mitochondrial
* TR: T-cell receptor
* TEC: to be experimentally confirmed (low confidence seq., usu noise or artifacts)
* prior to ensembl 99, biotype antisense was used for some lncRNAs, now merged into lncRNAs
