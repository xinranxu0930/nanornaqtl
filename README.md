# nanornaqtl

[![PyPI version](https://badge.fury.io/py/nanornaqtl.svg)](https://badge.fury.io/py/nanornaqtl)
[![Python versions](https://img.shields.io/pypi/pyversions/nanornaqtl.svg)](https://pypi.org/project/nanornaqtl/)

A comprehensive toolkit for molecular phenotyping and QTL analysis using Nanopore direct RNA sequencing data


## Introduction

**nanornaqtl** is an analysis tool specifically designed for Nanopore direct RNA sequencing data, used to identify multiple RNA molecular phenotypes and perform population-level QTL (Quantitative Trait Loci) analysis. This tool is particularly suitable for analyzing pooled samples (sequencing of multiple mixed individuals), with all analyses performed at the read level.

### Main Features

- **8 Types of Molecular Phenotype Identification** 📊:
  - RNA modification sites: m6A, m5C, pseudouridine (pseU), inosine
  - polyA tail length
  - Intron retention rate
  - Alternative polyadenylation site (APA)
  - Transcript isoform

- **8 Types of QTL Analysis** 🧪:
  - Modification QTL (m6A, m5C, pseU, inosine)
  - APA usage pattern QTL
  - Isoform usage pattern QTL
  - PolyA tail length QTL
  - Intron retention rate QTL

### Features ✨

- **Read-level analysis**: Suitable for population-level analysis of pooled samples
- **Bayesian statistical methods**: More robust for uneven coverage in Nanopore data
- **Parallel processing**: Supports multi-threading acceleration with chromosome-based parallel processing
- **Flexible parameter settings**: Customizable quality thresholds, coverage requirements, etc.

---

## Installation

### Install via pip 📦

```bash
pip install nanornaqtl
```

### Install from source

```bash
git clone https://github.com/xinranxu0930/nanornaqtl.git
cd nanornaqtl
poetry install
```

### Dependencies

#### External tools:

- **bedtools** (v2.30.0+)
- **samtools** (v1.15+)
- **IsoQuant** (for isoform analysis): [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant)

#### Python packages:

- pysam
- pymc
- arviz
- statsmodels
- pandas
- numpy
- scipy
- logomaker
- biopython

---

## Data Preprocessing

Before using nanornaqtl, Nanopore sequencing data needs to be preprocessed.

### Step 1: Data Format Conversion

If the sequencing data is in FAST5 format, it needs to be converted to POD5 format first (skip this step if already in POD5):

```bash
# Convert using pod5 tool
pod5 convert fast5 input.fast5 --output output.pod5
```

### Step 2: Basecalling and Alignment 🔧

Use [Dorado](https://github.com/nanoporetech/dorado) for basecalling and alignment:

#### Basic Alignment (Required)

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --reference fasta.fa \
  --mm2-opts "-x splice -k 14" \
  > output.bam
```

**Parameter explanation**:
- `--mm2-opts "-x splice -k 14"`: RNA splice alignment parameters

#### Optional Features

**If modification site analysis is needed**:

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --reference fasta.fa \
  --mm2-opts "-x splice -k 14" \
  --modified-bases-models model_name \
  > output.bam
```

Download modification recognition models:
```bash
# Download all models
dorado download --model all

# Or download specific model
dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.1.0
```

**If polyA tail length analysis is needed**:

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --reference fasta.fa \
  --mm2-opts "-x splice -k 14" \
  --estimate-poly-a \
  > output.bam
```

---

## Quick Start

Here is a complete analysis workflow example:

```bash
# 1. Data preprocessing
nanornaqtl prep \
  -b raw_output.bam \
  -p sample01 \
  -o ./work_dir \
  -t 8 \
  -q 10

# 2. Identify m6A modification sites
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -p sample01_m6A \
  -o ./work_dir \
  -t 20 \
  -f hg19.fa \
  --motif \
  --metaPlotR

# 3. Perform m6A QTL analysis
nanornaqtl qtl m6A \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_m6A_qtl \
  -o ./work_dir \
  --modification_csv sample01_m6A_sites_result.csv \
  --read_mod_dict sample01_m6A_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

## Detailed Usage Instructions

### 1. prep - Data Preprocessing

**Function**: Filter unmapped reads from raw BAM file, generate FASTQ and strand-separated BAM files for mapped reads.

**Input**:
- Raw BAM file output from Dorado (no preprocessing required)

**Command**:

```bash
nanornaqtl prep \
  -b <basecall_bam> \
  -p <output_prefix> \
  -o <output_dir> \
  -t <threads> \
  -q <min_mapq>
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `-b, --bam` | Input raw BAM file path | Required |
| `-p, --prefix` | Output file prefix name | Required |
| `-o, --output_dir` | Directory for all output files, will be created if not exists | Required |
| `-t, --threads` | Number of threads | 4 |
| `-q, --min_mapq` | Minimum read mapping quality score | 0 |

**Output files**:

- `<prefix>_calls_sorted_map.bam`: Merged BAM of mapped reads
- `<prefix>_calls_sorted_map0.bam`: Plus strand BAM (flag=0)
- `<prefix>_calls_sorted_map16.bam`: Minus strand BAM (flag=16)
- `<prefix>_calls_sorted_map.fastq`: FASTQ of mapped reads
- Corresponding `.bai` index files

**Example**:

```bash
nanornaqtl prep \
  -b raw_output.bam \
  -p sample01 \
  -o ./work_dir \
  -t 8 \
  -q 10
```

---

### 2. pheno - Molecular Phenotype Identification

#### 2.1 m6A Modification Site Identification

**Function**: Identify m6A (N6-methyladenosine) modification sites.

**Command**:

```bash
nanornaqtl pheno m6A \
  -b <map_bam> \
  -o <output_dir> \
  -p <output_prefix> \
  -t <threads> \
  -f <fasta.fa> \
  [--motif] \
  [--metaPlotR]
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `-b, --bam` | map.bam file path from prep step (`*_map.bam`) | Required |
| `-o, --output_dir` | Output file directory, will be created if not exists | Required |
| `-p, --prefix` | Output file prefix name | Required |
| `-t, --threads` | Number of threads | 4 |
| `--mod_threshold` | Modification probability threshold | 0.75 |
| `--min_qscore` | Minimum base quality score for modification site | 10 |
| `--min_mapq` | Minimum read mapping quality (MAPQ) | 0 |
| `-f, --fasta` | Reference genome fasta file path | Required |
| `--min_rate` | Minimum modification rate | 0.1 |
| `--min_cov` | Minimum coverage for modification site | 5 |
| `--motif` | Enable motif filtering (DRACH: `[GAT][GA]AC[ATC]`) | False |
| `--motifPaint` | Enable motif plot generation | False |
| `--metaPlotR` | Generate bed file for [metaPlotR](https://github.com/olarerin/metaPlotR) | False |

**Output files**:

- `<prefix>_m6A_sites_result.csv`: m6A modification site information

**Output column description**:

| Column name | Description |
|------|------|
| `chrom` | Chromosome |
| `pos_1base` | Position (1-based) |
| `strand` | Strand direction (+/-) |
| `mod_num` | Number of modified reads |
| `cov` | Total coverage |
| `mod_rate` | Modification rate |
| `motif` | Motif sequence (if --motif enabled) |

**Example**:

```bash
# With motif filtering
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -p sample01_m6A \
  -o ./work_dir \
  -t 20 \
  -f hg19.fa \
  --motif \
  --metaPlotR

# Without motif filtering
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -p sample01_m6A \
  -o ./work_dir \
  -f hg19.fa \
  -t 20
```

---

#### 2.2 m5C Modification Site Identification

**Function**: Identify m5C (5-methylcytosine) modification sites.

**Command**:

```bash
nanornaqtl pheno m5C \
  -b <map_bam> \
  -p <output_prefix> \
  -o <output_dir> \
  -t <threads> \
  -f <fasta.fa> \
  [--motif] \
  [--metaPlotR]
```

**Parameters**: Same as m6A

**Motif classification** (if --motif enabled):
- **CG**: CpG sites
- **CHG**: CHG sites
- **CHH**: CHH sites

**Output files**:

- `<prefix>_m5C_sites_result.csv`: m5C modification site information

**Output column description**: Same as m6A, with additional `motif_classification` column (CG/CHG/CHH)

**Example**:

```bash
nanornaqtl pheno m5C \
  -b sample01_calls_sorted_map.bam \
  -p sample01_m5C \
  -o ./work_dir \
  -t 20 \
  -f <fasta.fa> \
  --motif
```

---

#### 2.3 pseudouridine (pseU) Modification Site Identification

**Function**: Identify pseudouridine modification sites.

**Command**:

```bash
nanornaqtl pheno pseU \
  -b <map_bam> \
  -p <output_prefix> \
  -o <output_dir> \
  -t <threads> \
  -f <fasta.fa> \
  [--motif] \
  [--metaPlotR]
```

**Parameters**: Same as m6A

**Motif classification** (if --motif enabled):
- **pus1**: `[ACT][AG]T`
- **pus4**: `GTTC[ATCG]A`
- **pus7**: `TGTA[AG]`

**Output files**:

- `<prefix>_pseU_sites_result.csv`: pseU modification site information

**Output column description**: Same as m5C, includes motif classification (pus1/pus4/pus7)

**Example**:

```bash
nanornaqtl pheno pseU \
  -b sample01_calls_sorted_map.bam \
  -p sample01_pseU \
  -o ./work_dir \
  -f hg19.fa \
  -t 20 \
  --motif
```

---

#### 2.4 inosine Modification Site Identification

**Function**: Identify inosine (A-to-I editing) sites.

**Command**:

```bash
nanornaqtl pheno inosine \
  -b <map_bam> \
  -p <output_prefix> \
  -o <output_dir> \
  -f <fasta.fa> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**Parameters**: Same as m6A

**Motif** (if --motif enabled): `TA[GT]`

**Output files**:

- `<prefix>_inosine_sites_result.csv`: inosine site information

**Example**:

```bash
nanornaqtl pheno inosine \
  -b sample01_calls_sorted_map.bam \
  -p sample01_inosine \
  -o ./work_dir \
  -f hg19.fa \
  -t 20 \
  --motif
```

---

#### 2.5 polyA Tail Length Identification

**Function**: Extract polyA tail length for each read.

**Prerequisites**: BAM file must contain pt tag (generated using `--estimate-poly-a` during Dorado basecall)

**Command**:

```bash
nanornaqtl pheno polyA_tail \
  -b <map_bam> \
  -p <output_prefix> \
  -o <output_dir> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `-b, --bam` | map.bam file path from prep step | Required |
| `-p, --prefix` | Output file prefix name | Required |
| `-o, --output_dir` | Output directory path, will be created if not exists | Required |
| `-t, --threads` | Number of threads | 4 |

**Output files**:

- `<prefix>_polyAlen_result.csv`: polyA tail length information

**Output column description**:

| Column name | Description |
|------|------|
| `readID` | Read identifier |
| `polyA_length` | polyA tail length |

**Example**:

```bash
nanornaqtl pheno polyA_tail \
  -b sample01_calls_sorted_map.bam \
  -p sample01_polyA \
  -o ./work_dir \
  -t 20
```

---

#### 2.6 Intron Retention Rate Identification

**Function**: Calculate intron retention rate for each read.

**Command**:

```bash
nanornaqtl pheno intron_retention \
  -g <gtf_file> \
  -b <map_bam> \
  -o <output_dir> \
  -p <output_prefix>
```

**Parameter explanation**:

| Parameter | Description |
|------|------|
| `-g, --gtf` | GTF annotation file path |
| `-b, --bam` | map.bam file path from prep step |
| `-p, --prefix` | Output file prefix name |
| `-o, --output_dir` | Output directory path, will be created if not exists |

**Output files**:

- `<prefix>_intronRetention_result.csv`: Intron retention rate information

**Output column description**:

| Column name | Description |
|------|------|
| `readID` | Read identifier |
| `exon_len` | Exon length |
| `intron_len` | Intron length |
| `len` | Total length |
| `IntronRetentionRate` | Intron retention rate |

**Example**:

```bash
nanornaqtl pheno intron_retention \
  -g gencode.v47lift37.annotation.gtf \
  -b sample01_calls_sorted_map.bam \
  -o ./results \
  -p sample01
```

---

#### 2.7 Alternative polyadenylation (APA) Site Identification

**Function**: Identify alternative polyadenylation sites used by each read. Supports two modes: annotation based on APAdb database, or de novo APA site identification.

**Command**:

Mode 1: Using APAdb database
```bash
nanornaqtl pheno APA \
  -b  <map_bam> \
  -p  <output_prefix> \
  -o  <output_dir> \
  --apadb <apadb_file>
```

Mode 2: De novo APA site identification
```bash
nanornaqtl pheno APA \
  -b  <map_bam> \
  -p  <output_prefix> \
  -o  <output_dir> \
  -f  <fasta_file> \
  -d  <distance> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description |
|------|------|
| `-b, --bam` | map.bam file path from prep step |
| `-p, --prefix` | Output file prefix |
| `-o, --output_dir` | Output directory |
| `--apadb` | APAdb database file (**mutually exclusive with `-f`**, for APA annotation using APAdb database) |
| `-f, --fasta` | Reference genome fasta file (**mutually exclusive with `--apadb`**, for de novo APA identification) |
| `-d, --distance` | APA site merging window size (default: 50bp, only effective when using `-f`) |
| `-t, --threads` | Number of threads (default: 4, maximum: 44, only effective when using `-f`) |

**Differences between the two modes**:

- **APAdb mode** (`--apadb`): Uses known APA site database for annotation, fast, depends on database quality
- **De novo identification mode** (`-f`): Identifies PAS sites de novo by detecting poly(A) tail signals at read ends, then clusters them into APA intervals within specified windows, independent of external databases

**APAdb data source**:

Recommended to download from [APAdb](https://ngdc.cncb.ac.cn/databasecommons/database/id/853), or use custom APA site data.

**APAdb file format requirements** (BED format, 10 columns):
```
chr1	16442	16450	WASH7P.1:16442-16450	43	-	Intron	16443.0	16443	-
chr1	134934	134953	LOC729737.1:134934-134953	26	-	UTR3	134944.0	134944	-
```

**Important**: APAdb file needs to be sorted:
```bash
sort -k1,1 -k2,2n apa_raw.bed > apa_sorted.bed
```

**Output files**:

- `<prefix>_APA_result.csv`: APA usage information

**Output column description**:

| Column name | Description |
|------|------|
| `readID` | Read identifier |
| `APA_type` | APA site type |
| `PAS_site` | PAS site position (de novo identification mode only) |

**Example**:

Using APAdb:
```bash
nanornaqtl pheno APA \
  -b sample01_calls_sorted_map.bam \
  -o ./work_dir \
  -p sample01_APA \
  --apadb apa_sorted.bed
```

De novo identification:
```bash
nanornaqtl pheno APA \
  -b sample01_calls_sorted_map.bam \
  -o ./work_dir \
  -p sample01_APA \
  -f hg38.fa \
  -d 35 \
  -t 20
```

---

#### 2.8 isoform Identification

**Function**: Identify and quantify transcript isoforms.

**Tool**: Uses [IsoQuant](https://github.com/ablab/IsoQuant) tool for isoform analysis.

**Command example**:

```bash
isoquant.py \
  --reference hg38.fa \
  --genedb annotation.gtf \
  --bam map.bam \
  --data_type nanopore \
  -o ./work_dir \
  -t 60 \
  --complete_genedb
```

**Parameter explanation**:

| Parameter | Description |
|------|------|
| `--reference` | Reference genome FASTA file |
| `--genedb` | Gene annotation GTF file |
| `--bam` | Input BAM file |
| `--data_type` | Data type (nanopore/pacbio_ccs) |
| `-o` | Output directory |
| `-t` | Number of threads |
| `--complete_genedb` | Use complete gene database mode |

**Output files**:

IsoQuant will generate multiple output files, the key file is:
- `OUT.transcript_model_reads.tsv.gz`: Read-to-isoform mapping relationship

This file can be directly used for subsequent QTL analysis.

---

### 3. qtl - QTL Analysis

All QTL analysis modules require a **SNP information file** in the following format:

**SNP information file format** (`snp_info.txt`):

```
chrom    pos    A1    A2    AF
chr1    123456    A    G    0.35
chr1    789012    C    T    0.42
```

**Column description**:

| Column name | Description |
|------|------|
| `chrom` | Chromosome name (must match BAM file) |
| `pos` | SNP position (1-based) |
| `A1` | Alternative allele |
| `A2` | Reference allele |
| `AF` | Alternative allele frequency |

**Genome size file format** (`genome.chrom.sizes`):

```
chr1    249250621
chr2    243199373
chr3    198022430
```

**Format description**: Tab-separated, two columns: chromosome name and chromosome length

---

#### 3.1 m6A QTL Analysis

**Function**: Identify genetic variants associated with m6A modification.

**Command**:

```bash
nanornaqtl qtl m6A \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --modification_csv <m6A_sites_csv> \
  --read_mod_dict <m6A_reads_pkl> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `-b, --bam` | map.bam file path from prep step | Required |
| `--snp_info` | SNP information file path | Required |
| `-p, --prefix` | Output file prefix name | Required |
| `-o, --output_dir` | Output directory path, will be created if not exists | Required |
| `--modification_csv` | m6A site result file (`*_m6A_sites_result.csv`) | Required |
| `--read_mod_dict` | m6A read dictionary file (`*_m6A_reads_final.pkl`) | Required |
| `--geno_size` | Genome size file | Required |
| `-q, --min_qscore` | Minimum base quality score for variant site | 10 |
| `-c, --min_coverage` | Minimum coverage for variant site | 8 |
| `--mcmc_samples` | Number of MCMC samples | 2000 |
| `-t, --threads` | Number of threads | 4 |
| `--keep_tmp` | Keep temporary files | False |

**Output files**:

- `<prefix>_m6A_QTLs_result.csv`: m6A QTL results

**Output column description**:

| Column name | Description |
|------|------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `mod_pos_1base` | Modification site position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `A1_mod`, `A1_unmod` | Number of modified/unmodified reads for A1 |
| `A2_mod`, `A2_unmod` | Number of modified/unmodified reads for A2 |
| `BF` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `fisher_pvalue` | Fisher's exact test p-value |

**Statistical method explanation**:

1. **Bayesian method** (recommended): BF and posterior probability
2. **Frequentist method** (reference): Fisher's exact test

**Example**:

```bash
nanornaqtl qtl m6A \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_m6A_qtl \
  -o ./work_dir \
  --modification_csv sample01_m6A_sites_result.csv \
  --read_mod_dict sample01_m6A_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.2 m5C QTL Analysis

**Function**: Identify genetic variants associated with m5C modification.

**Command**:

```bash
nanornaqtl qtl m5C \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --modification_csv <m5C_sites_csv> \
  --read_mod_dict <m5C_reads_pkl> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameters**: Same as m6A QTL

**Output files**:

- `<prefix>_m5C_QTLs_result.csv`: m5C QTL results

**Output columns**: Same as m6A QTL

**Example**:

```bash
nanornaqtl qtl m5C \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_m5C_qtl \
  -o ./work_dir \
  --modification_csv sample01_m5C_sites_result.csv \
  --read_mod_dict sample01_m5C_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.3 pseudouridine (pseU) QTL Analysis

**Function**: Identify genetic variants associated with pseudouridine modification.

**Command**:

```bash
nanornaqtl qtl pseU \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --modification_csv <pseU_sites_csv> \
  --read_mod_dict <pseU_reads_pkl> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameters**: Same as m6A QTL

**Output files**:

- `<prefix>_pseU_QTLs_result.csv`: pseU QTL results

**Example**:

```bash
nanornaqtl qtl pseU \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_pseU_qtl \
  -o ./work_dir \
  --modification_csv sample01_pseU_sites_result.csv \
  --read_mod_dict sample01_pseU_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.4 inosine QTL Analysis

**Function**: Identify genetic variants associated with inosine modification.

**Command**:

```bash
nanornaqtl qtl inosine \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --modification_csv <inosine_sites_csv> \
  --read_mod_dict <inosine_reads_pkl> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameters**: Same as m6A QTL

**Output files**:

- `<prefix>_inosine_QTLs_result.csv`: inosine QTL results

**Example**:

```bash
nanornaqtl qtl inosine \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_inosine_qtl \
  -o ./work_dir \
  --modification_csv sample01_inosine_sites_result.csv \
  --read_mod_dict sample01_inosine_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.5 APA QTL Analysis

**Function**: Identify genetic variants affecting alternative polyadenylation site usage.

**Command**:

```bash
nanornaqtl qtl APA \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --read_overlap_file <apa_result_csv> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `-b, --bam` | map.bam file path from prep step | Required |
| `--snp_info` | SNP information file path | Required |
| `-p, --prefix` | Output file prefix name | Required |
| `-o, --output_dir` | Output directory path, will be created if not exists | Required |
| `--read_overlap_file` | APA result file (`*_APA_result.csv`) | Required |
| `--geno_size` | Genome size file | Required |
| `-q, --min_qscore` | Minimum base quality score for variant site | 10 |
| `-c, --min_coverage` | Minimum coverage for variant site | 8 |
| `--mcmc_samples` | Number of MCMC samples | 2000 |
| `-t, --threads` | Number of threads | 4 |
| `--keep_tmp` | Keep temporary files | False |

**Output files**:

- `<prefix>_APA_QTLs_result.csv`: APA QTL results

**Output column description**:

| Column name | Description |
|------|------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `BF` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `chi2_pvalue` | Chi-square test p-value |
| `TVD` | Total Variation Distance (range 0-1) |
| `dominant_shift` | Main APA type usage change direction and magnitude |

**Statistical method explanation**:

1. **Bayesian method** (recommended): BF and posterior probability
2. **Frequentist method** (reference): Chi-square test

**Effect size indicator explanation**:

- **TVD** (Total Variation Distance): Measures the degree of difference in APA usage distribution between A1 and A2 alleles
  - 0: Two distributions are completely identical
  - 1: Two distributions are completely different
  
- **dominant_shift**: Describes the change in usage proportion of each APA type between A1 and A2
  - Format example: `type1:A1↓(-0.56)` indicates reads carrying A1 allele use this APA type less (56% relative decrease)

**Example**:

```bash
nanornaqtl qtl APA \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_APA_qtl \
  -o ./work_dir \
  --read_overlap_file sample01_APA_result.csv \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.6 isoform QTL Analysis

**Function**: Identify genetic variants affecting transcript isoform usage.

**Command**:

```bash
nanornaqtl qtl isoform \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --read_overlap_file <isoquant_output> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description |
|------|------|
| `--read_overlap_file` | IsoQuant output file (`OUT.transcript_model_reads.tsv.gz`) |
| Other parameters | Same as APA QTL |

**Output files**:

- `<prefix>_isoform_QTLs_result.csv`: isoform QTL results

**Output columns**: Same as APA QTL, where:
- `chi2_pvalue`: Chi-square test p-value
- `TVD` and `dominant_shift` describe isoform usage pattern differences

**Example**:

```bash
nanornaqtl qtl isoform \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_isoform_qtl \
  -o ./work_dir \
  --read_overlap_file OUT.transcript_model_reads.tsv.gz \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.7 polyA Tail Length QTL Analysis

**Function**: Identify genetic variants affecting polyA tail length.

**Command**:

```bash
nanornaqtl qtl polyA_tail \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --polya_csv <polyA_csv> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `-b, --bam` | map.bam file path from prep step | Required |
| `--snp_info` | SNP information file path | Required |
| `-p, --prefix` | Output file prefix name | Required |
| `--polya_csv` | polyA tail length result CSV file path (`*_polyAlen_result.csv`) | Required |
| `--geno_size` | Genome size file | Required |
| `-q, --min_qscore` | Minimum base quality score for variant site | 10 |
| `-c, --min_coverage` | Minimum coverage for variant site | 8 |
| `--mcmc_samples` | Number of MCMC samples | 2000 |
| `-t, --threads` | Number of threads | 4 |
| `--keep_tmp` | Keep temporary files | False |

**Output files**:

- `<prefix>_polyA_tail_length_QTLs_result.csv`: polyA tail length QTL results

**Output column description**:

| Column name | Description |
|------|------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `A1_len` | polyA tail length list for A1 allele |
| `A2_len` | polyA tail length list for A2 allele |
| `beta` | Effect size (log ratio of mean lengths) |
| `SE` | Standard error |
| `KS_stat` | Kolmogorov-Smirnov statistic |
| `bayes_factor` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `p_welch` | Welch's t-test p-value |
| `p_mw` | Mann-Whitney U test p-value |
| `p_ks` | Kolmogorov-Smirnov test p-value |

**Frequentist method explanation** (three p-values):

- **p_welch**: Welch's t-test, suitable for unequal variances
- **p_mw**: Mann-Whitney U test, non-parametric test, robust to outliers
- **p_ks**: Kolmogorov-Smirnov test, tests distribution differences

**Example**:

```bash
nanornaqtl qtl polyA_tail \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_polyA_qtl \
  -o ./work_dir \
  --polya_csv sample01_polyAlen_result.csv \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

#### 3.8 Intron Retention Rate QTL Analysis

**Function**: Identify genetic variants affecting intron retention rate.

**Command**:

```bash
nanornaqtl qtl intron_retention \
  -b <map_bam> \
  --snp_info <snp_file> \
  -p <output_prefix> \
  -o <output_dir> \
  --ir_csv <ir_csv> \
  --geno_size <genome_size_file> \
  -t <threads>
```

**Parameter explanation**:

| Parameter | Description | Default |
|------|------|--------|
| `--ir_csv` | Intron retention rate result file (`*_intronRetention_result.csv`) | Required |
| Other parameters | Same as polyA tail length QTL |

**Output files**:

- `<prefix>_intron_retention_QTLs_result.csv`: Intron retention rate QTL results

**Output column description**:

| Column name | Description |
|------|------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `A1_IR` | Intron retention rate list for A1 allele |
| `A2_IR` | Intron retention rate list for A2 allele |
| `beta` | Effect size (logit difference) |
| `SE` | Standard error |
| `KS_stat` | Kolmogorov-Smirnov statistic |
| `bayes_factor` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `p_welch` | Welch's t-test p-value |
| `p_mw` | Mann-Whitney U test p-value |
| `p_ks` | Kolmogorov-Smirnov test p-value |

**Note**: Automatically filters low-quality data (e.g., IR=0 for more than 90% of reads)

**Example**:

```bash
nanornaqtl qtl intron_retention \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -p sample01_IR_qtl \
  -o ./work_dir \
  --ir_csv sample01_intronRetention_result.csv \
  --geno_size hg19.chrom.sizes \
  -t 20
```

---

## Output File Description

### File Naming Rules

- **prep module**: `<prefix>_calls_sorted_map[0/16].bam`, `<prefix>_calls_sorted_map.fastq`
- **pheno module**: `<prefix>_<phenotype>_result.csv`
- **qtl module**: `<prefix>_<phenotype>_QTLs_result.csv`

### Common Output Format

All CSV output files are standard tab-separated or comma-separated files that can be opened and analyzed with Excel, R, Python, and other tools.

### Key Statistical Indicator Interpretation 📈

#### Bayesian Methods

- **BF (Bayes Factor)**:
  - BF > 3: Significant evidence supporting association
  - BF > 10: Strong evidence
  - BF > 30: Very strong evidence
  
- **Posterior probability**: Posterior probability of null hypothesis being true; smaller values indicate higher likelihood of association

#### Frequentist Methods

- **p-value**: Traditional significance test p-value, typically p < 0.05 is considered significant
- **Note**: Frequentist methods may produce more false positives with uneven coverage; Bayesian methods are recommended

---

## Citation

If you use nanornaqtl in your research, please cite our GitHub repository:

```
nanornaqtl: A comprehensive toolkit for molecular phenotyping and QTL analysis using Nanopore direct RNA sequencing
https://github.com/xinranxu0930/nanornaqtl
```

---

## Contact 📧

- **Issue reporting**: [GitHub Issues](https://github.com/xinranxu0930/nanornaqtl/issues)
- **Email**: xinranxu0930@gmail.com
- **GitHub**: [https://github.com/xinranxu0930/nanornaqtl](https://github.com/xinranxu0930/nanornaqtl)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

---

## Changelog

### v1.0.0 (2026-01-18)

- Initial release
- Support for 7 types of molecular phenotype identification
- Support for 8 types of QTL analysis
- Implementation of Bayesian statistical methods
- Support for multi-threading parallel processing

### v1.0.3 (2026-01-20)

- Fixed coverage calculation
- Added logging system
- Optimized base mode support issues
- Added version query

### v1.0.4 (2026-01-27)

- Added new de novo APA identification feature
- Fixed parameter description errors in README