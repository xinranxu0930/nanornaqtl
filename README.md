# nanornaqtl

[![PyPI version](https://badge.fury.io/py/nanornaqtl.svg)](https://badge.fury.io/py/nanornaqtl)
[![Python versions](https://img.shields.io/pypi/pyversions/nanornaqtl.svg)](https://pypi.org/project/nanornaqtl/)

A comprehensive toolkit for multi-molecular phenotyping and QTL analysis based on Nanopore direct RNA sequencing data


## Introduction

**nanornaqtl** is a specialized analysis tool designed for Nanopore direct RNA sequencing data, used to identify various RNA molecular phenotypes and perform population-level QTL (Quantitative Trait Loci) analysis. This tool is particularly suitable for analyzing pooling samples (mixed sequencing of multiple individuals), with all analyses conducted at the read level.

### Key Features

- **7 Molecular Phenotype Identification** 📊:
  - RNA modification sites: m6A, m5C, pseudouridine (pseU), inosine
  - polyA tail length
  - Intron retention rate
  - Alternative polyA sites (APA)
  - Transcript isoforms

- **8 QTL Analyses** 🧪:
  - Modification QTLs (m6A, m5C, pseU, inosine)
  - APA usage pattern QTL
  - Isoform usage pattern QTL
  - polyA tail length QTL
  - Intron retention rate QTL

### Features ✨

- **Read-level analysis**: Suitable for population-level analysis of pooling samples
- **Bayesian statistical methods**: More robust for unevenly covered Nanopore data
- **Parallel processing**: Supports multi-threading acceleration, processing by chromosome in parallel
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

#### Basic alignment (required)

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --mm2-opts "-x splice -k 14" \
  > output.bam
```

**Parameter description**:
- `--mm2-opts "-x splice -k 14"`: RNA splice alignment parameters

#### Optional features

**If modification site analysis is needed**:

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
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
  -b output.bam \
  -p output_prefix \
  -t 4 \
  -q 0

# 2. Identify m6A modification sites
nanornaqtl pheno m6A \
  -b output_prefix_calls_sorted_map.bam \
  -o result_prefix \
  -t 20 \
  --motif

# 3. Perform m6A QTL analysis
nanornaqtl qtl m6A \
  -b output_prefix_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o qtl_result_prefix \
  -csv result_prefix_m6A_sites_result.csv \
  -pkl result_prefix_m6A_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -m m6A \
  --threads 20
```

---

## Detailed Usage Instructions

### 1. prep - Data Preprocessing

**Function**: Filter unmapped reads from raw BAM file, generate FASTQ of mapped reads and strand-separated BAM files.

**Input**:
- Raw BAM file output by Dorado (no preprocessing needed)

**Command**:

```bash
nanornaqtl prep \
  -b <basecall_bam> \
  -p <output_prefix> \
  -t <threads> \
  -q <min_mapq>
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-b, --basecall_bam` | Input raw BAM file path | Required |
| `-p, --dir_pre` | Output file prefix | Required |
| `-t, --threads` | Number of threads | 4 |
| `-q, --min_mapq` | Minimum MAPQ threshold | 0 |

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
  -o <output_prefix> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-b, --bam` | Input BAM file(`*_map.bam`) | Required |
| `-o, --output_prefix` | Output file prefix | Required |
| `-t, --threads` | Number of threads (max effective value 22) | 4 |
| `-f, --mod_threshold` | Modification probability threshold | 0.75 |
| `-q, --min_qscore` | Minimum base quality | 10 |
| `--min_mapq` | Minimum MAPQ | 0 |
| `-r, --min_rate` | Minimum modification rate | 0.1 |
| `-c, --min_cov` | Minimum coverage | 5 |
| `--motif` | Enable motif filtering (DRACH: `[GAT][GA]AC[ATC]`) | False |
| `--metaPlotR` | Generate bed file for [metaPlotR](https://github.com/olarerin/metaPlotR) | False |

**Output files**:

- `<prefix>_m6A_sites_result.csv`: m6A modification site information

**Output column description**:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `pos_1base` | Position (1-based) |
| `strand` | Strand (+/-) |
| `mod_num` | Number of modified reads |
| `cov` | Total coverage |
| `mod_rate` | Modification rate |
| `motif` | Motif sequence (if --motif enabled) |

**Example**:

```bash
# Using motif filtering
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -o sample01_m6A \
  -t 20 \
  --motif \
  --metaPlotR

# Without motif filtering
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -o sample01_m6A \
  -t 20
```

---

#### 2.2 m5C Modification Site Identification

**Function**: Identify m5C (5-methylcytosine) modification sites.

**Command**:

```bash
nanornaqtl pheno m5C \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
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
  -o sample01_m5C \
  -t 20 \
  --motif
```

---

#### 2.3 Pseudouridine (pseU) Modification Site Identification

**Function**: Identify pseudouridine modification sites.

**Command**:

```bash
nanornaqtl pheno pseU \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
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
  -o sample01_pseU \
  -t 20 \
  --motif
```

---

#### 2.4 Inosine Modification Site Identification

**Function**: Identify inosine (A-to-I editing) sites.

**Command**:

```bash
nanornaqtl pheno inosine \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**Parameters**: Same as m6A

**Motif** (if --motif enabled): `TA[GT]`

**Output files**:

- `<prefix>_inosine_sites_result.csv`: Inosine site information

**Example**:

```bash
nanornaqtl pheno inosine \
  -b sample01_calls_sorted_map.bam \
  -o sample01_inosine \
  -t 20 \
  --motif
```

---

#### 2.5 PolyA Tail Length Identification

**Function**: Extract polyA tail length for each read.

**Prerequisites**: BAM file must contain pt tag (generated by Dorado basecall with `--estimate-poly-a`)

**Command**:

```bash
nanornaqtl pheno polyA_tail \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads>
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-b, --bamfile` | Input BAM file | Required |
| `-o, --output_prefix` | Output prefix | Required |
| `-t, --threads` | Number of threads (max effective value 22) | 4 |

**Output files**:

- `<prefix>_polyAlen_result.csv`: PolyA tail length information

**Output column description**:

| Column | Description |
|--------|-------------|
| `readID` | Read identifier |
| `polyA_length` | PolyA tail length |

**Example**:

```bash
nanornaqtl pheno polyA_tail \
  -b sample01_calls_sorted_map.bam \
  -o sample01_polyA \
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

**Parameter description**:

| Parameter | Description |
|-----------|-------------|
| `-g, --gtf` | GTF annotation file path |
| `-b, --bam` | Input BAM file |
| `-o, --output` | Output directory |
| `-p, --output_prefix` | Output file prefix |

**Output files**:

- `<prefix>_intronRetention_result.csv`: Intron retention rate information

**Output column description**:

| Column | Description |
|--------|-------------|
| `readID` | Read identifier |
| `intron_retention_rate` | Intron retention rate |

**Example**:

```bash
nanornaqtl pheno intron_retention \
  -g gencode.v38.annotation.gtf \
  -b sample01_calls_sorted_map.bam \
  -o output_dir \
  -p sample01_IR
```

---

#### 2.7 Alternative PolyA Site (APA) Identification

**Function**: Identify alternative polyA sites and classify read usage patterns.

**Command**:

```bash
nanornaqtl pheno APA \
  -g <gtf_file> \
  -b <map_bam> \
  -o <output_dir> \
  -p <output_prefix>
```

**Parameter description**:

| Parameter | Description |
|-----------|-------------|
| `-g, --gtf` | GTF annotation file path |
| `-b, --bam` | Input BAM file |
| `-o, --output` | Output directory |
| `-p, --output_prefix` | Output file prefix |

**Output files**:

- `<prefix>_APA_result.csv`: APA site information

**Output column description**:

| Column | Description |
|--------|-------------|
| `readID` | Read identifier |
| `APA_type` | APA type classification |

**Example**:

```bash
nanornaqtl pheno APA \
  -g gencode.v38.annotation.gtf \
  -b sample01_calls_sorted_map.bam \
  -o output_dir \
  -p sample01_APA
```

---

### 3. qtl - QTL Analysis

**Prerequisite**: All QTL analyses require a SNP information file.

#### SNP Information File Format

The SNP information file must be a tab-separated text file containing the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| `chrom` | Chromosome | chr1 |
| `pos` | Position (1-based) | 12345 |
| `SNP` | SNP identifier | rs123456 |
| `A1` | Allele 1 | A |
| `A2` | Allele 2 | G |

**Example**:

```
chrom	pos	SNP	A1	A2
chr1	12345	rs123456	A	G
chr1	67890	rs789012	C	T
```

---

#### 3.1 m6A QTL Analysis

**Function**: Identify genetic variants associated with m6A modification.

**Command**:

```bash
nanornaqtl qtl m6A \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <m6A_sites_csv> \
  -pkl <m6A_reads_pkl> \
  --geno_size <genome_size_file> \
  -m m6A \
  --threads <threads>
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-b, --bam` | Input BAM file | Required |
| `--snp_info` | Variant site file | Required |
| `-o, --output_prefix` | Output prefix | Required |
| `-csv, --m6A_sites` | m6A site result file (`*_sites_result.csv`) | Required |
| `-pkl, --m6A_reads` | m6A read information file (`*_reads_final.pkl`) | Required |
| `--geno_size` | Genome size file | Required |
| `-m, --molecular_type` | Molecular phenotype type (m6A) | Required |
| `-q, --min_qscore` | Minimum base quality | 10 |
| `--min_coverage` | Minimum coverage | 8 |
| `--mcmc_samples` | MCMC sampling number | 1000 |
| `--threads` | Number of threads | 4 |
| `--keep_tmp` | Keep temporary files | False |

**Output files**:

- `<prefix>_m6A_QTLs_result.csv`: m6A QTL results

**Output column description**:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `mod_site_pos_1base` | Modification site position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `A1_mod`, `A1_unmod` | A1 allele modification counts |
| `A2_mod`, `A2_unmod` | A2 allele modification counts |
| `bayes_factor` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `p_fisher` | Fisher's exact test p-value |

**Statistical method description**:

1. **Bayesian method** (recommended):
   - **Bayes Factor (BF)**: Quantifies evidence for association
   - **Posterior probability**: Probability that null hypothesis is true
   
2. **Frequentist method** (reference):
   - **Fisher's exact test**: More aggressive, higher false positive rate

**Note**: Output includes all statistical results (significant and non-significant), users can filter based on their own thresholds (e.g., BF > 3).

**Example**:

```bash
nanornaqtl qtl m6A \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_m6A_qtl \
  -csv sample01_m6A_sites_result.csv \
  -pkl sample01_m6A_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -m m6A \
  --threads 20 \
  --mcmc_samples 2000
```

---

#### 3.2 m5C QTL Analysis

**Function**: Identify genetic variants associated with m5C modification.

**Command**:

```bash
nanornaqtl qtl m5C \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <m5C_sites_csv> \
  -pkl <m5C_reads_pkl> \
  --geno_size <genome_size_file> \
  -m m5C \
  --threads <threads>
```

**Parameters**: Same as m6A QTL, change `-m` parameter to `m5C`

**Output files**:

- `<prefix>_m5C_QTLs_result.csv`: m5C QTL results

**Output columns**: Same as m6A QTL

**Example**:

```bash
nanornaqtl qtl m5C \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_m5C_qtl \
  -csv sample01_m5C_sites_result.csv \
  -pkl sample01_m5C_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -m m5C \
  --threads 20
```

---

#### 3.3 Pseudouridine (pseU) QTL Analysis

**Function**: Identify genetic variants associated with pseudouridine modification.

**Command**:

```bash
nanornaqtl qtl pseU \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <pseU_sites_csv> \
  -pkl <pseU_reads_pkl> \
  --geno_size <genome_size_file> \
  -m pseU \
  --threads <threads>
```

**Parameters**: Same as m6A QTL, change `-m` parameter to `pseU`

**Output files**:

- `<prefix>_pseU_QTLs_result.csv`: pseU QTL results

**Example**:

```bash
nanornaqtl qtl pseU \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_pseU_qtl \
  -csv sample01_pseU_sites_result.csv \
  -pkl sample01_pseU_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -m pseU \
  --threads 20
```

---

#### 3.4 Inosine QTL Analysis

**Function**: Identify genetic variants associated with inosine modification.

**Command**:

```bash
nanornaqtl qtl inosine \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <inosine_sites_csv> \
  -pkl <inosine_reads_pkl> \
  --geno_size <genome_size_file> \
  -m inosine \
  --threads <threads>
```

**Parameters**: Same as m6A QTL, change `-m` parameter to `inosine`

**Output files**:

- `<prefix>_inosine_QTLs_result.csv`: Inosine QTL results

**Example**:

```bash
nanornaqtl qtl inosine \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_inosine_qtl \
  -csv sample01_inosine_sites_result.csv \
  -pkl sample01_inosine_reads_final.pkl \
  --geno_size hg19.chrom.sizes \
  -m inosine \
  --threads 20
```

---

#### 3.5 APA QTL Analysis

**Function**: Identify genetic variants affecting alternative polyA site usage.

**Command**:

```bash
nanornaqtl qtl APA \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -f <apa_result_csv> \
  --geno_size <genome_size_file> \
  -m APA \
  -t <threads>
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-b, --bam` | Input BAM file | Required |
| `--snp_info` | Variant site file | Required |
| `-o, --output_prefix` | Output prefix | Required |
| `-f, --read_overlap_file` | APA result file (`*_APA_result.csv`) | Required |
| `--geno_size` | Genome size file | Required |
| `-m, --molecular_type` | Molecular phenotype type (APA) | Required |
| `-q, --min_qscore` | Minimum base quality | 10 |
| `--min_coverage` | Minimum coverage | 8 |
| `--mcmc_samples` | MCMC sampling number | 1000 |
| `-t, --threads` | Number of threads | 4 |
| `--keep_tmp` | Keep temporary files | False |

**Output files**:

- `<prefix>_APA_QTLs_result.csv`: APA QTL results

**Output column description**:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `BF` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `chi2_pvalue` | Chi-square test p-value |
| `TVD` | Total Variation Distance (range 0-1) |
| `dominant_shift` | Direction and magnitude of major APA type usage change |

**Statistical method description**:

1. **Bayesian method** (recommended): BF and posterior probability
2. **Frequentist method** (reference): Chi-square test

**Effect size metric description**:

- **TVD** (Total Variation Distance): Measures the degree of difference in APA usage distribution between A1 and A2 alleles
  - 0: Two distributions are identical
  - 1: Two distributions are completely different
  
- **dominant_shift**: Describes the change in usage proportion of each APA type between A1 and A2
  - Example format: `type1:A1↓(-0.56)` indicates reads carrying A1 allele use this APA type less frequently (relative decrease of 56%)

**Example**:

```bash
nanornaqtl qtl APA \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_APA_qtl \
  -f sample01_APA_result.csv \
  --geno_size hg19.chrom.sizes \
  -m APA \
  -t 20
```

---

#### 3.6 Isoform QTL Analysis

**Function**: Identify genetic variants affecting transcript isoform usage.

**Command**:

```bash
nanornaqtl qtl isoform \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -f <isoquant_output> \
  --geno_size <genome_size_file> \
  -m isoform \
  -t <threads>
```

**Parameter description**:

| Parameter | Description |
|-----------|-------------|
| `-f, --read_overlap_file` | IsoQuant output file (`OUT.transcript_model_reads.tsv.gz`) |
| Other parameters | Same as APA QTL |

**Output files**:

- `<prefix>_isoform_QTLs_result.csv`: Isoform QTL results

**Output columns**: Same as APA QTL, where:
- `chi2_pvalue`: Chi-square test p-value
- `TVD` and `dominant_shift` describe isoform usage pattern differences

**Example**:

```bash
nanornaqtl qtl isoform \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_isoform_qtl \
  -f OUT.transcript_model_reads.tsv.gz \
  --geno_size hg19.chrom.sizes \
  -m isoform \
  -t 20
```

---

#### 3.7 PolyA Tail Length QTL Analysis

**Function**: Identify genetic variants affecting polyA tail length.

**Command**:

```bash
nanornaqtl qtl polyA_tail \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <polyA_csv> \
  --geno_size <genome_size_file> \
  --threads <threads>
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-b, --bam` | Input BAM file | Required |
| `--snp_info` | Variant site file | Required |
| `-o, --output_prefix` | Output prefix | Required |
| `-csv, --polya_csv` | PolyA result file (`*_polyAlen_result.csv`) | Required |
| `--geno_size` | Genome size file | Required |
| `-q, --min_qscore` | Minimum base quality | 10 |
| `--min_coverage` | Minimum coverage | 8 |
| `--mcmc_samples` | MCMC sampling number | 2000 |
| `--threads` | Number of threads | 4 |
| `--keep_tmp` | Keep temporary files | False |

**Output files**:

- `<prefix>_polyA_tail_length_QTLs_result.csv`: PolyA tail length QTL results

**Output column description**:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `SNP` | Variant ID |
| `snp_pos_1base` | SNP position |
| `A1`, `A2` | Alleles |
| `MAF` | Minor allele frequency |
| `A1_len` | PolyA tail length list for A1 allele |
| `A2_len` | PolyA tail length list for A2 allele |
| `beta` | Effect size (log ratio of mean lengths) |
| `SE` | Standard error |
| `KS_stat` | Kolmogorov-Smirnov statistic |
| `bayes_factor` | Bayes Factor |
| `posterior_prob` | Posterior probability |
| `p_welch` | Welch's t-test p-value |
| `p_mw` | Mann-Whitney U test p-value |
| `p_ks` | Kolmogorov-Smirnov test p-value |

**Frequentist method description** (three p-values):

- **p_welch**: Welch's t-test, suitable for unequal variances
- **p_mw**: Mann-Whitney U test, non-parametric test, robust to outliers
- **p_ks**: Kolmogorov-Smirnov test, tests distribution differences

**Example**:

```bash
nanornaqtl qtl polyA_tail \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_polyA_qtl \
  -csv sample01_polyAlen_result.csv \
  --geno_size hg19.chrom.sizes \
  --threads 20 \
  --mcmc_samples 2000
```

---

#### 3.8 Intron Retention Rate QTL Analysis

**Function**: Identify genetic variants affecting intron retention rate.

**Command**:

```bash
nanornaqtl qtl intron_retention \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <ir_csv> \
  --geno_size <genome_size_file> \
  --threads <threads>
```

**Parameter description**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-csv, --ir_csv` | Intron retention rate result file (`*_intronRetention_result.csv`) | Required |
| Other parameters | Same as polyA tail length QTL |

**Output files**:

- `<prefix>_intron_retention_QTLs_result.csv`: Intron retention rate QTL results

**Output column description**:

| Column | Description |
|--------|-------------|
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

**Note**: Automatically filters low-quality data (e.g., >90% reads with IR=0)

**Example**:

```bash
nanornaqtl qtl intron_retention \
  -b sample01_calls_sorted_map.bam \
  --snp_info snp_info.txt \
  -o sample01_IR_qtl \
  -csv sample01_intronRetention_result.csv \
  --geno_size hg19.chrom.sizes \
  --threads 20
```

---

## Output File Description

### File Naming Rules

- **prep module**: `<prefix>_calls_sorted_map[0/16].bam`, `<prefix>_calls_sorted_map.fastq`
- **pheno module**: `<prefix>_<phenotype>_result.csv`
- **qtl module**: `<prefix>_<phenotype>_QTLs_result.csv`

### General Output Format

All CSV output files are standard tab-separated or comma-separated files that can be opened and analyzed using Excel, R, Python, and other tools.

### Key Statistical Metrics Interpretation 📈

#### Bayesian Methods

- **BF (Bayes Factor)**:
  - BF > 3: Substantial evidence supporting association
  - BF > 10: Strong evidence
  - BF > 30: Very strong evidence
  
- **Posterior probability**: Posterior probability that null hypothesis is true, smaller values indicate stronger evidence for association

#### Frequentist Methods

- **p-value**: Traditional significance test p-value, typically p < 0.05 considered significant
- **Note**: Frequentist methods may produce more false positives with uneven coverage, Bayesian methods recommended as primary approach

---

## Citation

If you use nanornaqtl in your research, please cite our GitHub repository:

```
nanornaqtl: A comprehensive toolkit for molecular phenotyping and QTL analysis using Nanopore direct RNA sequencing
https://github.com/xinranxu0930/nanornaqtl
```

---

## Contact 📧

- **Issue Reporting**: [GitHub Issues](https://github.com/xinranxu0930/nanornaqtl/issues)
- **Email**: xinranxu0930@gmail.com
- **GitHub**: [https://github.com/xinranxu0930/nanornaqtl](https://github.com/xinranxu0930/nanornaqtl)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

---

## Changelog

### v1.0.0 (2026-01-18)

- Initial release
- Support for 7 molecular phenotype identification methods
- Support for 8 QTL analysis types
- Implementation of Bayesian statistical methods
- Support for multi-threaded parallel processing

### v1.0.3 (2026-01-20)

- Fix the coverage calculation
- Add a logging system
- Optimize the base mode support issue
- Add version query