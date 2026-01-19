# nanornaqtl

[![PyPI version](https://badge.fury.io/py/nanornaqtl.svg)](https://badge.fury.io/py/nanornaqtl)
[![Python versions](https://img.shields.io/pypi/pyversions/nanornaqtl.svg)](https://pypi.org/project/nanornaqtl/)

基于Nanopore direct RNA测序数据的多分子表型识别与QTL分析工具包


## 简介

**nanornaqtl** 是一个专门为Nanopore direct RNA测序数据设计的分析工具,用于识别多种RNA分子表型并进行群体水平的QTL (Quantitative Trait Loci) 分析。该工具特别适用于pooling样本(多个个体混合测序)的分析,所有分析均在read水平进行。

### 主要功能

- **7种分子表型识别** 📊:
  - RNA修饰位点: m6A, m5C, pseudouridine (pseU), inosine
  - polyA尾长
  - 内含子滞留率
  - 可变polyA位点 (APA)
  - 转录本异构体 (isoform)

- **8种QTL分析** 🧪:
  - 修饰QTL (m6A, m5C, pseU, inosine)
  - APA使用模式QTL
  - isoform使用模式QTL
  - polyA尾长QTL
  - 内含子滞留率QTL

### 特点 ✨

- **read水平分析**: 适用于pooling样本的群体水平分析
- **贝叶斯统计方法**: 对覆盖度不均匀的Nanopore数据更稳健
- **并行处理**: 支持多线程加速,按染色体并行处理
- **灵活的参数设置**: 可自定义质量阈值、覆盖度要求等

---

## 安装

### 通过pip安装 📦

```bash
pip install nanornaqtl
```

### 从源码安装

```bash
git clone https://github.com/xinranxu0930/nanornaqtl.git
cd nanornaqtl
poetry install
```

### 依赖软件

#### 外部工具:

- **bedtools** (v2.30.0+)
- **samtools** (v1.15+)
- **IsoQuant** (用于isoform分析): [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant)

#### Python包:

- pysam
- pymc
- arviz
- statsmodels
- pandas
- numpy
- scipy

---

## 数据预处理

在使用nanornaqtl之前,需要对Nanopore下机数据进行预处理。

### 步骤1: 数据格式转换

如果下机数据是FAST5格式,需要先转换为POD5格式(如果已是POD5可跳过此步):

```bash
# 使用pod5工具转换
pod5 convert fast5 input.fast5 --output output.pod5
```

### 步骤2: Basecalling和比对 🔧

使用[Dorado](https://github.com/nanoporetech/dorado)进行basecalling和比对:

#### 基础比对 (必需)

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --mm2-opts "-x splice -k 14" \
  > output.bam
```

**参数说明**:
- `--mm2-opts "-x splice -k 14"`: RNA splice比对参数

#### 可选功能

**如果需要进行修饰位点分析**:

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --mm2-opts "-x splice -k 14" \
  --modified-bases-models model_name \
  > output.bam
```

下载修饰识别模型:
```bash
# 下载所有模型
dorado download --model all

# 或下载特定模型
dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.1.0
```

**如果需要进行polyA尾长分析**:

```bash
dorado basecaller \
  model_path \
  pod5_directory/ \
  --mm2-opts "-x splice -k 14" \
  --estimate-poly-a \
  > output.bam
```

---

## 快速开始

以下是一个完整的分析流程示例:

```bash
# 1. 数据预处理
nanornaqtl prep \
  -b output.bam \
  -p output_prefix \
  -t 4 \
  -q 0

# 2. 识别m6A修饰位点
nanornaqtl pheno m6A \
  -b output_prefix_calls_sorted_map.bam \
  -o result_prefix \
  -t 20 \
  --motif

# 3. 进行m6A QTL分析
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

## 详细使用说明

### 1. prep - 数据预处理

**功能**: 从原始BAM文件中过滤未比对reads,生成比对reads的FASTQ和分链BAM文件。

**输入**:
- Dorado输出的原始BAM文件(无需任何预处理)

**命令**:

```bash
nanornaqtl prep \
  -b <basecall_bam> \
  -p <output_prefix> \
  -t <threads> \
  -q <min_mapq>
```

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-b, --basecall_bam` | 输入的原始BAM文件路径 | 必需 |
| `-p, --dir_pre` | 输出文件前缀 | 必需 |
| `-t, --threads` | 线程数 | 4 |
| `-q, --min_mapq` | 最小MAPQ阈值 | 0 |

**输出文件**:

- `<prefix>_calls_sorted_map.bam`: 比对reads的合并BAM
- `<prefix>_calls_sorted_map0.bam`: 正链BAM (flag=0)
- `<prefix>_calls_sorted_map16.bam`: 负链BAM (flag=16)
- `<prefix>_calls_sorted_map.fastq`: 比对reads的FASTQ
- 对应的`.bai`索引文件

**示例**:

```bash
nanornaqtl prep \
  -b raw_output.bam \
  -p sample01 \
  -t 8 \
  -q 10
```

---

### 2. pheno - 分子表型识别

#### 2.1 m6A修饰位点识别

**功能**: 识别m6A (N6-methyladenosine)修饰位点。

**命令**:

```bash
nanornaqtl pheno m6A \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-b, --bam` | 输入BAM文件(`*_map.bam`) | 必需 |
| `-o, --output_prefix` | 输出文件前缀 | 必需 |
| `-t, --threads` | 线程数(最大有效值22) | 4 |
| `-f, --mod_threshold` | 修饰概率阈值 | 0.75 |
| `-q, --min_qscore` | 最小碱基质量 | 10 |
| `--min_mapq` | 最小MAPQ | 0 |
| `-r, --min_rate` | 最小修饰率 | 0.1 |
| `-c, --min_cov` | 最小覆盖度 | 5 |
| `--motif` | 启用motif过滤(DRACH: `[GAT][GA]AC[ATC]`) | False |
| `--metaPlotR` | 生成用于[metaPlotR](https://github.com/olarerin/metaPlotR)的bed文件 | False |

**输出文件**:

- `<prefix>_m6A_sites_result.csv`: m6A修饰位点信息

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `chrom` | 染色体 |
| `pos_1base` | 位置(1-based) |
| `strand` | 链方向(+/-) |
| `mod_num` | 修饰reads数 |
| `cov` | 总覆盖度 |
| `mod_rate` | 修饰率 |
| `motif` | motif序列(如启用--motif) |

**示例**:

```bash
# 使用motif过滤
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -o sample01_m6A \
  -t 20 \
  --motif \
  --metaPlotR

# 不使用motif过滤
nanornaqtl pheno m6A \
  -b sample01_calls_sorted_map.bam \
  -o sample01_m6A \
  -t 20
```

---

#### 2.2 m5C修饰位点识别

**功能**: 识别m5C (5-methylcytosine)修饰位点。

**命令**:

```bash
nanornaqtl pheno m5C \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**参数**: 与m6A相同

**motif分类** (如启用--motif):
- **CG**: CpG位点
- **CHG**: CHG位点
- **CHH**: CHH位点

**输出文件**:

- `<prefix>_m5C_sites_result.csv`: m5C修饰位点信息

**输出列说明**: 与m6A相同,额外包含`motif_classification`列(CG/CHG/CHH)

**示例**:

```bash
nanornaqtl pheno m5C \
  -b sample01_calls_sorted_map.bam \
  -o sample01_m5C \
  -t 20 \
  --motif
```

---

#### 2.3 pseudouridine (pseU)修饰位点识别

**功能**: 识别pseudouridine修饰位点。

**命令**:

```bash
nanornaqtl pheno pseU \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**参数**: 与m6A相同

**motif分类** (如启用--motif):
- **pus1**: `[ACT][AG]T`
- **pus4**: `GTTC[ATCG]A`
- **pus7**: `TGTA[AG]`

**输出文件**:

- `<prefix>_pseU_sites_result.csv`: pseU修饰位点信息

**输出列说明**: 与m5C相同,包含motif分类(pus1/pus4/pus7)

**示例**:

```bash
nanornaqtl pheno pseU \
  -b sample01_calls_sorted_map.bam \
  -o sample01_pseU \
  -t 20 \
  --motif
```

---

#### 2.4 inosine修饰位点识别

**功能**: 识别inosine (A-to-I编辑)位点。

**命令**:

```bash
nanornaqtl pheno inosine \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads> \
  [--motif] \
  [--metaPlotR]
```

**参数**: 与m6A相同

**motif** (如启用--motif): `TA[GT]`

**输出文件**:

- `<prefix>_inosine_sites_result.csv`: inosine位点信息

**示例**:

```bash
nanornaqtl pheno inosine \
  -b sample01_calls_sorted_map.bam \
  -o sample01_inosine \
  -t 20 \
  --motif
```

---

#### 2.5 polyA尾长识别

**功能**: 提取每条read的polyA尾长。

**前提条件**: BAM文件需要包含pt标签(Dorado basecall时使用`--estimate-poly-a`生成)

**命令**:

```bash
nanornaqtl pheno polyA_tail \
  -b <map_bam> \
  -o <output_prefix> \
  -t <threads>
```

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-b, --bamfile` | 输入BAM文件 | 必需 |
| `-o, --output_prefix` | 输出前缀 | 必需 |
| `-t, --threads` | 线程数(最大有效值22) | 4 |

**输出文件**:

- `<prefix>_polyAlen_result.csv`: polyA尾长信息

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `readID` | read标识符 |
| `polyA_length` | polyA尾长度 |

**示例**:

```bash
nanornaqtl pheno polyA_tail \
  -b sample01_calls_sorted_map.bam \
  -o sample01_polyA \
  -t 20
```

---

#### 2.6 内含子滞留率识别

**功能**: 计算每条read的内含子滞留率。

**命令**:

```bash
nanornaqtl pheno intron_retention \
  -g <gtf_file> \
  -b <map_bam> \
  -o <output_dir> \
  -p <output_prefix>
```

**参数说明**:

| 参数 | 说明 |
|------|------|
| `-g, --gtf` | GTF注释文件路径 |
| `-b, --bam` | 输入BAM文件 |
| `-o, --output` | 输出目录 |
| `-p, --output_prefix` | 输出文件前缀 |

**输出文件**:

- `<prefix>_intronRetention_result.csv`: 内含子滞留率信息

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `readID` | read标识符 |
| `exon_len` | 外显子长度 |
| `intron_len` | 内含子长度 |
| `len` | 总长度 |
| `IntronRetentionRate` | 内含子滞留率 |

**示例**:

```bash
nanornaqtl pheno intron_retention \
  -g gencode.v47lift37.annotation.gtf \
  -b sample01_calls_sorted_map.bam \
  -o ./results \
  -p sample01
```

---

#### 2.7 可变polyA位点 (APA)识别

**功能**: 识别每条read使用的可变polyA位点。

**命令**:

```bash
nanornaqtl pheno APA \
  -b <map_bam> \
  -o <output_prefix> \
  -a <apa_file>
```

**参数说明**:

| 参数 | 说明 |
|------|------|
| `-b, --bam` | 输入BAM文件 |
| `-o, --output_prefix` | 输出前缀 |
| `-a, --apa_file` | APA位点数据文件(BED格式) |

**APA数据来源**:

推荐从[APAdb](https://ngdc.cncb.ac.cn/databasecommons/database/id/853)下载,或使用自定义APA位点数据。

**APA文件格式要求** (BED格式,10列):

```
chr1	16442	16450	WASH7P.1:16442-16450	43	-	Intron	16443.0	16443	-
chr1	134934	134953	LOC729737.1:134934-134953	26	-	UTR3	134944.0	134944	-
```

**重要**: APA文件需要排序:

```bash
sort -k1,1 -k2,2n apa_raw.bed > apa_sorted.bed
```

**输出文件**:

- `<prefix>_APA_result.csv`: APA使用信息

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `readID` | read标识符 |
| `APA_type` | APA位点类型 |

**示例**:

```bash
nanornaqtl pheno APA \
  -b sample01_calls_sorted_map.bam \
  -o sample01_APA \
  -a apa_sorted.bed
```

---

#### 2.8 isoform识别

**功能**: 识别和定量转录本异构体。

**工具**: 使用[IsoQuant](https://github.com/ablab/IsoQuant)工具进行isoform分析。

**命令示例**:

```bash
isoquant.py \
  --reference /path/to/reference.fa \
  --genedb /path/to/annotation.gtf \
  --bam /path/to/sample_calls_sorted_map.bam \
  --data_type nanopore \
  -o output_directory \
  -t 60 \
  --complete_genedb
```

**参数说明**:

| 参数 | 说明 |
|------|------|
| `--reference` | 参考基因组FASTA文件 |
| `--genedb` | GTF注释文件 |
| `--bam` | 输入BAM文件 |
| `--data_type` | 数据类型(nanopore) |
| `-o` | 输出目录 |
| `-t` | 线程数 |
| `--complete_genedb` | 使用完整基因组数据库 |

**输出**: IsoQuant会生成多个文件,其中`OUT.transcript_model_reads.tsv.gz`用于后续QTL分析。

详细使用方法请参考[IsoQuant官方文档](https://github.com/ablab/IsoQuant)。

---

### 3. qtl - QTL分析

所有QTL分析都需要以下共同输入:

#### 变异位点文件准备 📋

**变异位点文件格式** (`--snp_info`参数):

**必需列** (列名必须完全一致,顺序不限):

| 列名 | 说明 |
|------|------|
| `CHR` | 染色体 |
| `SNP` | 变异ID(即使是Indel/SV也用此列名) |
| `POS` | 位置(1-based) |
| `A1` | 等位基因1 |
| `A2` | 等位基因2 |
| `MAF` | 最小等位基因频率(A1的频率) |

**文件示例**:

```
CHR	SNP	POS	A1	A2	MAF	NCHROBS
chr1	chr1_51898_C_A	51898	A	C	0.1082	194
chr1	chr1_51928_G_A	51928	A	G	0.02041	196
```

**生成方法**:

使用plink生成:
```bash
plink --freq --bfile your_data --out snp_info
```

或自行构建,只需保证必需列存在且列名一致。可以包含其他列(如`NCHROBS`),不影响分析。

**重要说明**:

- **支持的变异类型**: SNV, Indel, SV
  - ⚠️ **注意**: 代码中对Indel和SV没有专门的算法设置,在分析时将Indel和SV当作SNV处理
- **MAF阈值**: 建议对变异进行筛选,处理MAF > 0.05的变异
  - MAF < 0.05的变异统计power较低,可能结果不佳
  - MAF > 0.05是标准的QTL过滤阈值
  - 非强制要求,用户可根据研究需求调整
- 建议对文件按CHR和POS排序

---

#### 3.1 m6A QTL分析

**功能**: 识别与m6A修饰相关的遗传变异。

**命令**:

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

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-b, --bam` | 输入BAM文件 | 必需 |
| `--snp_info` | 变异位点文件 | 必需 |
| `-o, --output_prefix` | 输出前缀 | 必需 |
| `-csv, --modification` | m6A位点结果文件(`*_m6A_sites_result.csv`) | 必需 |
| `-pkl, --read_mod_dict` | m6A reads字典文件(`*_m6A_reads_final.pkl`) | 必需 |
| `--geno_size` | 基因组大小文件(如hg19.chrom.sizes) | 必需 |
| `-m, --modification_type` | 修饰类型(m6A) | 必需 |
| `-q, --min_qscore` | 最小碱基质量 | 10 |
| `-c, --min_coverage` | 最小总覆盖度 | 8 |
| `--mcmc_samples` | MCMC采样数 | 1000 |
| `--threads` | 线程数 | 4 |
| `--keep_tmp` | 保留临时文件 | False |

**输出文件**:

- `<prefix>_m6A_QTLs_result.csv`: m6A QTL结果

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `chrom` | 染色体 |
| `SNP` | 变异ID |
| `snp_pos_1base` | SNP位置 |
| `A1`, `A2` | 等位基因 |
| `MAF` | 最小等位基因频率 |
| `mod_pos_1base` | 修饰位点位置 |
| `BF` | Bayes Factor(贝叶斯因子) |
| `posterior_prob` | 后验概率 |
| `fisher_pvalue` | Fisher精确检验p值 |

**统计方法说明**:

1. **贝叶斯方法**(推荐):
   - **BF (Bayes Factor)**: 显著性指标,BF > 3一般认为有显著关联
   - **后验概率**: 零假设的后验概率
   - **优势**: 对Nanopore数据覆盖度不均匀(最大/最小可差100倍)更稳健

2. **频率派方法**(参考):
   - **Fisher精确检验**: 更激进,假阳性率较高

**注意**: 输出包含所有统计结果(显著和不显著),用户可根据自己的阈值(如BF > 3)筛选。

**示例**:

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

#### 3.2 m5C QTL分析

**功能**: 识别与m5C修饰相关的遗传变异。

**命令**:

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

**参数**: 与m6A QTL相同,将`-m`参数改为`m5C`

**输出文件**:

- `<prefix>_m5C_QTLs_result.csv`: m5C QTL结果

**输出列**: 与m6A QTL相同

**示例**:

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

#### 3.3 pseudouridine (pseU) QTL分析

**功能**: 识别与pseudouridine修饰相关的遗传变异。

**命令**:

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

**参数**: 与m6A QTL相同,将`-m`参数改为`pseU`

**输出文件**:

- `<prefix>_pseU_QTLs_result.csv`: pseU QTL结果

**示例**:

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

#### 3.4 inosine QTL分析

**功能**: 识别与inosine修饰相关的遗传变异。

**命令**:

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

**参数**: 与m6A QTL相同,将`-m`参数改为`inosine`

**输出文件**:

- `<prefix>_inosine_QTLs_result.csv`: inosine QTL结果

**示例**:

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

#### 3.5 APA QTL分析

**功能**: 识别影响可变polyA位点使用的遗传变异。

**命令**:

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

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-b, --bam` | 输入BAM文件 | 必需 |
| `--snp_info` | 变异位点文件 | 必需 |
| `-o, --output_prefix` | 输出前缀 | 必需 |
| `-f, --read_overlap_file` | APA结果文件(`*_APA_result.csv`) | 必需 |
| `--geno_size` | 基因组大小文件 | 必需 |
| `-m, --molecular_type` | 分子表型类型(APA) | 必需 |
| `-q, --min_qscore` | 最小碱基质量 | 10 |
| `--min_coverage` | 最小覆盖度 | 8 |
| `--mcmc_samples` | MCMC采样数 | 1000 |
| `-t, --threads` | 线程数 | 4 |
| `--keep_tmp` | 保留临时文件 | False |

**输出文件**:

- `<prefix>_APA_QTLs_result.csv`: APA QTL结果

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `chrom` | 染色体 |
| `SNP` | 变异ID |
| `snp_pos_1base` | SNP位置 |
| `A1`, `A2` | 等位基因 |
| `MAF` | 最小等位基因频率 |
| `BF` | Bayes Factor |
| `posterior_prob` | 后验概率 |
| `chi2_pvalue` | Chi-square test(卡方检验) p值 |
| `TVD` | Total Variation Distance(总变差距离,范围0-1) |
| `dominant_shift` | 主要APA类型使用变化方向和大小 |

**统计方法说明**:

1. **贝叶斯方法**(推荐): BF和后验概率
2. **频率派方法**(参考): Chi-square test (卡方检验)

**效应量指标说明**:

- **TVD** (Total Variation Distance): 衡量A1和A2等位基因对应的APA使用分布差异程度
  - 0: 两个分布完全相同
  - 1: 两个分布完全不同
  
- **dominant_shift**: 描述每种APA类型在A1和A2之间的使用比例变化
  - 格式示例: `type1:A1↓(-0.56)` 表示携带A1等位基因的reads较少使用该APA类型(相对减少56%)

**示例**:

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

#### 3.6 isoform QTL分析

**功能**: 识别影响转录本异构体使用的遗传变异。

**命令**:

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

**参数说明**:

| 参数 | 说明 |
|------|------|
| `-f, --read_overlap_file` | IsoQuant输出文件(`OUT.transcript_model_reads.tsv.gz`) |
| 其他参数 | 与APA QTL相同 |

**输出文件**:

- `<prefix>_isoform_QTLs_result.csv`: isoform QTL结果

**输出列**: 与APA QTL相同,其中:
- `chi2_pvalue`: Chi-square test (卡方检验) p值
- `TVD`和`dominant_shift`用于描述isoform使用模式差异

**示例**:

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

#### 3.7 polyA尾长QTL分析

**功能**: 识别影响polyA尾长的遗传变异。

**命令**:

```bash
nanornaqtl qtl polyA_tail \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <polyA_csv> \
  --geno_size <genome_size_file> \
  --threads <threads>
```

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-b, --bam` | 输入BAM文件 | 必需 |
| `--snp_info` | 变异位点文件 | 必需 |
| `-o, --output_prefix` | 输出前缀 | 必需 |
| `-csv, --polya_csv` | polyA结果文件(`*_polyAlen_result.csv`) | 必需 |
| `--geno_size` | 基因组大小文件 | 必需 |
| `-q, --min_qscore` | 最小碱基质量 | 10 |
| `--min_coverage` | 最小覆盖度 | 8 |
| `--mcmc_samples` | MCMC采样数 | 2000 |
| `--threads` | 线程数 | 4 |
| `--keep_tmp` | 保留临时文件 | False |

**输出文件**:

- `<prefix>_polyA_tail_length_QTLs_result.csv`: polyA尾长QTL结果

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `chrom` | 染色体 |
| `SNP` | 变异ID |
| `snp_pos_1base` | SNP位置 |
| `A1`, `A2` | 等位基因 |
| `MAF` | 最小等位基因频率 |
| `A1_len` | A1等位基因的polyA尾长列表 |
| `A2_len` | A2等位基因的polyA尾长列表 |
| `beta` | 效应量(log ratio of mean lengths) |
| `SE` | 标准误 |
| `KS_stat` | Kolmogorov-Smirnov统计量 |
| `bayes_factor` | Bayes Factor |
| `posterior_prob` | 后验概率 |
| `p_welch` | Welch's t-test p值 |
| `p_mw` | Mann-Whitney U test p值 |
| `p_ks` | Kolmogorov-Smirnov test p值 |

**频率派方法说明**(三种p值):

- **p_welch**: Welch's t-test,适用于方差不齐的情况
- **p_mw**: Mann-Whitney U test,非参数检验,对异常值稳健
- **p_ks**: Kolmogorov-Smirnov test,检验分布差异

**示例**:

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

#### 3.8 内含子滞留率QTL分析

**功能**: 识别影响内含子滞留率的遗传变异。

**命令**:

```bash
nanornaqtl qtl intron_retention \
  -b <map_bam> \
  --snp_info <snp_file> \
  -o <output_prefix> \
  -csv <ir_csv> \
  --geno_size <genome_size_file> \
  --threads <threads>
```

**参数说明**:

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-csv, --ir_csv` | 内含子滞留率结果文件(`*_intronRetention_result.csv`) | 必需 |
| 其他参数 | 与polyA尾长QTL相同 |

**输出文件**:

- `<prefix>_intron_retention_QTLs_result.csv`: 内含子滞留率QTL结果

**输出列说明**:

| 列名 | 说明 |
|------|------|
| `chrom` | 染色体 |
| `SNP` | 变异ID |
| `snp_pos_1base` | SNP位置 |
| `A1`, `A2` | 等位基因 |
| `MAF` | 最小等位基因频率 |
| `A1_IR` | A1等位基因的内含子滞留率列表 |
| `A2_IR` | A2等位基因的内含子滞留率列表 |
| `beta` | 效应量(logit差异) |
| `SE` | 标准误 |
| `KS_stat` | Kolmogorov-Smirnov统计量 |
| `bayes_factor` | Bayes Factor |
| `posterior_prob` | 后验概率 |
| `p_welch` | Welch's t-test p值 |
| `p_mw` | Mann-Whitney U test p值 |
| `p_ks` | Kolmogorov-Smirnov test p值 |

**注意**: 自动过滤低质量数据(如90%以上reads的IR=0)

**示例**:

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

## 输出文件说明

### 文件命名规则

- **prep模块**: `<prefix>_calls_sorted_map[0/16].bam`, `<prefix>_calls_sorted_map.fastq`
- **pheno模块**: `<prefix>_<phenotype>_result.csv`
- **qtl模块**: `<prefix>_<phenotype>_QTLs_result.csv`

### 通用输出格式

所有CSV输出文件均为标准的制表符或逗号分隔文件,可用Excel、R、Python等工具打开和分析。

### 关键统计指标解读 📈

#### 贝叶斯方法

- **BF (Bayes Factor)**:
  - BF > 3: 有显著证据支持关联
  - BF > 10: 强证据
  - BF > 30: 非常强证据
  
- **后验概率**: 零假设为真的后验概率,越小表示越可能存在关联

#### 频率派方法

- **p值**: 传统显著性检验p值,通常p < 0.05认为显著
- **注意**: 频率派方法在覆盖度不均匀时可能产生较多假阳性,建议优先使用贝叶斯方法

---

## 引用

如果您在研究中使用了nanornaqtl,请引用我们的GitHub仓库:

```
nanornaqtl: A comprehensive toolkit for molecular phenotyping and QTL analysis using Nanopore direct RNA sequencing
https://github.com/xinranxu0930/nanornaqtl
```

---

## 联系方式 📧

- **问题反馈**: [GitHub Issues](https://github.com/xinranxu0930/nanornaqtl/issues)
- **邮箱**: xinranxu0930@gmail.com
- **GitHub**: [https://github.com/xinranxu0930/nanornaqtl](https://github.com/xinranxu0930/nanornaqtl)

---

## 许可证

本项目采用MIT许可证。详见[LICENSE](LICENSE)文件。

---

## 更新日志

### v1.0.0 (2026-01-18)

- 初始版本发布
- 支持7种分子表型识别
- 支持8种QTL分析
- 实现贝叶斯统计方法
- 支持多线程并行处理
