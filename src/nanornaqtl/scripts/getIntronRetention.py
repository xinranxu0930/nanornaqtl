#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import sys
from subprocess import call


def process_gtf(gtf_file, output_dir):
    """
    处理GTF文件，生成基因和外显子的BED文件
    """
    print("开始处理GTF文件...")

    # 读取GTF文件
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        names=["chr", "type", "start", "end", "strand", "info"],
    )

    # 处理基因信息
    print("处理基因信息...")
    gene_gtf = gtf[gtf["type"] == "gene"].copy()
    gene_gtf.loc[:, "geneID"] = gene_gtf["info"].str.extract(r'gene_id "(.*?)";')
    gene_gtf.loc[:, "start"] = gene_gtf["start"] - 1
    gene_gtf.loc[:, "score"] = 0
    gene_gtf = gene_gtf[["chr", "start", "end", "geneID", "score", "strand"]]
    gene_gtf = gene_gtf[~gene_gtf["chr"].isin(["chrX", "chrY", "chrM"])]

    # 保存并合并基因
    gene_bed = os.path.join(output_dir, "gene.bed")
    gene_gtf.to_csv(gene_bed, sep="\t", index=False, header=False)
    call(
        f"bedtools merge -i {output_dir}/gene.bed -s -c 4,5,6 -o collapse > {output_dir}/gene_merge.bed",
        shell=True,
    )
    gene_gtf = pd.read_csv(
        f"{output_dir}/gene_merge.bed",
        sep="\t",
        header=None,
        names=["chr", "start", "end", "geneID", "score", "strand"],
    )
    gene_gtf.loc[:, "score"] = 0
    gene_gtf.loc[:, "strand"] = (
        gene_gtf["strand"]
        .astype(str)
        .apply(lambda x: x.split(",")[0] if "," in x else x)
    )
    gene_gtf.to_csv(f"{output_dir}/gene_merge.bed", sep="\t", index=False, header=False)

    # 处理外显子信息
    print("处理外显子信息...")
    exon_gtf = gtf[gtf["type"] == "exon"].copy()
    exon_gtf.loc[:, "geneID"] = exon_gtf["info"].str.extract(r'gene_id "(.*?)";')
    exon_gtf.loc[:, "start"] = exon_gtf["start"] - 1
    exon_gtf.loc[:, "score"] = 0
    exon_gtf = exon_gtf[["chr", "start", "end", "geneID", "score", "strand"]]
    exon_gtf = exon_gtf[~exon_gtf["chr"].isin(["chrX", "chrY", "chrM"])]
    exon_gtf.to_csv(f"{output_dir}/exon.bed", sep="\t", index=False, header=False)
    call(
        f"sort -k 1,1 -k2,2n {output_dir}/exon.bed > {output_dir}/exon_sorted.bed",
        shell=True,
    )
    call(
        f"bedtools merge -i {output_dir}/exon_sorted.bed -s -c 4,5,6 -o collapse > {output_dir}/exon_merge.bed",
        shell=True,
    )
    exon_merged_df = pd.read_csv(
        f"{output_dir}/exon_merge.bed",
        sep="\t",
        header=None,
        names=["chr", "start", "end", "geneID", "score", "strand"],
    )
    exon_merged_df.loc[:, "score"] = 0
    exon_merged_df.loc[:, "strand"] = (
        exon_merged_df["strand"]
        .astype(str)
        .apply(lambda x: x.split(",")[0] if "," in x else x)
    )
    exon_merged_df.to_csv(
        f"{output_dir}/exon_merge.bed", sep="\t", index=False, header=False
    )
    call(f"rm {output_dir}/exon_sorted.bed", shell=True)

    # 生成内含子文件
    print("生成内含子信息...")
    call(
        f"bash -c 'bedtools subtract -a {output_dir}/gene_merge.bed -b {output_dir}/exon_merge.bed -s > {output_dir}/intron.bed'",
        shell=True,
    )
    print("GTF处理完成")


def identify_valid_reads(bam_file, output_dir):
    print("开始识别有效reads...")

    # 生成必要的基因/外显子文件
    call(
        f"bedtools intersect -a {output_dir}/gene_merge.bed -b {output_dir}/exon_merge.bed -s -wa -wb > {output_dir}/fully_overlapped_genes.bed",
        shell=True,
    )
    df = pd.read_csv(
        f"{output_dir}/fully_overlapped_genes.bed",
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3, 4, 5, 7, 8],
        names=["chr", "start", "end", "gene", "score", "strand", "s_t", "e_t"],
    )
    df_no_intron = df[(df["start"] == df["s_t"]) & (df["end"] == df["e_t"])].copy()
    df_no_intron = df_no_intron.drop(["s_t", "e_t"], axis=1)
    df_no_intron.to_csv(
        f"{output_dir}/genes_no_intron.bed", sep="\t", index=False, header=False
    )
    call(f"rm {output_dir}/fully_overlapped_genes.bed", shell=True)
    
    # 使用管道和命令行工具处理readID
    # 1. 获取所有reads的ID，并排序去重
    print("获取所有read IDs...")
    call(
        f"bedtools bamtobed -i {bam_file} | cut -f 4 | sort -u > {output_dir}/all_reads.txt",
        shell=True,
    )

    # 2. 获取没有比对到基因上的reads的ID
    print("获取未比对到基因的read IDs...")
    call(
        f"bash -c 'bedtools subtract -a <(bedtools bamtobed -i {bam_file}) -b {output_dir}/gene_merge.bed -A -s | cut -f 4 | sort -u > {output_dir}/read_no_gene.txt'",
        shell=True,
    )

    # 3. 获取比对到无内含子基因上的reads的ID
    print("获取比对到无内含子基因的read IDs...")
    call(
        f"bedtools intersect -a <(bedtools bamtobed -i {bam_file}) -b {output_dir}/genes_no_intron.bed -s -wa | cut -f 4 | sort -u > {output_dir}/read_gene_no_intron.txt",
        shell=True,
    )
    
    # 4. 使用comm命令高效地找出差集
    # 先合并两个要排除的ID集合
    print("合并无效read IDs...")
    call(
        f"cat {output_dir}/read_no_gene.txt {output_dir}/read_gene_no_intron.txt | sort -u > {output_dir}/invalid_reads.txt",
        shell=True
    )
    
    # comm -23 file1 file2 会打印出只在file1中存在的行
    print("计算有效read IDs...")
    call(
        f"comm -23 {output_dir}/all_reads.txt {output_dir}/invalid_reads.txt > {output_dir}/valid_readID.txt",
        shell=True
    )
    
    # 清理中间的ID文件
    call(f"rm {output_dir}/all_reads.txt {output_dir}/read_no_gene.txt {output_dir}/read_gene_no_intron.txt {output_dir}/invalid_reads.txt {output_dir}/genes_no_intron.bed", shell=True)
    print("有效reads识别完成!")


def get_valid_bam(bam_file, output_dir, output_prefix):
    call(f"samtools view -N {output_dir}/valid_readID.txt {bam_file} -b > {output_dir}/{output_prefix}_calls_sorted_mod_map_valid.bam", shell=True)
    call(f"samtools index {output_dir}/{output_prefix}_calls_sorted_mod_map_valid.bam", shell=True)
    call(f"bedtools bamtobed -splitD -i {output_dir}/{output_prefix}_calls_sorted_mod_map_valid.bam > {output_dir}/{output_prefix}_calls_sort_map_valid_bam2bed.bed", shell=True)
    call(f"sort -k 1,1 -k2,2n {output_dir}/{output_prefix}_calls_sort_map_valid_bam2bed.bed > {output_dir}/valid_bam2bed_sort.bed", shell=True)
    call(f"bedtools intersect -s -a {output_dir}/valid_bam2bed_sort.bed -b {output_dir}/intron.bed > {output_dir}/{output_prefix}_intron.bed", shell=True)
    call(f"bedtools intersect -s -a {output_dir}/valid_bam2bed_sort.bed -b {output_dir}/exon_merge.bed > {output_dir}/{output_prefix}_exon.bed", shell=True)


def calculate_intron_retention(output_dir, output_prefix):
    """
    优化版本的内含子滞留率计算函数
    """
    print("开始计算内含子滞留率...")
    
    # 首先获取所有唯一的readID
    print("读取所有readID...")
    all_reads_df = pd.read_csv(
        f"{output_dir}/valid_bam2bed_sort.bed",
        sep="\t",
        header=None,
        usecols=[3],
        names=["readID"],
    )
    unique_reads = set(all_reads_df["readID"])
    print(f"总共有 {len(unique_reads)} 个unique reads")
    
    # 使用字典来存储结果，避免多次merge
    read_stats = {read_id: {"exon_len": 0, "intron_len": 0} for read_id in unique_reads}
    
    # 处理exon重叠
    print("处理exon重叠...")
    try:
        exon_df = pd.read_csv(
            f"{output_dir}/{output_prefix}_exon.bed",
            sep="\t",
            header=None,
            usecols=[1, 2, 3],
            names=["s", "e", "readID"],
            chunksize=10000  # 分块读取以减少内存使用
        )
        
        for chunk in exon_df:
            chunk["exon_overlap"] = chunk["e"] - chunk["s"]
            exon_grouped = chunk.groupby("readID")["exon_overlap"].sum()
            
            for read_id, overlap in exon_grouped.items():
                if read_id in read_stats:
                    read_stats[read_id]["exon_len"] += overlap
                    
    except FileNotFoundError:
        print("警告: exon.bed 文件不存在，所有reads的exon长度将为0")
    except pd.errors.EmptyDataError:
        print("警告: exon.bed 文件为空")
    
    # 处理intron重叠
    print("处理intron重叠...")
    try:
        intron_df = pd.read_csv(
            f"{output_dir}/{output_prefix}_intron.bed",
            sep="\t",
            header=None,
            usecols=[1, 2, 3],
            names=["s", "e", "readID"],
            chunksize=10000  # 分块读取
        )
        
        for chunk in intron_df:
            chunk["intron_overlap"] = chunk["e"] - chunk["s"]
            intron_grouped = chunk.groupby("readID")["intron_overlap"].sum()
            
            for read_id, overlap in intron_grouped.items():
                if read_id in read_stats:
                    read_stats[read_id]["intron_len"] += overlap
                    
    except FileNotFoundError:
        print("警告: intron.bed 文件不存在，所有reads的intron长度将为0")
    except pd.errors.EmptyDataError:
        print("警告: intron.bed 文件为空")
    
    # 转换为DataFrame并计算最终结果
    print("计算最终结果...")
    results_data = []
    for read_id, stats in read_stats.items():
        exon_len = stats["exon_len"]
        intron_len = stats["intron_len"]
        
        # 只保留有exon覆盖的reads
        if exon_len > 0:
            total_len = exon_len + intron_len
            stability = intron_len / total_len if total_len > 0 else 0
            
            results_data.append({
                "readID": read_id,
                "exon_len": int(exon_len),
                "intron_len": int(intron_len),
                "len": int(total_len),
                "IntronRetentionRate": stability
            })
    
    # 创建最终DataFrame
    df = pd.DataFrame(results_data)
    print(f"最终结果包含 {len(df)} 个reads")
    
    # 保存结果
    df.to_csv(f"{output_dir}/{output_prefix}_intronRetention_result.csv", index=False)
    
    # 清理临时文件
    cleanup_files = [
        f"{output_dir}/valid_bam2bed_sort.bed",
        f"{output_dir}/{output_prefix}_calls_sorted_mod_map_valid.bam",
        f"{output_dir}/{output_prefix}_calls_sorted_mod_map_valid.bam.bai",
        f"{output_dir}/{output_prefix}_calls_sort_map_valid_bam2bed.bed",
        f"{output_dir}/{output_prefix}_intron.bed",
        f"{output_dir}/{output_prefix}_exon.bed"
    ]
    
    for file_path in cleanup_files:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
        except OSError:
            print(f"警告: 无法删除文件 {file_path}")
    
    # 清理其他临时文件
    other_cleanup_files = [
        f"{output_dir}/valid_readID.txt",
        f"{output_dir}/exon_merge.bed", 
        f"{output_dir}/gene_merge.bed"
    ]
    
    for file_path in other_cleanup_files:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
        except OSError:
            print(f"警告: 无法删除文件 {file_path}")
    
    print("内含子滞留率计算完成!")


def main():
    parser = argparse.ArgumentParser(
        description="内含子滞留率分析流程 - 整合GTF处理、BAM过滤和内含子滞留率计算",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python getIntronRetention.py \\
    --gtf gencode.annotation.gtf \\
    --bam calls_sorted_map.bam \\
    --output ./results \\
    --output_prefix sample1
        """,
    )

    parser.add_argument("-g", "--gtf", required=True, help="GTF注释文件路径")
    parser.add_argument("-b", "--bam", required=True, help="输入BAM文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出目录路径")
    parser.add_argument("-p", "--output_prefix", required=True, help="输出文件前缀")
    args = parser.parse_args()

    # 检查输入文件
    if not os.path.exists(args.gtf):
        print(f"GTF文件不存在: {args.gtf}")
        sys.exit(1)
    if not os.path.exists(args.bam):
        print(f"BAM文件不存在: {args.bam}")
        sys.exit(1)

    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)

    try:
        # 步骤1: 处理GTF文件
        print("1.处理GTF文件...", flush=True)
        process_gtf(args.gtf, args.output)
        print("GTF文件处理完成!")

        # 步骤2: 识别有效reads
        print("2.识别有效reads...", flush=True)
        identify_valid_reads(args.bam, args.output)
        print("有效reads识别完成!")

        # 步骤3: 筛选出只含有有效reads的bam文件
        print("3.筛选出只包含有效reads的bam文件...", flush=True)
        get_valid_bam(args.bam, args.output, args.output_prefix)
        print("筛选bam文件完成!")

        # 步骤4: 计算内含子滞留率
        print("4.计算内含子滞留率...", flush=True)
        calculate_intron_retention(args.output, args.output_prefix)
        print("内含子滞留率计算完成!")

    except Exception as e:
        print(f"分析过程中出现错误: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
