#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
APA (Alternative Polyadenylation) Analysis Script - Parallel Version
从BAM文件中识别PAS位点，并聚类为APA类型
支持按染色体多线程并行处理
"""

import pandas as pd
import pysam
import argparse
import subprocess
import os
import tempfile
from datetime import datetime
from multiprocessing import Pool, cpu_count

def log_print(message):
    """打印带时间戳的日志信息"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def get_strand_bam_paths(bam_path):
    """
    根据用户输入的BAM路径，推断正链和负链BAM文件路径
    """
    if bam_path.endswith('.bam'):
        base_path = bam_path[:-4]
        bam0 = f"{base_path}0.bam"
        bam16 = f"{base_path}16.bam"
    else:
        bam0 = f"{bam_path}0.bam"
        bam16 = f"{bam_path}16.bam"
    
    return bam0, bam16

def get_chromosomes_from_bam(bamfile):
    """
    从BAM文件获取所有染色体列表
    """
    samfile = pysam.AlignmentFile(bamfile, "rb")
    chroms = list(samfile.references)
    samfile.close()
    return chroms

def get_pas_from_bam_by_chrom(args_tuple):
    """
    从BAM文件中提取指定染色体的PAS位点（用于并行处理）
    
    Parameters:
    -----------
    args_tuple : tuple
        (bamfile, chrom, strand)
    
    Returns:
    --------
    dict : {read_id: (chrom, pos)}
    """
    bamfile, chrom, strand = args_tuple
    res_dict = dict()
    
    try:
        samfile = pysam.AlignmentFile(bamfile, "rb")
        
        for read in samfile.fetch(chrom):
            if read.is_unmapped or read.cigar is None:
                continue
            
            if strand == '+':
                if read.cigar[-1][0] == 4:  # soft-clip
                    polyA_cigar_len = read.cigar[-1][1]
                    polyA_cigar_seq = read.query_sequence[-polyA_cigar_len:]
                    end_pos = read.reference_end
                    base = 'A'
                else:
                    continue
            else:  # strand == '-'
                if read.cigar[0][0] == 4:  # soft-clip
                    polyA_cigar_len = read.cigar[0][1]
                    polyA_cigar_seq = read.query_sequence[:polyA_cigar_len]
                    end_pos = read.reference_start
                    base = 'T'
                else:
                    continue
            
            # 判断是否为poly(A)尾
            if polyA_cigar_len >= 20:
                ratio_front = polyA_cigar_seq[:20].count(base) / 20
                ratio_end = polyA_cigar_seq[-20:].count(base) / 20
                if strand == '+':
                    condition = (ratio_front > 0.7 and ratio_end > 0.8)
                else:
                    condition = (ratio_front > 0.8 and ratio_end > 0.7)
            else:
                condition = polyA_cigar_seq.count(base) / polyA_cigar_len > 0.8 if polyA_cigar_len > 0 else False
            
            if condition:
                res_dict[read.query_name] = (chrom, end_pos)
        
        samfile.close()
    except Exception as e:
        print(f"Warning: 处理 {chrom} ({strand}) 时出错: {e}")
    
    return res_dict

def get_pas_parallel(bam0, bam16, threads):
    """
    并行处理正负链BAM文件，按染色体分配任务
    
    Parameters:
    -----------
    bam0 : str
        正链BAM文件路径
    bam16 : str
        负链BAM文件路径
    threads : int
        线程数
    
    Returns:
    --------
    tuple : (pas_plus_dict, pas_minus_dict)
    """
    # 获取染色体列表
    chroms_plus = get_chromosomes_from_bam(bam0)
    chroms_minus = get_chromosomes_from_bam(bam16)
    
    log_print(f"  正链染色体数: {len(chroms_plus)}")
    log_print(f"  负链染色体数: {len(chroms_minus)}")
    
    # 构建任务列表：(bamfile, chrom, strand)
    tasks = []
    for chrom in chroms_plus:
        tasks.append((bam0, chrom, '+'))
    for chrom in chroms_minus:
        tasks.append((bam16, chrom, '-'))
    
    log_print(f"  总任务数: {len(tasks)}, 使用线程数: {threads}")
    
    # 并行处理
    with Pool(processes=threads) as pool:
        results = pool.map(get_pas_from_bam_by_chrom, tasks)
    
    # 合并结果
    pas_plus = {}
    pas_minus = {}
    
    for i, task in enumerate(tasks):
        _, chrom, strand = task
        if strand == '+':
            pas_plus.update(results[i])
        else:
            pas_minus.update(results[i])
    
    return pas_plus, pas_minus

def filter_internal_priming(pas_dict, fasta_file, strand):
    """
    过滤内部priming假阳性
    """
    if len(pas_dict) == 0:
        return pas_dict
    
    fasta = pysam.FastaFile(fasta_file)
    filtered_dict = dict()
    base = 'A' if strand == '+' else 'T'
    
    # 按染色体分组处理
    chrom_dict = {}
    for read_id, (chrom, pos) in pas_dict.items():
        if chrom not in chrom_dict:
            chrom_dict[chrom] = []
        chrom_dict[chrom].append((read_id, pos))
    
    for chrom, read_list in chrom_dict.items():
        try:
            chrom_seq = fasta.fetch(chrom).upper()
            chrom_len = len(chrom_seq)
            
            for read_id, pos in read_list:
                start = max(0, pos - 10)
                end = min(chrom_len, pos + 11)
                local_seq = chrom_seq[start:end]
                
                if local_seq.count(base) <= 10:
                    filtered_dict[read_id] = (chrom, pos)
        except Exception as e:
            log_print(f"Warning: 无法获取染色体 {chrom} 的序列: {e}")
            for read_id, pos in read_list:
                filtered_dict[read_id] = (chrom, pos)
    
    fasta.close()
    return filtered_dict

def create_bed_df(pas_dict, strand):
    """
    将PAS位点转换为BED格式的DataFrame
    """
    if len(pas_dict) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'readID', 'score', 'strand'])
    
    rows = []
    for read_id, (chrom, pos) in pas_dict.items():
        start = pos
        end = pos + 1
        rows.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'readID': read_id,
            'score': 0,
            'strand': strand
        })
    
    df = pd.DataFrame(rows)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    return df

def run_bedtools_merge(bed_file, output_file, distance):
    """
    运行bedtools merge合并相近的PAS位点
    """
    sorted_bed = bed_file + ".sorted"
    cmd_sort = f"bedtools sort -i {bed_file} > {sorted_bed}"
    subprocess.call(cmd_sort, shell=True)
    
    cmd_merge = f"bedtools merge -i {sorted_bed} -s -d {distance} -c 4,6 -o collapse,distinct > {output_file}"
    subprocess.call(cmd_merge, shell=True)
    
    os.remove(sorted_bed)

def map_reads_to_apa(merged_file, pas_df):
    """
    将每条read映射到对应的APA区间
    """
    merged_df = pd.read_csv(merged_file, sep='\t', header=None,
                            names=['chrom', 'start', 'end', 'readIDs', 'strand'])
    
    # 构建readID到PAS_site的映射字典，避免重复查询DataFrame
    read_to_pas = {}
    for _, row in pas_df.iterrows():
        read_to_pas[row['readID']] = f"{row['chrom']}_{row['start']}_{row['strand']}"
    
    results = []
    for _, row in merged_df.iterrows():
        apa_type = f"{row['chrom']}:{row['start']}:{row['end']}:{row['strand']}"
        read_ids = row['readIDs'].split(',')
        
        for read_id in read_ids:
            if read_id in read_to_pas:
                results.append({
                    'readID': read_id,
                    'APA_type': apa_type,
                    'PAS_site': read_to_pas[read_id]
                })
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="APA Analysis: 从BAM文件识别PAS位点并聚类为APA类型 (并行版本)")
    parser.add_argument("-b", "--bam", type=str, required=True, 
                        help="Mapped BAM file (会自动推断*0.bam和*16.bam)")
    parser.add_argument("-f", "--fasta", type=str, required=True,
                        help="Reference genome fasta file")
    parser.add_argument("-o", "--output_prefix", type=str, required=True,
                        help="Output prefix")
    parser.add_argument("-d", "--distance", type=int, default=35,
                        help="Merge window size for bedtools merge (default: 35)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Number of threads for parallel processing (default: 4, max: 44)")
    args = parser.parse_args()
    
    # 限制最大线程数为44
    threads = min(args.threads, 44)
    
    # ========== Step 1: 推断正负链BAM文件路径 ==========
    log_print("Step 1: 推断正负链BAM文件路径...")
    bam0, bam16 = get_strand_bam_paths(args.bam)
    
    if not os.path.exists(bam0):
        raise FileNotFoundError(f"正链BAM文件不存在: {bam0}")
    if not os.path.exists(bam16):
        raise FileNotFoundError(f"负链BAM文件不存在: {bam16}")
    
    log_print(f"  正链BAM: {bam0}")
    log_print(f"  负链BAM: {bam16}")
    
    # ========== Step 2: 并行提取PAS位点 ==========
    log_print(f"Step 2: 并行提取PAS位点 (threads={threads})...")
    pas_plus, pas_minus = get_pas_parallel(bam0, bam16, threads)
    log_print(f"  正链原始PAS位点数: {len(pas_plus)}")
    log_print(f"  负链原始PAS位点数: {len(pas_minus)}")
    
    # ========== Step 3: 过滤内部priming假阳性 ==========
    log_print("Step 3: 过滤内部priming假阳性...")
    pas_plus_filtered = filter_internal_priming(pas_plus, args.fasta, '+')
    pas_minus_filtered = filter_internal_priming(pas_minus, args.fasta, '-')
    log_print(f"  正链过滤后PAS位点数: {len(pas_plus_filtered)}")
    log_print(f"  负链过滤后PAS位点数: {len(pas_minus_filtered)}")
    
    # ========== Step 4: 转换为BED格式 ==========
    log_print("Step 4: 转换为BED格式...")
    bed_plus = create_bed_df(pas_plus_filtered, '+')
    bed_minus = create_bed_df(pas_minus_filtered, '-')
    
    bed_all = pd.concat([bed_plus, bed_minus], ignore_index=True)
    log_print(f"  总PAS位点数: {len(bed_all)}")
    
    if len(bed_all) == 0:
        log_print("Warning: 未找到任何PAS位点，程序退出")
        return
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
        tmp_bed_path = tmp_bed.name
        bed_all[['chrom', 'start', 'end', 'readID', 'score', 'strand']].to_csv(
            tmp_bed_path, sep='\t', index=False, header=False)
    
    # ========== Step 5: 运行bedtools merge ==========
    log_print(f"Step 5: 运行bedtools merge (distance={args.distance}bp)...")
    tmp_merged_path = tmp_bed_path + ".merged"
    run_bedtools_merge(tmp_bed_path, tmp_merged_path, args.distance)
    
    merged_count = sum(1 for _ in open(tmp_merged_path))
    log_print(f"  合并后APA区间数: {merged_count}")
    
    # ========== Step 6: 映射reads到APA区间 ==========
    log_print("Step 6: 映射reads到APA区间...")
    result_df = map_reads_to_apa(tmp_merged_path, bed_all)
    log_print(f"  最终结果行数: {len(result_df)}")
    
    # ========== Step 7: 保存结果 ==========
    output_file = f"{args.output_prefix}_APA_result.csv"
    log_print(f"Step 7: 保存结果到 {output_file}...")
    result_df.to_csv(output_file, index=False)
    
    # 清理临时文件
    os.remove(tmp_bed_path)
    os.remove(tmp_merged_path)
    
    log_print("完成!")
    log_print(f"  输出文件: {output_file}")
    log_print(f"  总reads数: {len(result_df)}")
    log_print(f"  总APA类型数: {result_df['APA_type'].nunique()}")

if __name__ == "__main__":
    main()