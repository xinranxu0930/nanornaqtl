import pysam
import pandas as pd
import re
import argparse
import sys
import os
from typing import List, Tuple, Optional, Set
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm



# ============ Motif提取函数  ============

def extract_inosine_motif(seq, idx, strand):
    """提取inosine motif"""
    if strand == "+":
        return seq[idx - 1 : idx + 2]
    else:
        return get_reverse_complementary_sequence(seq[idx - 1 : idx + 2])


def extract_m5C_motif_3bp(seq, idx, strand):
    """提取m5C 3bp motif (CHH/CHG)"""
    if strand == "+":
        return seq[idx : idx + 3]
    else:
        return get_reverse_complementary_sequence(seq[idx - 2 : idx + 1])


def extract_m5C_motif_2bp(seq, idx, strand):
    """提取m5C 2bp motif (CG)"""
    if strand == "+":
        return seq[idx : idx + 2]
    else:
        return get_reverse_complementary_sequence(seq[idx - 1 : idx + 1])


def extract_pseU_pus1_motif(seq, idx, strand):
    """提取pseU pus1 motif"""
    if strand == "+":
        return seq[idx - 2 : idx + 1]
    else:
        return get_reverse_complementary_sequence(seq[idx : idx + 3])


def extract_pseU_pus4_motif(seq, idx, strand):
    """提取pseU pus4 motif"""
    if strand == "+":
        return seq[idx - 2 : idx + 4]
    else:
        return get_reverse_complementary_sequence(seq[idx - 3 : idx + 3])


def extract_pseU_pus7_motif(seq, idx, strand):
    """提取pseU pus7 motif"""
    if strand == "+":
        return seq[idx - 2 : idx + 3]
    else:
        return get_reverse_complementary_sequence(seq[idx - 2 : idx + 3])


def extract_m6A_motif(seq, idx, strand):
    """提取m6A motif"""
    if strand == "+":
        return seq[idx - 2 : idx + 3]
    else:
        return get_reverse_complementary_sequence(seq[idx - 2 : idx + 3])


def classify_motif_simple(match, region):
    """简单的motif分类 - 返回匹配的序列"""
    return match.group()


def classify_m5C_motif(match, region, pattern_index):
    """m5C motif分类"""
    if len(match.group()) == 3 and match.group()[2] in "ATC":
        return "CHH"
    elif len(match.group()) == 3 and match.group()[2] == "G":
        return "CHG"
    else:
        return "CG"


def classify_pseU_motif(match, region, pattern_index):
    """pseU motif分类"""
    if pattern_index == 0:
        return "pus4"
    elif pattern_index == 1:
        return "pus7"
    else:
        return "pus1"


def find_position(cigar_tuples, read_start, target_index):
    # op:M0 I1 N3 D2 S4
    # return 0-base position
    current_index = 0
    current_position = read_start
    for cigar in cigar_tuples:
        op, length = cigar
        if op in [0, 1, 4]:
            current_index += length
        if current_index <= target_index:
            if op in [0, 2, 3]:
                current_position += length
        else:
            if op == 0:
                return current_position + (target_index - (current_index - length))
            else:
                return None


def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans("ATGCatgc", "TACGtacg")
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq


def process_single_read(read, mod_threshold, min_qscore, min_mapq, strand, mod_tuple, base_dict, motif_info):
    """处理单个read,返回修饰位点集合"""
    results = set()
    
    read_id = read.query_name
    read_chr = read.reference_name
    read_seq = read.query_sequence
    read_qscore = read.query_qualities

    # read筛选:确保read中含有修饰信息
    if read.modified_bases is None:
        return results
    
    # read筛选:MAPQ质量筛选
    if read.mapping_quality < min_mapq:
        return results
    
    # read筛选:read都位于常染色体
    if not (
        read_chr.startswith("chr")
        and read_chr[3:].isdigit()
        and 1 <= int(read_chr[3:]) <= 22
    ):
        return results

    # 获取MM标签,提取修饰位点的位置和概率信息
    mm_tag = read.modified_bases.get(mod_tuple, [])
    if len(mm_tag) == 0:
        return results
    
    for tup in mm_tag:
        if read_qscore[tup[0]] < min_qscore:
            continue
        if motif_info is None:
            # 碱基模式
            if base_dict and read_seq[tup[0]] != base_dict[strand]:
                continue
        else:
            # motif模式
            motif = None
            motif_classification = ""
            for i, (pattern, extractor) in enumerate(
                zip(motif_info["patterns"], motif_info["extractors"])
            ):
                motif_region = extractor(read_seq, tup[0], strand)
                match = pattern.match(motif_region)
                if match:
                    motif = motif_info["classifier"](match, motif_region)
                    if "classifier2" in motif_info:
                        motif_classification = motif_info["classifier2"](
                            match, motif_region, i
                        )
                    break
            if not motif:
                continue
        
        if tup[1] >= (mod_threshold * 256):
            target_index = tup[0]
            geno_pos = find_position(
                read.cigartuples, read.reference_start, target_index
            )
            if geno_pos:
                if motif_info is None:
                    results.add(
                        f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}"
                    )
                else:
                    if motif_classification:
                        results.add(
                            f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}\t{motif}\t{motif_classification}"
                        )
                    else:
                        results.add(
                            f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}\t{motif}"
                        )
    
    return results


def process_chromosome_chunk(args):
    """
    处理单个染色体的修饰位点
    返回该染色体的修饰位点集合
    """
    chrom, map_bam_file, mod_threshold, min_qscore, min_mapq, strand, mod_tuple, base_dict, motif_info = args
    
    haplotype_set = set()
    
    # 重定向stderr到devnull
    original_stderr = sys.stderr
    sys.stderr = open(os.devnull, "w")
    
    try:
        with pysam.AlignmentFile(map_bam_file, "rb") as bamfile:
            # 只处理指定染色体
            for read in bamfile.fetch(chrom):
                read_results = process_single_read(
                    read, mod_threshold, min_qscore, min_mapq, 
                    strand, mod_tuple, base_dict, motif_info
                )
                haplotype_set.update(read_results)
    
    except Exception as e:
        sys.stderr = original_stderr
        print(f"Warning: Error processing {chrom}: {e}", file=sys.stderr)
        return set()
    
    finally:
        sys.stderr.close()
        sys.stderr = original_stderr
    
    return haplotype_set


def get_read_modify_parallel(
    map_bam_file,
    mod_threshold,
    min_qscore,
    min_mapq,
    threads,
    strand,
    mod_tuple,
    base_dict=None,
    motif_info=None,
):

    valid_chroms = [f"chr{i}" for i in range(1, 23)]
    num_tasks = len(valid_chroms)  # 22个任务
    
    # 智能线程数调整
    effective_threads = threads
    if threads > num_tasks:
        print(f"\n⚠️  警告: 指定线程数({threads})大于任务所需，自动调整为{num_tasks}个线程以避免资源浪费)")
        effective_threads = num_tasks
    
    print(f"开始扫描修饰位点 (链: {strand})")
    print(f"  - 指定线程数: {threads}")
    print(f"  - 实际使用线程数: {effective_threads}")
    print(f"  - 处理染色体: chr1-chr22 (共{num_tasks}个任务)")
    print(f"  - 修饰阈值: {mod_threshold}")
    print(f"  - 最小碱基质量: {min_qscore}")
    print(f"  - 最小MAPQ: {min_mapq}")
    
    # 准备参数列表
    args_list = [
        (chrom, map_bam_file, mod_threshold, min_qscore, min_mapq,
         strand, mod_tuple, base_dict, motif_info)
        for chrom in valid_chroms
    ]
    
    # 并行处理染色体
    final_set = set()
    
    with Pool(processes=effective_threads) as pool:
        # 使用imap_unordered并添加进度条
        results = list(tqdm(
            pool.imap_unordered(process_chromosome_chunk, args_list),
            total=len(args_list),
            desc="  扫描染色体",
            unit="chr",
            ncols=80
        ))
    
    # 合并所有染色体的结果
    for result_set in results:
        final_set.update(result_set)
    
    print(f"\n✓ 扫描完成! 检测到 {len(final_set):,} 个修饰位点")
    
    return final_set


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="区分正负链,获取修饰位点",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-o", "--output_prefix", type=str, help="Output folder and file prefix"
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        help="BAM file path.(上一步生成的*_calls_sorted_map0.bam or *_calls_sorted_map16.bam)",
    )
    parser.add_argument(
        "-m",
        "--mod_type",
        type=str,
        choices=["m6A", "inosine", "m5C", "pseU"],
        help="Modification Type.(4 types: m6A, inosine(A to I), m5C, pseU(pseudouridine))",
    )
    parser.add_argument(
        "-f",
        "--mod_threshold",
        type=float,
        default=0.75,
        help="Modification Threshold(default=0.75)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use.(default=4, max effective=22 for single strand)",
    )
    parser.add_argument(
        "-q",
        "--min_qscore",
        type=int,
        default=10,
        help="Base Min Query Quality(default=10)",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=0,
        help="Minimum Read MAPQ score (default=0)",
    )
    parser.add_argument("-s", "--strand", type=str, help="Strand")
    parser.add_argument(
        "--motif",
        help="Motif Context.\nIf not provided, the base context will be used.",
        action="store_true",
    )
    args = parser.parse_args()

    try:
        # 修饰类型配置 - 使用普通函数以支持多进程pickle
        inosine_base_dict = {"+": "A", "-": "T"}
        inosine_motif_info = {
            "patterns": [re.compile(r"TA[GT]")],
            "extractors": [extract_inosine_motif],
            "classifier": classify_motif_simple,
        }

        m5C_base_dict = {"+": "C", "-": "G"}
        m5C_motif_info = {
            "patterns": [
                re.compile(r"C[ATC][ATC]"),  # CHH
                re.compile(r"C[ATC]G"),  # CHG
                re.compile(r"CG"),  # CG
            ],
            "extractors": [
                extract_m5C_motif_3bp,  # CHH
                extract_m5C_motif_3bp,  # CHG
                extract_m5C_motif_2bp,  # CG
            ],
            "classifier": classify_motif_simple,
            "classifier2": classify_m5C_motif,
        }

        pseU_base_dict = {"+": "T", "-": "A"}
        pseU_motif_info = {
            "patterns": [
                re.compile(r"GTTC[ATCG]A"),  # pus4 (索引 0)
                re.compile(r"TGTA[AG]"),  # pus7 (索引 1)
                re.compile(r"[ACT][AG]T"),  # pus1 (索引 2)
            ],
            "extractors": [
                extract_pseU_pus1_motif,
                extract_pseU_pus4_motif,
                extract_pseU_pus7_motif,
            ],
            "classifier": classify_motif_simple,
            "classifier2": classify_pseU_motif,
        }

        m6A_base_dict = {"+": "A", "-": "T"}
        m6A_motif_info = {
            "patterns": [re.compile(r"[GAT][GA]AC[ATC]")],
            "extractors": [extract_m6A_motif],
            "classifier": classify_motif_simple,
        }

        strand_dict = {"+": 0, "-": 1}
        mod_dict = {
            "m6A": ("A", strand_dict[args.strand], "a"),
            "inosine": ("A", strand_dict[args.strand], 17596),
            "m5C": ("C", strand_dict[args.strand], "m"),
            "pseU": ("T", strand_dict[args.strand], 17802),
        }
        mod_info_summary = {
            "m6A": m6A_motif_info,
            "inosine": inosine_motif_info,
            "m5C": m5C_motif_info,
            "pseU": pseU_motif_info,
        }
        mod_base_dict_summary = {
            "m6A": m6A_base_dict,
            "inosine": inosine_base_dict,
            "m5C": m5C_base_dict,
            "pseU": pseU_base_dict,
        }

        # 执行并行处理
        if args.motif:
            n = f"_{args.mod_type}_motif"
            read_mod_set = get_read_modify_parallel(
                args.bam,
                args.mod_threshold,
                args.min_qscore,
                args.min_mapq,
                args.threads,
                args.strand,
                mod_dict[args.mod_type],
                motif_info=mod_info_summary[args.mod_type],
            )
        else:
            n = f"_{args.mod_type}_base"
            read_mod_set = get_read_modify_parallel(
                args.bam,
                args.mod_threshold,
                args.min_qscore,
                args.min_mapq,
                args.threads,
                args.strand,
                mod_dict[args.mod_type],
                base_dict=mod_base_dict_summary[args.mod_type],
            )

        read_mod_file = (
            f"{args.output_prefix}_{'f' if args.strand == '+' else 'r'}{n}_tmp.txt"
        )
        
        print(f"\n正在保存结果到: {read_mod_file}")
        with open(read_mod_file, "w") as file:
            for item in read_mod_set:
                file.write(f"{item}\n")
        
        print(f"✓ 完成! 结果已保存")
        print("="*60)
        pass

    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)