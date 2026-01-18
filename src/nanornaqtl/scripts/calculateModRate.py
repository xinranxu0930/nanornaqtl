import pandas as pd
import pickle
import argparse
import pysam
from multiprocessing import Pool
from collections import Counter
from tqdm import tqdm
import sys
import os


def generate_reads_dict(pileup_df):
    """
    第一步: 直接从pileup数据生成reads字典
    """
    print("\n[1/4] 生成reads字典...")
    
    reads_dict = {}
    
    # 按(chrom, pos, strand)分组
    grouped = pileup_df.groupby(['chrom', 'pos_0base', 'strand'])['read_id']
    
    for (chrom, pos, strand), read_ids in tqdm(grouped, desc="  处理位点", unit="site"):
        key = f"{chrom}_{pos}_{strand}"
        reads_dict[key] = ';'.join(read_ids.tolist())
    
    print(f"✓ 生成完成: {len(reads_dict):,} 个位点")
    return reads_dict


def aggregate_modification_sites(pileup_df, motif_cols):
    """
    第二步: 统计每个位点的修饰信息
    """
    print("\n[2/4] 统计修饰位点信息...")
    
    # 构建聚合字典
    agg_dict = {
        'read_id': lambda x: list(set(x))  # 去重read_id
    }
    
    # 处理motif列 - 选择出现最多的
    for col in motif_cols:
        if col in pileup_df.columns:
            agg_dict[col] = lambda x, col=col: Counter(x).most_common(1)[0][0] if len(x) > 0 else ''
    
    # 按(chrom, pos, strand)分组聚合
    print("  按位点分组统计...")
    grouped = pileup_df.groupby(['chrom', 'pos_0base', 'strand']).agg(agg_dict).reset_index()
    
    # 计算修饰read数量
    grouped['mod_num'] = grouped['read_id'].apply(len)
    
    # 不再需要read_id列表
    grouped = grouped.drop(columns=['read_id'])
    
    print(f"✓ 统计完成: {len(grouped):,} 个唯一位点")
    return grouped


def process_chromosome_coverage(args):
    """
    处理单个染色体的覆盖度计算
    """
    chrom, positions, bam_file, min_mapq = args
    
    coverage_dict = {}
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            position_set = set(positions)
            
            # 获取染色体范围
            if not positions:
                return {chrom: {}}
            
            min_pos = min(positions)
            max_pos = max(positions)
            
            # pileup扫描
            for pileup_column in bam.pileup(
                chrom, min_pos, max_pos + 1, 
                min_mapping_quality=min_mapq,
                stepper='nofilter',
                ignore_overlaps=False,
                ignore_orphans=False
            ):
                if pileup_column.pos in position_set:
                    # 统计该位点通过MAPQ阈值的read数量
                    coverage = sum(1 for read in pileup_column.pileups 
                                 if not read.is_del and not read.is_refskip)
                    coverage_dict[pileup_column.pos] = coverage
    
    except Exception as e:
        print(f"警告: 处理染色体 {chrom} 时出错: {e}", file=sys.stderr)
        return {chrom: {}}
    
    return {chrom: coverage_dict}


def calculate_coverage_parallel(site_stats_df, bam_file, min_mapq, threads):
    """
    第三步: 并行计算每个位点的覆盖度
    """
    print("\n[3/4] 计算覆盖度...")
    
    # 按染色体收集需要计算覆盖度的位点
    positions_by_chrom = {}
    for chrom in site_stats_df['chrom'].unique():
        positions = site_stats_df[site_stats_df['chrom'] == chrom]['pos_0base'].unique().tolist()
        positions_by_chrom[chrom] = sorted(positions)
    
    num_chroms = len(positions_by_chrom)
    
    # 智能线程数管理
    effective_threads = threads
    if threads > num_chroms:
        print(f"\n⚠️  警告: 指定线程数({threads})大于染色体数({num_chroms})")
        print(f"    自动调整为{num_chroms}个线程以避免资源浪费\n")
        effective_threads = num_chroms
    
    print(f"  - 染色体数量: {num_chroms}")
    print(f"  - 指定线程数: {threads}")
    print(f"  - 实际使用线程数: {effective_threads}")
    print(f"  - 总位点数: {len(site_stats_df):,}")
    
    # 准备参数
    args_list = [
        (chrom, positions, bam_file, min_mapq)
        for chrom, positions in positions_by_chrom.items()
    ]
    
    # 并行计算覆盖度
    coverage_map = {}
    with Pool(processes=effective_threads) as pool:
        results = list(tqdm(
            pool.imap_unordered(process_chromosome_coverage, args_list),
            total=len(args_list),
            desc="  处理染色体",
            unit="chr"
        ))
    
    # 合并结果
    for result in results:
        coverage_map.update(result)
    
    # 添加覆盖度到DataFrame
    def get_coverage(row):
        chrom_cov = coverage_map.get(row['chrom'], {})
        return chrom_cov.get(row['pos_0base'], 0)
    
    site_stats_df['cov'] = site_stats_df.apply(get_coverage, axis=1)
    
    print(f"✓ 覆盖度计算完成")
    return site_stats_df


def calculate_mod_rate_and_filter(site_stats_df, min_rate, min_cov):
    """
    第四步: 计算修饰率并过滤
    """
    print("\n[4/4] 计算修饰率并过滤...")
    
    # 计算修饰率
    site_stats_df['mod_rate'] = site_stats_df.apply(
        lambda row: round(row['mod_num'] / row['cov'], 4) if row['cov'] > 0 else 0.0,
        axis=1
    )
    
    # 数据验证
    invalid_sites = site_stats_df[site_stats_df['mod_num'] > site_stats_df['cov']]
    if len(invalid_sites) > 0:
        print(f"⚠️  警告: 发现 {len(invalid_sites)} 个异常位点(mod_num > cov),已自动修正")
        # 修正异常数据
        site_stats_df.loc[site_stats_df['mod_num'] > site_stats_df['cov'], 'mod_num'] = \
            site_stats_df.loc[site_stats_df['mod_num'] > site_stats_df['cov'], 'cov']
        # 重新计算修饰率
        site_stats_df['mod_rate'] = site_stats_df.apply(
            lambda row: round(row['mod_num'] / row['cov'], 4) if row['cov'] > 0 else 0.0,
            axis=1
        )
    
    # 过滤
    print(f"  过滤前: {len(site_stats_df):,} 个位点")
    filtered_df = site_stats_df[
        (site_stats_df['mod_rate'] >= min_rate) & 
        (site_stats_df['cov'] >= min_cov)
    ].copy()
    print(f"  过滤后: {len(filtered_df):,} 个位点")
    
    # 添加1-based位置
    filtered_df['pos_1base'] = filtered_df['pos_0base'] + 1
    
    # 检查修饰率范围
    if len(filtered_df) > 0:
        min_mod_rate = filtered_df['mod_rate'].min()
        max_mod_rate = filtered_df['mod_rate'].max()
        print(f"  修饰率范围: {min_mod_rate:.4f} - {max_mod_rate:.4f}")
        
        if max_mod_rate > 1.0:
            print(f"⚠️  错误: 存在修饰率>1的异常数据!", file=sys.stderr)
    
    print(f"✓ 过滤完成")
    return filtered_df


def determine_file_columns(pileup_file, modification_type):
    """
    检测pileup文件的列结构
    
    Returns:
    --------
    tuple : (input_cols, output_cols, motif_cols)
    """
    try:
        # 读取第一行检测列数
        with open(pileup_file, 'r') as f:
            first_line = f.readline().strip()
            n_cols = len(first_line.split('\t')) if first_line else 4
        
        # 基础列
        input_cols = ["chrom", "pos_0base", "strand", "read_id"]
        output_cols = ["chrom", "pos_1base", "strand", "mod_num", "cov", "mod_rate"]
        motif_cols = []
        
        # 根据修饰类型和列数确定motif列
        if n_cols == 5:
            input_cols.append('motif')
            output_cols.append('motif')
            motif_cols.append('motif')
        elif n_cols == 6:
            input_cols.extend(['motif', 'motif_classification'])
            output_cols.extend(['motif', 'motif_classification'])
            motif_cols.extend(['motif', 'motif_classification'])
        
        print(f"检测到 {n_cols} 列, motif列: {motif_cols if motif_cols else '无'}")
        return input_cols, output_cols, motif_cols
        
    except Exception as e:
        print(f"警告: 列结构检测失败 ({e}), 使用默认配置", file=sys.stderr)
        return ["chrom", "pos_0base", "strand", "read_id"], \
               ["chrom", "pos_1base", "strand", "mod_num", "cov", "mod_rate"], \
               []


def process_modifications(
    bam_file,
    pileup_file,
    output_prefix,
    strand,
    modification_type,
    min_rate,
    min_cov,
    min_mapq,
    threads
):
    """
    主处理函数
    """
    print("=" * 60)
    print(f"修饰位点分析: {modification_type}")
    print(f"  链方向: {strand}")
    print(f"  最小修饰率: {min_rate}")
    print(f"  最小覆盖度: {min_cov}")
    print(f"  最小MAPQ: {min_mapq}")
    print(f"  线程数: {threads}")
    print("=" * 60)
    
    # 检查输入文件
    for file_path, name in [(bam_file, "BAM"), (pileup_file, "Pileup")]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{name}文件不存在: {file_path}")
    
    # 检测文件列结构
    input_cols, output_cols, motif_cols = determine_file_columns(pileup_file, modification_type)
    
    # 加载pileup数据
    print(f"\n加载pileup文件: {pileup_file}")
    pileup_df = pd.read_csv(pileup_file, header=None, sep='\t', dtype=str)
    
    if pileup_df.empty:
        print("警告: pileup文件为空")
        return create_empty_outputs(output_prefix, modification_type, strand, output_cols)
    
    # 设置列名
    pileup_df.columns = input_cols[:len(pileup_df.columns)]
    pileup_df['pos_0base'] = pd.to_numeric(pileup_df['pos_0base'], errors='coerce')
    pileup_df = pileup_df.dropna(subset=['pos_0base'])
    pileup_df['pos_0base'] = pileup_df['pos_0base'].astype(int)
    
    print(f"总记录数: {len(pileup_df):,}")
    
    # 筛选指定链的数据
    strand_df = pileup_df[pileup_df['strand'] == strand].copy()
    
    if strand_df.empty:
        print(f"警告: 链 {strand} 没有数据")
        return create_empty_outputs(output_prefix, modification_type, strand, output_cols)
    
    print(f"链 {strand} 记录数: {len(strand_df):,}")
    
    # 第一步: 生成reads字典 (PKL文件)
    reads_dict = generate_reads_dict(strand_df)
    
    # 第二步: 统计修饰位点信息(去重+motif选择)
    site_stats_df = aggregate_modification_sites(strand_df, motif_cols)
    
    # 第三步: 计算覆盖度
    site_stats_df = calculate_coverage_parallel(site_stats_df, bam_file, min_mapq, threads)
    
    # 第四步: 计算修饰率并过滤
    final_df = calculate_mod_rate_and_filter(site_stats_df, min_rate, min_cov)
    
    # 保存结果
    sites_file, reads_dict_file = save_results(
        final_df, reads_dict, output_prefix, modification_type, strand, output_cols
    )
    
    print("\n" + "=" * 60)
    print("处理完成!")
    print(f"  修饰位点文件: {sites_file}")
    print(f"  Reads字典文件: {reads_dict_file}")
    print("=" * 60)
    
    return sites_file, reads_dict_file


def create_empty_outputs(output_prefix, modification_type, strand, output_cols):
    """创建空输出文件"""
    sites_file = f"{output_prefix}_{modification_type}_sites_{strand}_tmp.csv"
    reads_dict_file = f"{output_prefix}_{modification_type}_reads_{strand}_tmp.pkl"
    
    pd.DataFrame(columns=output_cols).to_csv(sites_file, index=False)
    with open(reads_dict_file, 'wb') as f:
        pickle.dump({}, f)
    
    return sites_file, reads_dict_file


def save_results(final_df, reads_dict, output_prefix, modification_type, strand, output_cols):
    """保存结果文件"""
    sites_file = f"{output_prefix}_{modification_type}_sites_{strand}_tmp.csv"
    reads_dict_file = f"{output_prefix}_{modification_type}_reads_{strand}_tmp.pkl"
    
    # 保存sites文件
    if not final_df.empty:
        final_df = final_df.reindex(columns=output_cols, fill_value='')
    else:
        final_df = pd.DataFrame(columns=output_cols)
    
    final_df.to_csv(sites_file, index=False)
    
    # 保存reads字典
    with open(reads_dict_file, 'wb') as f:
        pickle.dump(reads_dict, f)
    
    return sites_file, reads_dict_file


def main():
    parser = argparse.ArgumentParser(
        description="计算修饰位点的修饰率 - 按染色体并行优化版本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python calculateModRate.py -b calls_sorted_map0.bam -f f_m6A_base_tmp.txt -o output -s + -m m6A -t 20
  python calculateModRate.py -b calls_sorted_map16.bam -f r_m5C_motif_tmp.txt -o output -s - -m m5C -t 22

说明:
  - 最多使用22个线程(对应22个染色体)
  - 自动去重修饰位点
  - 保证 mod_num <= cov
  - Motif冲突时选择出现最多的
        """
    )
    
    parser.add_argument("-b", "--bam", required=True, help="BAM文件路径")
    parser.add_argument("-f", "--pileup_file", required=True, help="Pileup文件路径")
    parser.add_argument("-o", "--output_prefix", required=True, help="输出文件前缀")
    parser.add_argument("-s", "--strand", required=True, choices=["+", "-"], help="链方向")
    parser.add_argument("-m", "--modification_type", required=True,
                       choices=["m6A", "m5C", "pseU", "inosine"], help="修饰类型")
    parser.add_argument("-r", "--min_rate", type=float, default=0.1, 
                       help="最小修饰率 (默认: 0.1)")
    parser.add_argument("-c", "--min_cov", type=int, default=5, 
                       help="最小覆盖度 (默认: 5)")
    parser.add_argument("--min_mapq", type=int, default=0, 
                       help="最小MAPQ阈值 (默认: 0)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                       help="线程数 (默认: 4, 最大有效值: 22)")
    
    args = parser.parse_args()
    
    try:
        sites_file, reads_dict_file = process_modifications(
            bam_file=args.bam,
            pileup_file=args.pileup_file,
            output_prefix=args.output_prefix,
            strand=args.strand,
            modification_type=args.modification_type,
            min_rate=args.min_rate,
            min_cov=args.min_cov,
            min_mapq=args.min_mapq,
            threads=args.threads
        )
        
    except Exception as e:
        print(f"\n错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()