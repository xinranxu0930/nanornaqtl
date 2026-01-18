#!/usr/bin/env python3
"""
多线程并行处理polyA尾长QTL识别脚本
并行执行callplQTL.py处理所有染色体和链的组合（共44个任务）
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import time


def setup_argparse():
    """设置命令行参数"""
    parser = argparse.ArgumentParser(
        description="并行处理polyA尾长QTL识别任务"
    )

    # 必需参数
    parser.add_argument(
        "-b", "--bam", type=str, required=True,
        help="BAM文件路径 (例如: *_calls_sorted_map.bam)"
    )
    parser.add_argument(
        "--snp_info", type=str, required=True, 
        help="SNP信息文件路径"
    )
    parser.add_argument(
        "-o", "--output_prefix", type=str, required=True, 
        help="输出文件前缀"
    )
    parser.add_argument(
        "-csv", "--polya_csv", type=str, required=True,
        help="*_polyAlen_result.csv路径"
    )
    parser.add_argument(
        "--geno_size", type=str, required=True,
        help="基因组大小文件路径 (例如：hg19.chrom.sizes)"
    )

    # 可选参数
    parser.add_argument(
        "-q", "--min_qscore", type=int, default=10, 
        help="最小碱基质量分数 (默认: 10)"
    )
    parser.add_argument(
        "--min_coverage", type=int, default=8, 
        help="最小总覆盖度 (默认: 8)"
    )
    parser.add_argument(
        "--mcmc_samples", type=int, default=2000, 
        help="MCMC采样数 (默认: 2000)"
    )
    parser.add_argument(
        "--threads", type=int, default=4, 
        help="并行线程数 (默认: 4)"
    )
    parser.add_argument(
        "--keep_tmp", action="store_true", 
        help="保留临时文件"
    )

    return parser


def get_bam_files(base_bam_path):
    """获取正负链BAM文件路径"""
    plus_strand_bam = base_bam_path.replace(".bam", "0.bam")
    minus_strand_bam = base_bam_path.replace(".bam", "16.bam")

    if not os.path.exists(plus_strand_bam):
        raise FileNotFoundError(f"正链BAM文件不存在: {plus_strand_bam}")
    if not os.path.exists(minus_strand_bam):
        raise FileNotFoundError(f"负链BAM文件不存在: {minus_strand_bam}")

    return plus_strand_bam, minus_strand_bam


def generate_tasks(args):
    """
    生成所有polyA QTL分析任务（22条染色体 × 2条链 = 44个任务）
    
    Returns:
        list: 任务列表，每个元素为 (chrom, strand, bam_file, args)
    """
    plus_strand_bam, minus_strand_bam = get_bam_files(args.bam)
    
    tasks = []
    chromosomes = [f"chr{i}" for i in range(1, 23)]
    
    for chrom in chromosomes:
        # 正链任务
        tasks.append((chrom, "+", plus_strand_bam, args))
        # 负链任务
        tasks.append((chrom, "-", minus_strand_bam, args))
    
    return tasks


def run_single_task(task_info):
    """
    运行单个polyA QTL分析任务
    
    Args:
        task_info: (chrom, strand, bam_file, args)
    
    Returns:
        tuple: (task_id, success, output_file_or_error)
    """
    chrom, strand, bam_file, args = task_info
    task_id = f"{chrom}_{strand}"
    
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        callplqtl_script = os.path.join(current_dir, "callplQTL.py")
        
        # 构建命令
        cmd = [
            "python", callplqtl_script,
            "-b", bam_file,
            "--snp_info", args.snp_info,
            "-o", args.output_prefix,
            "-csv", args.polya_csv,
            "-c", chrom,
            "-s", strand,
            "--geno_size", args.geno_size,
            "-q", str(args.min_qscore),
            "--min_coverage", str(args.min_coverage),
            "--mcmc_samples", str(args.mcmc_samples)
        ]
        
        # 执行命令
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        if result.returncode == 0:
            # 构建预期输出文件路径
            output_file = f"{args.output_prefix}_polyAlen_haplotype_{chrom}_{strand}_tmp.csv"
            
            if os.path.exists(output_file):
                return task_id, True, output_file
            else:
                return task_id, False, "输出文件不存在"
        else:
            error_msg = result.stderr[:200] if result.stderr else "未知错误"
            return task_id, False, error_msg
            
    except Exception as e:
        return task_id, False, str(e)


def parallel_run(tasks, n_threads):
    """
    并行执行所有任务
    
    Returns:
        tuple: (successful_files, failed_tasks)
    """
    successful_files = []
    failed_tasks = []
    
    print(f"使用 {n_threads} 个线程并行处理 {len(tasks)} 个任务...")
    print("-" * 50)
    
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        future_to_task = {
            executor.submit(run_single_task, task): f"{task[0]}_{task[1]}"
            for task in tasks
        }
        
        completed = 0
        for future in as_completed(future_to_task):
            task_id = future_to_task[future]
            completed += 1
            
            try:
                result_task_id, success, result = future.result()
                
                if success:
                    successful_files.append(result)
                    print(f"[{completed}/{len(tasks)}] {task_id} 完成")
                else:
                    failed_tasks.append((task_id, result))
                    print(f"[{completed}/{len(tasks)}] {task_id} 失败: {result}")
                    
            except Exception as e:
                failed_tasks.append((task_id, str(e)))
                print(f"[{completed}/{len(tasks)}] {task_id} 异常: {str(e)}")
    
    return successful_files, failed_tasks


def merge_results(successful_files, output_prefix, keep_tmp):
    """
    合并所有成功的结果文件
    
    Returns:
        str: 最终合并文件的路径，失败返回None
    """
    if not successful_files:
        print("没有成功的结果文件需要合并")
        return None
    
    print("-" * 50)
    print(f"合并 {len(successful_files)} 个结果文件...")
    
    all_dataframes = []
    
    for file_path in successful_files:
        try:
            df = pd.read_csv(file_path)
            if len(df) > 0:
                all_dataframes.append(df)
        except Exception as e:
            print(f"读取文件失败: {file_path}, 错误: {str(e)}")
    
    if not all_dataframes:
        print("没有有效的数据可以合并")
        return None
    
    # 合并所有数据
    merged_df = pd.concat(all_dataframes, ignore_index=True)
    
    # 按染色体和位置排序
    try:
        merged_df["chrom_num"] = merged_df["chrom"].str.replace("chr", "").astype(int)
        merged_df = merged_df.sort_values(
            by=["chrom_num", "snp_pos_1base"], ascending=[True, True]
        )
        merged_df = merged_df.drop("chrom_num", axis=1)
    except Exception:
        merged_df = merged_df.sort_values(by=["chrom", "snp_pos_1base"])
    
    # 生成最终文件名并保存
    final_output = f"{output_prefix}_polyA_tail_length_QTLs_result.csv"
    merged_df.to_csv(final_output, index=False)
    print(f"合并完成，总行数: {len(merged_df)}")
    print(f"结果保存至: {final_output}")
    
    # 删除临时文件
    if not keep_tmp:
        print("删除临时文件...")
        for file_path in successful_files:
            try:
                os.remove(file_path)
            except Exception:
                pass
    
    return final_output


def main():
    """主函数"""
    parser = setup_argparse()
    args = parser.parse_args()
    
    print("=" * 50)
    print("polyA尾长QTL并行分析")
    print("=" * 50)
    
    start_time = time.time()
    
    # 验证输入文件
    required_files = [args.snp_info, args.polya_csv, args.geno_size]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"错误: 输入文件不存在: {file_path}")
            sys.exit(1)
    
    # 验证callplQTL.py脚本
    current_dir = os.path.dirname(os.path.abspath(__file__))
    callplqtl_script = os.path.join(current_dir, "callplQTL.py")
    if not os.path.exists(callplqtl_script):
        print(f"错误: 脚本文件不存在: {callplqtl_script}")
        sys.exit(1)
    
    # 验证BAM文件
    try:
        plus_strand_bam, minus_strand_bam = get_bam_files(args.bam)
        print(f"正链BAM: {plus_strand_bam}")
        print(f"负链BAM: {minus_strand_bam}")
    except FileNotFoundError as e:
        print(f"错误: {str(e)}")
        sys.exit(1)
    
    # 生成任务
    tasks = generate_tasks(args)
    n_tasks = len(tasks)  # 44
    
    # 确定线程数
    n_threads = args.threads
    if n_threads > n_tasks:
        print(f"提示: 线程数({n_threads})超过任务数({n_tasks})，将使用{n_tasks}个线程")
        n_threads = n_tasks
    
    print(f"总任务数: {n_tasks}")
    print(f"MCMC采样数: {args.mcmc_samples}")
    print(f"最小覆盖度: {args.min_coverage}")
    print("-" * 50)
    
    # 并行执行
    successful_files, failed_tasks = parallel_run(tasks, n_threads)
    
    # 输出统计
    print("-" * 50)
    print(f"成功: {len(successful_files)}/{n_tasks}")
    print(f"失败: {len(failed_tasks)}/{n_tasks}")
    
    if failed_tasks:
        print("失败的任务:")
        for task_id, error in failed_tasks:
            print(f"  {task_id}: {error}")
    
    # 合并结果
    if successful_files:
        final_file = merge_results(successful_files, args.output_prefix, args.keep_tmp)
        
        if final_file:
            # 显示简要统计
            try:
                final_df = pd.read_csv(final_file)
                print("-" * 50)
                print(f"结果统计:")
                print(f"  总SNP数: {len(final_df)}")
                print(f"  染色体数: {final_df['chrom'].nunique()}")
                
                if 'BF' in final_df.columns:
                    sig_bf3 = len(final_df[final_df['BF'] > 3])
                    sig_bf10 = len(final_df[final_df['BF'] > 10])
                    print(f"  BF>3的候选QTL: {sig_bf3}")
                    print(f"  BF>10的候选QTL: {sig_bf10}")
                
            except Exception:
                pass
    else:
        print("错误: 没有成功的结果，无法生成最终文件")
        sys.exit(1)
    
    end_time = time.time()
    total_time = end_time - start_time
    print("=" * 50)
    print(f"总耗时: {total_time:.1f}秒 ({total_time/60:.1f}分钟)")
    print("=" * 50)


if __name__ == "__main__":
    main()