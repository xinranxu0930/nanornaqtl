#!/usr/bin/env python3
"""
多线程并行处理QTL识别脚本
用于并行执行calliu3auQTL.py脚本处理所有染色体和链的组合
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import time
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("qtl_processing.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def setup_argparse():
    """设置命令行参数"""
    parser = argparse.ArgumentParser(description="并行处理QTL识别任务")

    # 必需参数
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="BAM文件路径 (例如: *_calls_sorted_map.bam)",
    )
    parser.add_argument("--snp_info", type=str, required=True, help="SNP信息文件路径")
    parser.add_argument(
        "-o", "--output_prefix", type=str, required=True, help="输出文件前缀"
    )
    parser.add_argument(
        "-f",
        "--read_overlap_file",
        type=str,
        help="read2apa(*_APA_result.csv) or read2isoform(OUT.transcript_model_reads.tsv.gz) overlap file path"
    )
    parser.add_argument(
        "--geno_size",
        type=str,
        required=True,
        help="基因组大小文件路径 (例如：hg19.chrom.sizes)",
    )
    parser.add_argument(
        "-m",
        "--molecular_type",
        type=str,
        required=True,
        choices=["APA", "isoform"],
        help="分子表型类型")


    # 可选参数
    parser.add_argument(
        "-q", "--min_qscore", type=int, default=10, help="最小碱基质量分数 (默认: 10)"
    )
    parser.add_argument(
        "--min_coverage", type=int, default=8, help="最小总覆盖度(default=8)",
    )
    parser.add_argument(
        "--mcmc_samples", type=int, default=1000, help="MCMC采样数(default=1000)",
    )
    parser.add_argument("-t", "--threads", type=int, default=4, help="并行线程数 (默认: 4)")
    parser.add_argument("--keep_tmp", action="store_true", help="保留临时文件")

    return parser


def get_bam_files(base_bam_path):
    plus_strand_bam = base_bam_path.replace(".bam", "0.bam")
    minus_strand_bam = base_bam_path.replace(".bam", "16.bam")

    # 检查文件是否存在
    if not os.path.exists(plus_strand_bam):
        raise FileNotFoundError(f"正链BAM文件不存在: {plus_strand_bam}")
    if not os.path.exists(minus_strand_bam):
        raise FileNotFoundError(f"负链BAM文件不存在: {minus_strand_bam}")

    return plus_strand_bam, minus_strand_bam


def generate_tasks(args):
    """
    生成所有需要处理的任务列表

    Returns:
        list: 任务列表，每个任务包含 (染色体, 链, BAM文件路径)
    """
    # 定义所有染色体
    chromosomes = [f"chr{i}" for i in range(1, 23)]  # chr1-chr22

    # 获取对应的BAM文件
    plus_strand_bam, minus_strand_bam = get_bam_files(args.bam)

    tasks = []

    # 生成所有染色体和链的组合
    for chrom in chromosomes:
        # 正链任务
        tasks.append((chrom, "+", plus_strand_bam))
        # 负链任务
        tasks.append((chrom, "-", minus_strand_bam))

    logger.info(f"总共生成了 {len(tasks)} 个处理任务")
    return tasks


def run_single_task(task_info):
    """
    运行单个QTL处理任务

    Args:
        task_info: 包含任务信息和参数的元组

    Returns:
        tuple: (任务标识, 是否成功, 输出文件路径或错误信息)
    """
    (chrom, strand, bam_file), args = task_info

    task_id = f"{chrom}_{strand}"

    try:
        # 获取当前脚本所在目录
        current_dir = os.path.dirname(os.path.abspath(__file__))
        callmodqtl_script = os.path.join(current_dir, "calliu3auQTL.py")

        # 构建命令
        cmd = [
            "python",
            callmodqtl_script,
            "-b",
            bam_file,
            "--snp_info",
            args.snp_info,
            "-o",
            args.output_prefix,
            "-f",
            args.read_overlap_file,
            "-c",
            chrom,
            "-s",
            strand,
            "--geno_size",
            args.geno_size,
            "-q",
            str(args.min_qscore),
            "-m",
            args.molecular_type,
            "--min_coverage",
            str(args.min_coverage),
            "--mcmc_samples",
            str(args.mcmc_samples)
        ]

        logger.info(f"开始处理任务: {task_id}")

        # 执行命令
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode == 0:
            # 构建预期的输出文件路径
            output_file = f"{args.output_prefix}_{args.molecular_type}_haplotype_{chrom}_{strand}_tmp.csv"

            if os.path.exists(output_file):
                logger.info(f"任务 {task_id} 完成，生成文件: {output_file}")
                return task_id, True, output_file
            else:
                logger.warning(f"任务 {task_id} 执行成功但未生成输出文件")
                return task_id, False, "输出文件不存在"
        else:
            error_msg = f"命令执行失败: {result.stderr}"
            logger.error(f"任务 {task_id} 失败: {error_msg}")
            return task_id, False, error_msg

    except Exception as e:
        error_msg = f"任务执行异常: {str(e)}"
        logger.error(f"任务 {task_id} 异常: {error_msg}")
        return task_id, False, error_msg


def merge_results(output_prefix, molecular_type, successful_files, keep_tmp=False):
    """
    合并所有成功的结果文件

    Args:
        output_prefix: 输出文件前缀
        molecular_type: 修饰类型
        successful_files: 成功生成的文件列表
        keep_tmp: 是否保留临时文件

    Returns:
        str: 最终合并文件的路径
    """
    if not successful_files:
        logger.warning("没有成功的结果文件需要合并")
        return None

    logger.info(f"开始合并 {len(successful_files)} 个结果文件")

    all_dataframes = []

    for file_path in successful_files:
        try:
            df = pd.read_csv(file_path)
            if len(df) > 0:
                all_dataframes.append(df)
                logger.debug(f"成功读取文件: {file_path}, 行数: {len(df)}")
            else:
                logger.warning(f"文件为空: {file_path}")
        except Exception as e:
            logger.error(f"读取文件失败: {file_path}, 错误: {str(e)}")

    if not all_dataframes:
        logger.warning("没有有效的数据可以合并")
        return None

    # 合并所有数据
    merged_df = pd.concat(all_dataframes, ignore_index=True)

    # 排序：按染色体、SNP位置、修饰位置排序
    try:
        # 提取染色体编号用于排序
        merged_df["chrom_num"] = merged_df["chrom"].str.replace("chr", "").astype(int)
        merged_df = merged_df.sort_values(
            by=["chrom_num", "snp_pos_1base"],
            ascending=[True, True],
        )
        merged_df = merged_df.drop("chrom_num", axis=1)
    except Exception as e:
        logger.warning(f"排序失败，使用默认排序: {str(e)}")
        merged_df = merged_df.sort_values(
            by=["chrom", "snp_pos_1base"], ascending=[True, True]
        )

    # 生成最终文件名
    final_output = f"{output_prefix}_{molecular_type}_QTLs_result.csv"

    # 保存合并结果
    merged_df.to_csv(final_output, index=False)
    logger.info(f"合并完成，总行数: {len(merged_df)}, 输出文件: {final_output}")

    # 删除临时文件
    if not keep_tmp:
        logger.info("开始删除临时文件...")
        deleted_count = 0
        for file_path in successful_files:
            try:
                os.remove(file_path)
                deleted_count += 1
                logger.debug(f"删除临时文件: {file_path}")
            except Exception as e:
                logger.error(f"删除文件失败: {file_path}, 错误: {str(e)}")
        logger.info(f"已删除 {deleted_count} 个临时文件")
    else:
        logger.info("保留临时文件")

    return final_output


def main():
    """主函数"""
    parser = setup_argparse()
    args = parser.parse_args()

    logger.info("=" * 50)
    logger.info("开始并行QTL处理任务")
    logger.info("=" * 50)

    start_time = time.time()

    try:
        # 验证输入文件
        required_files = [
            args.snp_info,
            args.read_overlap_file,
            args.geno_size,
        ]
        for file_path in required_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"输入文件不存在: {file_path}")

        # 验证calliu3auQTL.py脚本文件
        current_dir = os.path.dirname(os.path.abspath(__file__))
        callmodqtl_script = os.path.join(current_dir, "calliu3auQTL.py")
        if not os.path.exists(callmodqtl_script):
            raise FileNotFoundError(f"QTL脚本不存在: {callmodqtl_script}")

        # 生成任务列表
        tasks = generate_tasks(args)

        # 准备任务信息
        task_infos = [(task, args) for task in tasks]

        # 并行执行任务
        successful_files = []
        failed_tasks = []

        logger.info(f"使用 {args.threads} 个线程并行处理...")

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            # 提交所有任务
            future_to_task = {
                executor.submit(run_single_task, task_info): task_info[0]
                for task_info in task_infos
            }

            # 收集结果
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    task_id, success, result = future.result()

                    if success:
                        successful_files.append(result)
                    else:
                        failed_tasks.append((task_id, result))

                except Exception as e:
                    task_id = f"{task[0]}_{task[1]}"
                    failed_tasks.append((task_id, str(e)))
                    logger.error(f"任务 {task_id} 处理异常: {str(e)}")

        # 输出统计信息
        logger.info("=" * 50)
        logger.info("任务执行完成统计:")
        logger.info(f"总任务数: {len(tasks)}")
        logger.info(f"成功任务数: {len(successful_files)}")
        logger.info(f"失败任务数: {len(failed_tasks)}")

        if failed_tasks:
            logger.warning("失败的任务:")
            for task_id, error in failed_tasks:
                logger.warning(f"  {task_id}: {error}")

        # 合并结果
        if successful_files:
            final_file = merge_results(
                args.output_prefix,
                args.molecular_type,
                successful_files,
                args.keep_tmp,
            )

            if final_file:
                logger.info("=" * 50)
                logger.info(f"处理完成！最终结果文件: {final_file}")

                # 显示文件统计信息
                try:
                    final_df = pd.read_csv(final_file)
                    logger.info(f"最终文件包含 {len(final_df)} 行数据")
                    logger.info(f"涉及 {final_df['chrom'].nunique()} 个染色体")
                    logger.info(f"涉及 {final_df['SNP'].nunique()} 个SNP位点")
                except Exception as e:
                    logger.error(f"读取最终文件统计信息失败: {str(e)}")
            else:
                logger.error("合并结果失败")
                sys.exit(1)
        else:
            logger.error("没有成功的结果文件，无法生成最终结果")
            sys.exit(1)

        end_time = time.time()
        total_time = end_time - start_time
        logger.info(f"总处理时间: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")

    except Exception as e:
        logger.error(f"程序执行失败: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
