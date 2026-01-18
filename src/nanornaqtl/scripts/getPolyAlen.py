import pysam
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Tuple


# 22条常染色体
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def process_chromosome(bamfile: str, chrom: str) -> List[Tuple[str, int]]:
    """
    处理单个染色体，提取所有read的ID和polyA长度
    
    Args:
        bamfile: BAM文件路径
        chrom: 染色体名称
    
    Returns:
        包含(readID, polyA_length)元组的列表
    """
    read_polyA_pairs = []

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        # 检查染色体是否存在于BAM文件中
        if chrom not in bam.references:
            return read_polyA_pairs
        
        for read in bam.fetch(chrom):
            try:
                pt_tag = read.get_tag("pt")
                if pt_tag is not None and pt_tag > 0:
                    read_polyA_pairs.append((read.query_name, pt_tag))
            except KeyError:
                # pt tag不存在，跳过
                continue

    return read_polyA_pairs


def get_polyA_len(bamfile: str, output_prefix: str, num_threads: int = 4) -> None:
    """
    从BAM文件提取polyA长度（按染色体并行处理）
    
    Args:
        bamfile: BAM文件路径
        output_prefix: 输出文件前缀
        num_threads: 用户指定的线程数
    """
    print(f"开始处理BAM文件: {bamfile}")
    
    # 处理线程数
    max_threads = len(AUTOSOMES)  # 22
    if num_threads >= max_threads:
        print(f"提示: 用户指定了{num_threads}个线程，但只有{max_threads}条常染色体需要处理")
        print(f"      实际将使用{max_threads}个线程（每条染色体一个线程）")
        actual_threads = max_threads
    else:
        print(f"使用{num_threads}个线程处理{max_threads}条染色体")
        actual_threads = num_threads
    
    # 检查BAM文件中存在哪些染色体
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        available_chroms = set(bam.references)
    
    chroms_to_process = [chrom for chrom in AUTOSOMES if chrom in available_chroms]
    missing_chroms = [chrom for chrom in AUTOSOMES if chrom not in available_chroms]
    
    if missing_chroms:
        print(f"警告: 以下染色体在BAM文件中不存在: {', '.join(missing_chroms)}")
    
    print(f"将处理{len(chroms_to_process)}条染色体: {', '.join(chroms_to_process)}")
    
    # 使用线程池并行处理各染色体
    all_results = []
    completed_count = 0
    
    with ThreadPoolExecutor(max_workers=actual_threads) as executor:
        # 提交所有染色体任务
        future_to_chrom = {
            executor.submit(process_chromosome, bamfile, chrom): chrom 
            for chrom in chroms_to_process
        }
        
        # 收集结果
        for future in as_completed(future_to_chrom):
            chrom = future_to_chrom[future]
            completed_count += 1
            
            try:
                result = future.result()
                read_count = len(result)
                all_results.extend(result)
                print(f"[{completed_count}/{len(chroms_to_process)}] {chrom} 完成，找到 {read_count} 条reads")
            except Exception as e:
                print(f"[{completed_count}/{len(chroms_to_process)}] {chrom} 处理出错: {e}")
    
    # 保存结果
    if all_results:
        df = pd.DataFrame(all_results, columns=['readID', 'polyA_length'])
        output_file = f"{output_prefix}_polyAlen_result.csv"
        df.to_csv(output_file, index=False)
        
        print(f"\n处理完成！")
        print(f"结果保存到: {output_file}")
        print(f"总共找到 {len(all_results)} 条带有polyA长度信息的reads")
    else:
        print("\n警告: 在BAM文件中未找到任何带有pt标签的reads")


def main():
    parser = argparse.ArgumentParser(description="从BAM文件提取polyA长度（按染色体并行处理）")
    parser.add_argument("-b", "--bamfile", type=str, required=True, 
                       help="BAM文件路径（需要有索引文件.bai）")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, 
                       help="输出文件前缀")
    parser.add_argument("-t", "--threads", type=int, default=4, 
                       help="线程数 (默认: 4，最大有效值: 22)")
    
    args = parser.parse_args()
    
    try:
        get_polyA_len(args.bamfile, args.output_prefix, args.threads)
    except Exception as e:
        print(f"错误: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())