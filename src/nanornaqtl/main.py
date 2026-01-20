import argparse
import sys
import os
import subprocess
import threading
import logging
from datetime import datetime

def setup_logging(output_dir, prefix):
    """
    设置日志系统，将子脚本输出写入日志文件
    """
    log_file = os.path.join(output_dir, f"{prefix}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    # 创建logger
    logger = logging.getLogger('nanornaqtl')
    logger.setLevel(logging.DEBUG)
    
    # 文件handler - 记录所有信息
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    
    # 清除已有的handlers
    logger.handlers.clear()
    logger.addHandler(file_handler)
    
    return logger, log_file


def run_command(command_list, description=None, working_dir=None, logger=None):
    """
    执行外部shell命令，捕获输出到日志文件
    - command_list: 命令及其参数的列表
    - description: 命令的描述
    - working_dir: 命令的执行目录
    - logger: 日志记录器
    """
    if description:
        print(f"执行命令: {description}")
    
    if logger:
        logger.info(f"执行命令: {' '.join(command_list)}")
    
    try:
        result = subprocess.run(
            command_list,
            check=True,
            text=True,
            capture_output=True,
            cwd=working_dir,
        )
        
        # 将子脚本的stdout和stderr写入日志
        if logger:
            if result.stdout:
                logger.info(f"[STDOUT]\n{result.stdout}")
            if result.stderr:
                logger.warning(f"[STDERR]\n{result.stderr}")
        
        print("命令执行成功。\n")

    except FileNotFoundError:
        error_msg = f"错误: 命令 '{command_list[0]}' 未找到。请确保所需软件已安装并在系统的PATH中。"
        print(error_msg, file=sys.stderr)
        if logger:
            logger.error(error_msg)
        sys.exit(1)
        
    except subprocess.CalledProcessError as e:
        error_msg = f"错误：命令执行失败！\n  命令: {' '.join(e.cmd)}\n  返回码: {e.returncode}"
        print(error_msg, file=sys.stderr)
        
        if logger:
            logger.error(error_msg)
            if e.stdout:
                logger.error(f"[STDOUT]\n{e.stdout}")
            if e.stderr:
                logger.error(f"[STDERR]\n{e.stderr}")
        
        # 打印简要错误信息到控制台
        print(f"  详细错误信息已写入日志文件", file=sys.stderr)
        if e.stderr:
            # 只打印stderr的最后几行到控制台
            stderr_lines = e.stderr.strip().split('\n')
            print(f"  错误摘要: {stderr_lines[-1] if stderr_lines else 'Unknown error'}", file=sys.stderr)
        
        sys.exit(1)


def run_prep_workflow(args):
    """
    prep 功能的核心工作流。
    根据用户输入，依次执行 getMapRead.py 和 samtools 命令，返回fastq、bam0、bam16、bam.bai。
    """
    print("==============================================")
    print(f"  nanornaqtl: 'prep' 工作流启动")
    print("==============================================")

    # 设置输出目录（提前，因为日志需要）
    output_dir_path = os.path.abspath(args.output_dir)
    os.makedirs(output_dir_path, exist_ok=True)
    
    # 设置日志系统
    logger, log_file = setup_logging(output_dir_path, args.prefix)
    print(f"日志文件: {log_file}")
    logger.info(f"开始 prep 工作流")
    logger.info(f"参数: {vars(args)}")

    # --- 1. 准备路径和文件名 ---

    # 获取用户输入的参数，并将路径转换为绝对路径
    input_bam_path = os.path.abspath(args.bam)
    prefix = args.prefix
    threads = str(args.threads)

    # print("--- 步骤 1: 预排序输入的BAM文件 ---")
    tmp_sorted_bam_name = f"{prefix}_tmp_sorted.bam"
    tmp_sorted_bam_path = os.path.join(output_dir_path, tmp_sorted_bam_name)

    cmd_sort = [
        "samtools",
        "sort",
        "-@",
        threads,
        "-o",
        tmp_sorted_bam_path,
        input_bam_path,
    ]
    run_command(
        cmd_sort, description="预排序输入的BAM文件...", working_dir=output_dir_path, logger=logger
    )

    # 为排序后的临时BAM建立索引，以供pysam使用
    run_command(
        ["samtools", "index", tmp_sorted_bam_path],
        description="为排序后的临时BAM建立索引...",
        working_dir=output_dir_path,
        logger=logger
    )

    # 构建所有需要用到的文件名
    map0_bam = f"{prefix}_calls_sorted_map0.bam"
    map16_bam = f"{prefix}_calls_sorted_map16.bam"
    map_bam = f"{prefix}_calls_sorted_map.bam"
    output_fastq = f"{prefix}.fastq"

    # 找到 getMapRead.py 脚本的绝对路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    getMapRead_script_path = os.path.join(base_dir, "scripts", "getMapRead.py")
    if not os.path.exists(getMapRead_script_path):
        print(f"错误：无法找到脚本 {getMapRead_script_path}", file=sys.stderr)
        print("请确认 'scripts/getMapRead.py' 文件在包内。", file=sys.stderr)
        sys.exit(1)

    print("准备工作完成:")
    print(f"  - 输入BAM: {input_bam_path}")
    print(f"  - 输出目录: {output_dir_path}")
    print(f"  - 输出文件前缀: {prefix}")

    # --- 2. 执行命令流程 ---

    # 步骤 1: 运行 getMapRead.py
    cmd_getMapRead = [
        "python",
        getMapRead_script_path,
        "-b",
        tmp_sorted_bam_path,
        "-p",
        prefix,
        "-t",
        threads,
        "-q",
        str(args.min_mapq),
    ]
    # 在指定的输出目录中执行此命令
    run_command(
        cmd_getMapRead,
        description="提取比对成功的reads...",
        working_dir=output_dir_path,
    )

    # 步骤 2: 为 getMapRead.py 生成的BAM文件建立索引
    run_command(
        ["samtools", "index", map0_bam],
        description="为正链BAM建立索引...",
        working_dir=output_dir_path,
        logger=logger
    )
    run_command(
        ["samtools", "index", map16_bam],
        description="为反链BAM建立索引...",
        working_dir=output_dir_path,
        logger=logger
    )
    run_command(
        ["samtools", "index", map_bam],
        description="为最终BAM建立索引...",
        working_dir=output_dir_path,
        logger=logger
    )

    # 步骤 3: 从最终的比对BAM中提取FASTQ
    cmd_fastq = ["samtools", "fastq", "-@", threads, "-0", output_fastq, map_bam]
    run_command(
        cmd_fastq,
        description="从最终的比对BAM中提取FASTQ...",
        working_dir=output_dir_path,
        logger=logger
    )

    # --- 3. 清理临时文件 ---
    print(f"\n--- 正在清理临时文件及其索引 ---")
    try:
        os.remove(tmp_sorted_bam_path)
        os.remove(f"{tmp_sorted_bam_path}.bai")  # 清理索引文件
        print("临时文件清理完毕。")
    except OSError as e:
        print(
            f"警告：无法删除临时文件 {tmp_sorted_bam_path} 或其索引。错误: {e}",
            file=sys.stderr,
        )

    print("\n==============================================")
    print("  'prep' 工作流全部执行完毕！")
    print(f"  所有输出文件已生成在: {output_dir_path}")
    print(f"  日志文件: {log_file}")
    print("==============================================")


def validate_phenotype_parameters(args):
    """
    验证每种分子表型分析所需的参数是否完整
    现在大部分验证已经通过argparse的required参数处理，这里只做额外的逻辑检查
    """
    # 由于现在使用子解析器，必需的参数已经通过argparse自动验证
    # 这里可以添加一些额外的逻辑验证，比如文件是否存在等
    pass


def run_pheno_workflow(args):
    """
    pheno 功能的核心工作流。
    """
    print("==============================================")
    print(f"  nanornaqtl: 'pheno' 工作流启动")
    print(f"  分析分子表型: {args.phenotype}")
    print("==============================================")

    # 支持的分子表型列表
    allowed_phenotypes = {
        "m6A", "m5C", "pseU", "inosine",
        "polyA_tail", "APA", "intron_retention",
    }

    if args.phenotype not in allowed_phenotypes:
        print(
            f"错误：不支持的分子表型 '{args.phenotype}'，仅支持: {', '.join(allowed_phenotypes)}",
            file=sys.stderr,
        )
        sys.exit(1)

    # 验证 motifPaint 和 motif 参数的组合
    if hasattr(args, 'motifPaint') and hasattr(args, 'motif'):
        if args.motifPaint and not args.motif:
            print("⚠️  警告: --motifPaint 需要配合 --motif 使用")
            print("   当前未指定 --motif，将按 base 模式运行，跳过 motifPaint")
            args.motifPaint = False

    # 验证参数
    validate_phenotype_parameters(args)

    # 获取脚本路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    scripts_dir = os.path.join(base_dir, "scripts")

    # 设置输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 设置日志系统
    logger, log_file = setup_logging(output_dir, args.prefix)
    print(f"日志文件: {log_file}")
    logger.info(f"开始 {args.phenotype} 分子表型分析")
    logger.info(f"参数: {vars(args)}")

    # 根据分子表型执行相应的分析流程
    if args.phenotype in ["m6A", "m5C", "pseU", "inosine"]:
        run_modification_analysis(args, scripts_dir, output_dir, logger)
    elif args.phenotype == "polyA_tail":
        run_polyA_analysis(args, scripts_dir, output_dir, logger)
    elif args.phenotype == "APA":
        run_APA_analysis(args, scripts_dir, output_dir, logger)
    elif args.phenotype == "intron_retention":
        run_intron_retention_analysis(args, scripts_dir, output_dir, logger)

    print("\n==============================================")
    print(f"  '{args.phenotype}' 分子表型分析全部执行完毕！")
    print(f"  所有输出文件已生成在: {output_dir}")
    print(f"  日志文件: {log_file}")
    print("==============================================")


def run_modification_analysis(args, scripts_dir, output_dir, logger):
    """
    运行修饰位点分析（m6A, m5C, pseU, inosine）
    需要3步：pileupRead.py → calculateModRate.py → mergeModSites.py
    """
    phenotype = args.phenotype

    # 从map.bam构建map0.bam和map16.bam的路径
    bam_path = os.path.abspath(args.bam)
    bam0_path = bam_path.replace(".bam", "0.bam")
    bam16_path = bam_path.replace(".bam", "16.bam")

    # 检查BAM文件是否存在
    if not os.path.exists(bam0_path):
        print(
            f"错误：找不到文件 {bam0_path}，请确保已运行prep步骤生成了map0.bam文件",
            file=sys.stderr,
        )
        sys.exit(1)
    if not os.path.exists(bam16_path):
        print(
            f"错误：找不到文件 {bam16_path}，请确保已运行prep步骤生成了map16.bam文件",
            file=sys.stderr,
        )
        sys.exit(1)

    # 第一步：运行pileupRead.py分别处理正负链
    print(f"--- 1. 获取{phenotype}修饰位点 ---")

    # 正链分析 (strand +)
    cmd_pos = [
        "python",
        os.path.join(scripts_dir, "pileupRead.py"),
        "-o", args.prefix,
        "-b", bam0_path,
        "-s", "+",
        "-m", phenotype,
        "-f", str(args.mod_threshold),
        "-q", str(args.min_qscore),
        "-t", str(args.threads),
        "--min_mapq", str(args.min_mapq),
    ]
    if args.motif:
        cmd_pos.append("--motif")

    run_command(
        cmd_pos, description=f"获取{phenotype}修饰位点（正链）", working_dir=output_dir, logger=logger
    )

    # 负链分析 (strand -)
    cmd_neg = [
        "python",
        os.path.join(scripts_dir, "pileupRead.py"),
        "-o", args.prefix,
        "-b", bam16_path,
        "-s", "-",
        "-m", phenotype,
        "-f", str(args.mod_threshold),
        "-q", str(args.min_qscore),
        "-t", str(args.threads),
        "--min_mapq", str(args.min_mapq),
    ]
    if args.motif:
        cmd_neg.append("--motif")

    run_command(
        cmd_neg, description=f"获取{phenotype}修饰位点（负链）", working_dir=output_dir, logger=logger
    )

    # 第二步：运行calculateModRate.py计算修饰率
    print(f"--- 2. 计算{phenotype}修饰率 ---")

    # 构建第一步输出文件名
    pos_file = (
        f"{args.prefix}_f_{phenotype}_motif_tmp.txt"
        if args.motif
        else f"{args.prefix}_f_{phenotype}_base_tmp.txt"
    )
    neg_file = (
        f"{args.prefix}_r_{phenotype}_motif_tmp.txt"
        if args.motif
        else f"{args.prefix}_r_{phenotype}_base_tmp.txt"
    )

    # 正链计算修饰率
    cmd_calc_pos = [
        "python", os.path.join(scripts_dir, "calculateModRate.py"),
        "-b", bam0_path,
        "-f", pos_file,
        "-o", args.prefix,
        "-s", "+",
        "-m", phenotype,
        "-r", str(args.min_rate),
        "-c", str(args.min_cov),
        "-t", str(args.threads),
        "--min_mapq", str(args.min_mapq),
    ]

    run_command(
        cmd_calc_pos,
        description=f"计算{phenotype}修饰率（正链）",
        working_dir=output_dir,
        logger=logger
    )

    # 负链计算修饰率
    cmd_calc_neg = [
        "python", os.path.join(scripts_dir, "calculateModRate.py"),
        "-b", bam16_path,
        "-f", neg_file,
        "-o", args.prefix,
        "-s", "-",
        "-m", phenotype,
        "-r", str(args.min_rate),
        "-c", str(args.min_cov),
        "-t", str(args.threads),
        "--min_mapq", str(args.min_mapq),
    ]

    run_command(
        cmd_calc_neg,
        description=f"计算{phenotype}修饰率（负链）",
        working_dir=output_dir,
        logger=logger
    )

    # 第三步：合并正负链结果并绘制motif图
    print(f"--- 第三步：合并{phenotype}结果 ---")

    cmd_merge = [
        "python",
        os.path.join(scripts_dir, "mergeModSites.py"),
        "-o",
        args.prefix,
        "-m",
        phenotype,
        "-f",
        os.path.abspath(args.fasta),
    ]

    if args.motifPaint and args.motif:
        cmd_merge.append("--motifPaint")
    elif args.motifPaint and not args.motif:
        print("⚠️  警告: --motifPaint 需要配合 --motif 使用")
        print("   当前未指定 --motif，将按 base 模式运行，跳过 motifPaint")

    if args.metaPlotR:
        cmd_merge.append("--metaPlotR")

    run_command(
        cmd_merge,
        description=f"合并{phenotype}结果并生成motif图",
        working_dir=output_dir,
        logger=logger
    )


def run_polyA_analysis(args, scripts_dir, output_dir, logger):
    """运行polyA长度分析"""
    print("--- 运行polyA长度分析 ---")

    bam_path = os.path.abspath(args.bam)

    cmd = [
        "python",
        os.path.join(scripts_dir, "getPolyAlen.py"),
        "-b",
        bam_path,
        "-o",
        args.prefix,
        "-t",
        str(args.threads),
    ]

    run_command(cmd, description="分析polyA长度", working_dir=output_dir, logger=logger)


def run_APA_analysis(args, scripts_dir, output_dir, logger):
    """运行APA分析"""
    print("--- 运行APA分析 ---")

    bam_path = os.path.abspath(args.bam)

    cmd = [
        "python",
        os.path.join(scripts_dir, "getAPA.py"),
        "-b",
        bam_path,
        "-o",
        args.prefix,
        "-a",
        os.path.abspath(args.apadb),
    ]

    run_command(cmd, description="分析APA", working_dir=output_dir, logger=logger)


def run_intron_retention_analysis(args, scripts_dir, output_dir, logger):
    """运行内含子滞留率分析"""
    print("--- 运行内含子滞留率分析 ---")

    bam_path = os.path.abspath(args.bam)

    cmd = [
        "python",
        os.path.join(scripts_dir, "getIntronRetention.py"),
        "-g",
        os.path.abspath(args.gtf),
        "-b",
        bam_path,
        "-o",
        output_dir,
        "-p",
        args.prefix,
    ]

    run_command(cmd, description="分析内含子滞留率", working_dir=output_dir, logger=logger)


def run_qtl_workflow(args):
    """
    QTL功能的核心工作流。
    根据用户指定的分子表型运行相应的QTL分析流程。
    """
    print("==============================================")
    print(f"  nanornaqtl: 'qtl' 工作流启动")
    print(f"  分析分子表型: {args.molecular_phenotype}")
    print("==============================================")

    # 根据分子表型执行相应的QTL分析
    if args.molecular_phenotype in ["m6A", "m5C", "pseU", "inosine"]:
        run_modification_qtl_analysis(args)
    elif args.molecular_phenotype == "polyA_tail":
        run_polya_qtl_analysis(args)
    elif args.molecular_phenotype == "APA":
        run_apa_qtl_analysis(args)
    elif args.molecular_phenotype == "isoform":
        run_isoform_qtl_analysis(args)
    elif args.molecular_phenotype == "intron_retention":
        run_intron_retention_qtl_analysis(args)

    print("\n==============================================")
    print(f"  '{args.molecular_phenotype}' QTL分析全部执行完毕！")
    print(f"  结果文件已生成")
    print("==============================================")


def run_modification_qtl_analysis(args):
    """
    运行修饰QTL分析（m6A, m5C, pseU, inosine）
    使用parallelModqtlProcessor.py的功能
    """
    print(f"--- 运行{args.molecular_phenotype}修饰QTL分析 ---")

    # 获取脚本路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    parallel_script = os.path.join(base_dir, "scripts", "parallelModqtlProcessor.py")

    if not os.path.exists(parallel_script):
        print(f"错误：找不到脚本 {parallel_script}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 构建命令
    cmd = [
        "python",
        parallel_script,
        "-b",
        os.path.abspath(args.bam),
        "--snp_info",
        os.path.abspath(args.snp_info),
        "-o",
        args.prefix,
        "-csv",
        os.path.abspath(args.modification_csv),
        "-pkl",
        os.path.abspath(args.read_mod_dict),
        "--geno_size",
        os.path.abspath(args.geno_size),
        "-m",
        args.molecular_phenotype,
        "-q",
        str(args.min_qscore),
        "--threads",
        str(args.threads),
        "-c",
        str(args.min_coverage),
        "--mcmc_samples",
        str(args.mcmc_samples)
    ]

    if args.keep_tmp:
        cmd.append("--keep_tmp")

    run_command(
        cmd,
        description=f"运行{args.molecular_phenotype}修饰QTL分析",
        working_dir=output_dir,
    )


def run_polya_qtl_analysis(args):
    """
    运行polyA尾长QTL分析
    使用parallelPlQTLProcessor.py的功能
    """
    print("--- 运行polyA尾长QTL分析 ---")

    # 获取脚本路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    parallel_script = os.path.join(base_dir, "scripts", "parallelPlQTLProcessor.py")

    if not os.path.exists(parallel_script):
        print(f"错误：找不到脚本 {parallel_script}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 构建命令
    cmd = [
        "python",
        parallel_script,
        "-b",
        os.path.abspath(args.bam),
        "--snp_info",
        os.path.abspath(args.snp_info),
        "-o",
        args.prefix,
        "-csv",
        os.path.abspath(args.polya_csv),
        "--geno_size",
        os.path.abspath(args.geno_size),
        "-q",
        str(args.min_qscore),
        "--min_coverage",
        str(args.min_coverage),
        "--mcmc_samples",
        str(args.mcmc_samples),
        "--threads",
        str(args.threads),
    ]
    if args.keep_tmp:
        cmd.append("--keep_tmp")

    run_command(cmd, description="运行polyA尾长QTL分析", working_dir=output_dir)


def run_apa_qtl_analysis(args):
    """
    运行APA QTL分析
    使用parallel_iu_3au_qtlProcessor.py的功能
    """
    print("——— 运行APA QTL分析 ———")

    # 获取脚本路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    parallel_script = os.path.join(
        base_dir, "scripts", "parallel_iu_3au_qtlProcessor.py"
    )

    if not os.path.exists(parallel_script):
        print(f"错误：找不到脚本 {parallel_script}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 构建命令
    cmd = [
        "python",
        parallel_script,
        "-b",
        os.path.abspath(args.bam),
        "--snp_info",
        os.path.abspath(args.snp_info),
        "-o",
        args.prefix,
        "-f",
        os.path.abspath(args.read_overlap_file),
        "--geno_size",
        os.path.abspath(args.geno_size),
        "-m",
        "APA",
        "-q",
        str(args.min_qscore),
        "-t",
        str(args.threads),
        "--min_coverage",
        str(args.min_coverage),
        "--mcmc_samples",
        str(args.mcmc_samples)
    ]

    if args.keep_tmp:
        cmd.append("--keep_tmp")

    run_command(cmd, description="运行APA QTL分析", working_dir=output_dir)


def run_isoform_qtl_analysis(args):
    """
    运行isoform QTL分析
    使用parallel_iu_3au_qtlProcessor.py的功能
    """
    print("——— 运行isoform QTL分析 ———")

    # 获取脚本路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    parallel_script = os.path.join(
        base_dir, "scripts", "parallel_iu_3au_qtlProcessor.py"
    )

    if not os.path.exists(parallel_script):
        print(f"错误：找不到脚本 {parallel_script}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 构建命令
    cmd = [
        "python",
        parallel_script,
        "-b",
        os.path.abspath(args.bam),
        "--snp_info",
        os.path.abspath(args.snp_info),
        "-o",
        args.prefix,
        "-f",
        os.path.abspath(args.read_overlap_file),
        "--geno_size",
        os.path.abspath(args.geno_size),
        "-m",
        "isoform",
        "-q",
        str(args.min_qscore),
        "-t",
        str(args.threads),
        "--min_coverage",
        str(args.min_coverage),
        "--mcmc_samples",
        str(args.mcmc_samples)
    ]

    if args.keep_tmp:
        cmd.append("--keep_tmp")

    run_command(cmd, description="运行isoform QTL分析", working_dir=output_dir)


def run_intron_retention_qtl_analysis(args):
    """
    运行内含子滞留QTL分析
    使用parallelirqtlProcessor.py的功能
    """
    print("——— 运行内含子滞留QTL分析 ———")

    # 获取脚本路径
    base_dir = os.path.dirname(os.path.realpath(__file__))
    parallel_script = os.path.join(base_dir, "scripts", "parallelirqtlProcessor.py")

    if not os.path.exists(parallel_script):
        print(f"错误：找不到脚本 {parallel_script}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 构建命令
    cmd = [
        "python",
        parallel_script,
        "-b",
        os.path.abspath(args.bam),
        "--snp_info",
        os.path.abspath(args.snp_info),
        "-o",
        args.prefix,
        "-csv",
        os.path.abspath(args.ir_csv),
        "--geno_size",
        os.path.abspath(args.geno_size),
        "-q",
        str(args.min_qscore),
        "--min_coverage",
        str(args.min_coverage),
        "--mcmc_samples",
        str(args.mcmc_samples),
        "--threads",
        str(args.threads),
    ]

    if args.keep_tmp:
        cmd.append("--keep_tmp")

    run_command(cmd, description="运行内含子滞留QTL分析", working_dir=output_dir)


def main():
    """程序的主入口，负责定义和解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="nanornaqtl: 一个用于分析Nanopore Direct RNA数据分析多种分子表型及xQTLs的工具。",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="可用的子命令")

    # ---- 定义 "prep" 子命令 ----
    parser_prep = subparsers.add_parser(
        "prep",
        help="数据预处理：利用dorado生成的BAM文件，生成后续分析所需的BAM和FASTQ文件。",
    )
    # 为 'prep' 添加它需要的所有参数
    parser_prep.add_argument(
        "-b", "--bam", required=True, help="输入的Dorado BAM文件路径。"
    )
    parser_prep.add_argument(
        "-p", "--prefix", required=True, help="所有输出文件的统一前缀名。"
    )
    parser_prep.add_argument(
        "-t", "--threads", type=int, default=4, help="使用的线程数 (默认: 4)。"
    )
    parser_prep.add_argument(
        "-o",
        "--output_dir",
        default=".",
        help="所有输出文件存放的目录，如果目录不存在则会新建 (默认: 当前目录)。",
    )
    parser_prep.add_argument(
        "-q",
        "--min_mapq",
        type=int,
        default=0,
        help="最小比对质量分数 (默认: 0)。",
    )
    # 将 prep 子命令与它的执行函数 `run_prep_workflow` 绑定起来
    parser_prep.set_defaults(func=run_prep_workflow)

    # ---- 定义 "pheno" 子命令 ----
    parser_pheno = subparsers.add_parser(
        "pheno",
        help="分子表型识别：选择具体的分子表型进行分析",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # 为pheno创建子子解析器
    pheno_subparsers = parser_pheno.add_subparsers(
        dest="phenotype",
        help="""
使用方法: nanornaqtl pheno <分子表型> [参数]
例如: nanornaqtl pheno m6A -b input.bam -p sample -o output/
请选择要分析的分子表型：

""",
        required=True,
    )

    # 定义通用参数添加函数
    def add_common_args(parser):
        """为子解析器添加通用参数"""
        parser.add_argument(
            "-b", "--bam", required=True, help="从prep步骤得到的map.bam文件路径"
        )
        parser.add_argument("-p", "--prefix", required=True, help="输出文件前缀名")
        parser.add_argument(
            "-o",
            "--output-dir",
            required=True,
            help="输出目录路径，如果目录不存在则会新建",
        )
        parser.add_argument(
            "-t",
            "--threads",
            type=int,
            default=4,
            help="使用的线程数，不建议设置很高（默认: 4）",
        )
        parser.add_argument(
            "--min_mapq",
            type=int,
            default=0,
            help="最小read比对质量(MAPQ)，用于统一筛选所有reads（默认: 0）"
        )

    # === m6A修饰位点分析 ===
    parser_m6a = pheno_subparsers.add_parser(
        "m6A",
        help="N6-甲基腺苷修饰位点识别",
        description="nanornaqtl pheno m6A: 进行m6A修饰位点的检测和分析",
    )
    add_common_args(parser_m6a)
    parser_m6a.add_argument(
        "-f", "--fasta", required=True, help="参考基因组fasta文件路径"
    )
    parser_m6a.add_argument(
        "--mod_threshold", type=float, default=0.75, help="修饰阈值（默认: 0.75）"
    )
    parser_m6a.add_argument(
        "--min_qscore",
        type=int,
        default=10,
        help="修饰位点最小碱基质量分数（默认: 10）",
    )
    parser_m6a.add_argument(
        "--min_rate", type=float, default=0.1, help="修饰位点最小修饰率（默认: 0.1）"
    )
    parser_m6a.add_argument(
        "--min_cov", type=int, default=5, help="修饰位点最小覆盖度（默认: 5）"
    )
    parser_m6a.add_argument(
        "--motif", action="store_true", help="可选，是否对修饰位点进行motif筛选"
    )
    parser_m6a.add_argument(
        "--motifPaint", action="store_true", help="可选，是否绘制motif图"
    )
    parser_m6a.add_argument(
        "--metaPlotR", action="store_true", help="可选，是否生成MetaPlotR所需的bed文件"
    )
    parser_m6a.set_defaults(func=run_pheno_workflow)

    # === m5C修饰位点分析 ===
    parser_m5c = pheno_subparsers.add_parser(
        "m5C",
        help="5-甲基胞嘧啶修饰位点识别",
        description="nanornaqtl pheno m5C: 进行m5C修饰位点的检测和分析",
    )
    add_common_args(parser_m5c)
    parser_m5c.add_argument(
        "-f", "--fasta", required=True, help="参考基因组fasta文件路径"
    )
    parser_m5c.add_argument(
        "--mod_threshold", type=float, default=0.75, help="修饰阈值（默认: 0.75）"
    )
    parser_m5c.add_argument(
        "--min_qscore",
        type=int,
        default=10,
        help="修饰位点最小碱基质量分数（默认: 10）",
    )
    parser_m5c.add_argument(
        "--min_rate", type=float, default=0.1, help="修饰位点最小修饰率（默认: 0.1）"
    )
    parser_m5c.add_argument(
        "--min_cov", type=int, default=5, help="修饰位点最小覆盖度（默认: 5）"
    )
    parser_m5c.add_argument(
        "--motif", action="store_true", help="可选，是否对修饰位点进行motif筛选"
    )
    parser_m5c.add_argument(
        "--motifPaint", action="store_true", help="可选，是否绘制motif图"
    )
    parser_m5c.add_argument(
        "--metaPlotR", action="store_true", help="可选，是否生成MetaPlotR所需的bed文件"
    )
    parser_m5c.set_defaults(func=run_pheno_workflow)

    # === pseU修饰位点分析 ===
    parser_pseu = pheno_subparsers.add_parser(
        "pseU",
        help="假尿嘧啶修饰位点识别",
        description="nanornaqtl pheno pseU: 进行pseU修饰位点的检测和分析",
    )
    add_common_args(parser_pseu)
    parser_pseu.add_argument(
        "-f", "--fasta", required=True, help="参考基因组fasta文件路径"
    )
    parser_pseu.add_argument(
        "--mod_threshold", type=float, default=0.75, help="修饰阈值（默认: 0.75）"
    )
    parser_pseu.add_argument(
        "--min_qscore",
        type=int,
        default=10,
        help="修饰位点最小碱基质量分数（默认: 10）",
    )
    parser_pseu.add_argument(
        "--min_rate", type=float, default=0.1, help="修饰位点最小修饰率（默认: 0.1）"
    )
    parser_pseu.add_argument(
        "--min_cov", type=int, default=5, help="修饰位点最小覆盖度（默认: 5）"
    )
    parser_pseu.add_argument(
        "--motif", action="store_true", help="可选，是否对修饰位点进行motif筛选"
    )
    parser_pseu.add_argument(
        "--motifPaint", action="store_true", help="可选，是否绘制motif图"
    )
    parser_pseu.add_argument(
        "--metaPlotR", action="store_true", help="可选，是否生成MetaPlotR所需的bed文件"
    )
    parser_pseu.set_defaults(func=run_pheno_workflow)

    # === inosine修饰位点分析 ===
    parser_inosine = pheno_subparsers.add_parser(
        "inosine",
        help="肌苷修饰位点识别",
        description="nanornaqtl pheno inosine: 进行inosine修饰位点的检测和分析",
    )
    add_common_args(parser_inosine)
    parser_inosine.add_argument(
        "-f", "--fasta", required=True, help="参考基因组fasta文件路径"
    )
    parser_inosine.add_argument(
        "--mod_threshold", type=float, default=0.75, help="修饰阈值（默认: 0.75）"
    )
    parser_inosine.add_argument(
        "--min_qscore",
        type=int,
        default=10,
        help="修饰位点最小碱基质量分数（默认: 10）",
    )
    parser_inosine.add_argument(
        "--min_rate", type=float, default=0.1, help="修饰位点最小修饰率（默认: 0.1）"
    )
    parser_inosine.add_argument(
        "--min_cov", type=int, default=5, help="修饰位点最小覆盖度（默认: 5）"
    )
    parser_inosine.add_argument(
        "--motif", action="store_true", help="可选，是否对修饰位点进行motif筛选"
    )
    parser_inosine.add_argument(
        "--motifPaint", action="store_true", help="可选，是否绘制motif图"
    )
    parser_inosine.add_argument(
        "--metaPlotR", action="store_true", help="可选，是否生成MetaPlotR所需的bed文件"
    )
    parser_inosine.set_defaults(func=run_pheno_workflow)

    # === polyA尾长度分析 ===
    parser_polya = pheno_subparsers.add_parser(
        "polyA_tail",
        help="polyA尾长度分析",
        description="nanornaqtl pheno polyA_tail: 分析转录本的polyA尾长度分布",
    )
    add_common_args(parser_polya)
    parser_polya.set_defaults(func=run_pheno_workflow)

    # === APA分析 ===
    parser_apa = pheno_subparsers.add_parser(
        "APA",
        help="多聚腺苷酸化位点分析",
        description="nanornaqtl pheno APA: 分析选择性多聚腺苷酸化位点",
    )
    add_common_args(parser_apa)
    parser_apa.add_argument("--apadb", required=True, help="APAdb bed文件路径")
    parser_apa.set_defaults(func=run_pheno_workflow)

    # === 内含子滞留率分析 ===
    parser_ir = pheno_subparsers.add_parser(
        "intron_retention",
        help="内含子滞留率分析",
        description="nanornaqtl pheno intron_retention: 分析转录本中内含子的滞留程度",
    )
    add_common_args(parser_ir)
    parser_ir.add_argument("-g", "--gtf", required=True, help="基因注释GTF文件路径")
    parser_ir.set_defaults(func=run_pheno_workflow)

    # ---- 定义 "qtl" 子命令 ----
    parser_qtl = subparsers.add_parser(
        "qtl",
        help="QTL识别：支持8种分子表型的QTL分析",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # 为qtl创建子子解析器
    qtl_subparsers = parser_qtl.add_subparsers(
        dest="molecular_phenotype",
        help="""
使用方法: nanornaqtl qtl <分子表型> [参数]
例如: nanornaqtl qtl m6A -b input.bam -p sample --snp_info snps.vcf
请选择要分析的分子表型：

""",
        required=True,
    )

    # 定义通用QTL参数添加函数
    def add_common_qtl_args(parser):
        """为QTL子解析器添加通用参数"""
        parser.add_argument(
            "-b", "--bam", required=True, help="从prep步骤得到的map.bam文件路径"
        )
        parser.add_argument("-p", "--prefix", required=True, help="输出文件前缀名")
        parser.add_argument(
            "-o",
            "--output_dir",
            required=True,
            help="输出目录路径，如果目录不存在则会新建",
        )
        parser.add_argument(
            "--snp_info", required=True, help="SNP信息文件路径（VCF格式或其他格式）"
        )
        parser.add_argument(
            "--geno_size",
            required=True,
            help="基因组大小文件路径（如hg19.chrom.sizes）",
        )
        parser.add_argument(
            "-c",
            "--min_coverage",
            type=int,
            default=8,
            help="最小总覆盖度（默认: 8）",
        )
        parser.add_argument(
            "--mcmc_samples",
            type=int,
            default=2000,
            help="MCMC采样数(default=2000)"
        )
        parser.add_argument(
            "-q",
            "--min_qscore",
            type=int,
            default=10,
            help="最小碱基质量分数（默认: 10）",
        )
        parser.add_argument(
            "-t", "--threads", type=int, default=4, help="使用的线程数（默认: 4）"
        )
        parser.add_argument("--keep_tmp", action="store_true", help="保留临时文件")

    # === m6A修饰QTL分析 ===
    parser_m6a_qtl = qtl_subparsers.add_parser(
        "m6A",
        help="N6-甲基腺苷修饰QTL识别",
        description="nanornaqtl qtl m6A: 进行m6A修饰位点的QTL分析",
    )
    add_common_qtl_args(parser_m6a_qtl)
    parser_m6a_qtl.add_argument(
        "--modification_csv",
        required=True,
        help="m6A修饰位点结果CSV文件路径（来自pheno步骤）",
    )
    parser_m6a_qtl.add_argument(
        "--read_mod_dict",
        required=True,
        help="修饰read字典PKL文件路径（来自pheno步骤）",
    )
    parser_m6a_qtl.set_defaults(func=run_qtl_workflow)

    # === m5C修饰QTL分析 ===
    parser_m5c_qtl = qtl_subparsers.add_parser(
        "m5C",
        help="5-甲基胞嘧啶修饰QTL识别",
        description="nanornaqtl qtl m5C: 进行m5C修饰位点的QTL分析",
    )
    add_common_qtl_args(parser_m5c_qtl)
    parser_m5c_qtl.add_argument(
        "--modification_csv",
        required=True,
        help="m5C修饰位点结果CSV文件路径（来自pheno步骤）",
    )
    parser_m5c_qtl.add_argument(
        "--read_mod_dict",
        required=True,
        help="修饰read字典PKL文件路径（来自pheno步骤）",
    )
    parser_m5c_qtl.set_defaults(func=run_qtl_workflow)

    # === pseU修饰QTL分析 ===
    parser_pseu_qtl = qtl_subparsers.add_parser(
        "pseU",
        help="假尿嘧啶修饰QTL识别",
        description="nanornaqtl qtl pseU: 进行pseU修饰位点的QTL分析",
    )
    add_common_qtl_args(parser_pseu_qtl)
    parser_pseu_qtl.add_argument(
        "--modification_csv",
        required=True,
        help="pseU修饰位点结果CSV文件路径（来自pheno步骤）",
    )
    parser_pseu_qtl.add_argument(
        "--read_mod_dict",
        required=True,
        help="修饰read字典PKL文件路径（来自pheno步骤）",
    )
    parser_pseu_qtl.set_defaults(func=run_qtl_workflow)

    # === inosine修饰QTL分析 ===
    parser_inosine_qtl = qtl_subparsers.add_parser(
        "inosine",
        help="肌苷修饰QTL识别",
        description="nanornaqtl qtl inosine: 进行inosine修饰位点的QTL分析",
    )
    add_common_qtl_args(parser_inosine_qtl)
    parser_inosine_qtl.add_argument(
        "--modification_csv",
        required=True,
        help="inosine修饰位点结果CSV文件路径（来自pheno步骤）",
    )
    parser_inosine_qtl.add_argument(
        "--read_mod_dict",
        required=True,
        help="修饰read字典PKL文件路径（来自pheno步骤）",
    )
    parser_inosine_qtl.set_defaults(func=run_qtl_workflow)

    # === polyA尾长QTL分析 ===
    parser_polya_qtl = qtl_subparsers.add_parser(
        "polyA_tail",
        help="polyA尾长QTL识别",
        description="nanornaqtl qtl polyA_tail: 进行polyA尾长的QTL分析",
    )
    add_common_qtl_args(parser_polya_qtl)
    parser_polya_qtl.add_argument(
        "--polya_csv", required=True, help="polyA尾长结果CSV文件路径（来自pheno步骤）"
    )
    parser_polya_qtl.set_defaults(func=run_qtl_workflow)

    # === APA QTL分析 ===
    parser_apa_qtl = qtl_subparsers.add_parser(
        "APA",
        help="多聚腺苷酸化位点QTL识别",
        description="nanornaqtl qtl APA: 进行选择性多聚腺苷酸化位点的QTL分析",
    )
    add_common_qtl_args(parser_apa_qtl)
    parser_apa_qtl.add_argument(
        "--read_overlap_file",
        required=True,
        help="read与APA位点的重叠文件路径（来自pheno步骤，*_APA_result.csv）",
    )
    parser_apa_qtl.set_defaults(func=run_qtl_workflow)

    # === isoform QTL分析 ===
    parser_isoform_qtl = qtl_subparsers.add_parser(
        "isoform",
        help="转录本异构体QTL识别",
        description="nanornaqtl qtl isoform: 进行转录本异构体使用的QTL分析",
    )
    add_common_qtl_args(parser_isoform_qtl)
    parser_isoform_qtl.add_argument(
        "--read_overlap_file",
        required=True,
        help="read与isoform的重叠文件路径（来自pheno步骤，OUT.transcript_model_reads.tsv.gz）",
    )
    parser_isoform_qtl.set_defaults(func=run_qtl_workflow)

    # === 内含子滞留QTL分析 ===
    parser_ir_qtl = qtl_subparsers.add_parser(
        "intron_retention",
        help="内含子滞留QTL识别",
        description="nanornaqtl qtl intron_retention: 进行内含子滞留率的QTL分析",
    )
    add_common_qtl_args(parser_ir_qtl)
    parser_ir_qtl.add_argument(
        "--ir_csv", required=True, help="内含子滞留结果CSV文件路径（来自pheno步骤）"
    )
    parser_ir_qtl.set_defaults(func=run_qtl_workflow)

    # 解析用户从命令行输入的参数
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # 根据用户选择的子命令，执行绑定的函数
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
