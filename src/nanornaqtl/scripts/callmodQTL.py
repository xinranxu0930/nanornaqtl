import pandas as pd
import pysam
import argparse
import pickle
import numpy as np
import pymc as pm
import arviz as az
from scipy import stats
import warnings
from scipy.stats import gaussian_kde
from statsmodels.stats.multitest import multipletests
from scipy.special import gammaln, betaln

warnings.filterwarnings("ignore")


def create_mod_dict(row):
    if row["mod_rate"] != 100:
        return row["pos_0base"], row["mod_rate"]


def count_haplotype(chrom, end, strand, bamfile, mod_dict, mod_id_dict, snp_file_dict, base_minQ):
    """
    统计haplotype信息
    
    参数:
        chrom: 染色体
        end: 染色体结束位置
        strand: 正负链
        bamfile: bam文件路径
        mod_dict: 修饰位点字典 {pos0: mod_rate}
        mod_id_dict: 修饰read字典 {chrom_pos0_strand: "read1;read2;..."}
        snp_file_dict: SNP信息字典 {pos0: "rsid;A1;A2;MAF"}
        base_minQ: 最小碱基质量
    
    返回:
        haplotype_df: 包含所有SNP-mod对的DataFrame
    """
    # 初始化三个字典
    mod_all_reads = {}   # {pos0: [read1, read2, read3, ...]} 所有覆盖该修饰位点的reads
    snp_bases = {}       # {pos0: {read1: "A", read2: "T", ...}} 非修饰位点的SNP碱基信息
    snp_mod_bases = {}   # {pos0: {read1: "A", read2: "T", ...}} 修饰位点的SNP碱基信息
    
    # 使用pysam遍历bam文件
    with pysam.AlignmentFile(bamfile, "rb") as samfile:
        for pileupcolumn in samfile.pileup(
            chrom, 0, end, 
            min_base_quality=base_minQ, 
            stepper="samtools", 
            max_depth=5000000
        ):
            base_pos = pileupcolumn.reference_pos
            
            # 步骤1：收集所有修饰位点覆盖的reads
            if base_pos in mod_dict:
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]
                # 过滤掉空碱基（N或deletion）
                mod_all_reads[base_pos] = [
                    read_name 
                    for read_name, seq in zip(pileupcolumn.get_query_names(), seqs) 
                    if seq  # 只保留非空碱基
                ]
            
            # 步骤2：收集SNP位点信息
            if base_pos in snp_file_dict:
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]
                # 判断是否是杂合位点（至少有两种不同的碱基）
                unique_seqs = set([seq for seq in seqs if seq])
                
                if len(unique_seqs) > 1:
                    # 判断该SNP位点是否本身就是修饰位点
                    if base_pos in mod_dict:
                        # SNP位点本身是修饰位点
                        snp_mod_bases[base_pos] = {
                            read_name: seq 
                            for read_name, seq in zip(pileupcolumn.get_query_names(), seqs) 
                            if seq
                        }
                    else:
                        # SNP位点不是修饰位点
                        snp_bases[base_pos] = {
                            read_name: seq 
                            for read_name, seq in zip(pileupcolumn.get_query_names(), seqs) 
                            if seq
                        }
    
    # 检查是否有有效的SNP位点
    if len(snp_mod_bases) == 0 and len(snp_bases) == 0:
        print(f"{chrom} {strand} 没有杂合的SNP位点")
        return None
    
    # 初始化结果DataFrame
    haplotype_df = pd.DataFrame(columns=[
        "chrom",
        "strand",
        "snp_pos_1base",
        "mod_pos_1base",
        "SNP",
        "A1",
        "A2",
        "MAF",
        "A1_unmod",
        "A2_unmod",
        "A1_mod",
        "A2_mod",
        "ismod?",
        "mod_rate",
    ])
    
    ith = 0
    
    # 分别处理非修饰位点的SNP和修饰位点的SNP
    bases_list = [snp_bases, snp_mod_bases]
    ismod_list = ["No", "Yes"]
    
    for bases, ismod in zip(bases_list, ismod_list):
        result = start_get_haplotypes(
            chrom, strand, mod_all_reads, bases, mod_id_dict, snp_file_dict
        )
        
        if result is None:
            continue
        
        # 将结果添加到DataFrame
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,
                strand,
                i[0] + 1,  # snp_pos (转换为1-based)
                i[1] + 1,  # mod_pos (转换为1-based)
                i[8],      # snpID
                i[2],      # A1
                i[3],      # A2
                i[9],      # maf
                i[4],      # A1_unmod
                i[5],      # A2_unmod
                i[6],      # A1_mod
                i[7],      # A2_mod
                ismod,
                mod_dict[i[1]],  # mod_rate
            )
            ith += 1
    
    if len(haplotype_df) != 0:
        return haplotype_df
    else:
        print(f"{chrom} {strand} 没有A1、A2和unmod、mod对应的read 无法创建列联表")
        return None


def start_get_haplotypes(chrom, strand, mod_all_reads, snp_bases, mod_id_dict, snp_file_dict):
    """
    开始获取haplotype信息
    
    参数:
        chrom: 染色体
        strand: 正负链
        mod_all_reads: 修饰位点的所有reads {pos0: [read1, read2, ...]}
        snp_bases: SNP位点的碱基信息 {pos0: {read1: "A", read2: "T", ...}}
        mod_id_dict: 被修饰的reads {chrom_pos0_strand: "read1;read2;..."}
        snp_file_dict: SNP信息 {pos0: "rsid;A1;A2;MAF"}
    
    返回:
        res_l: 列表，每个元素是一个SNP-mod对的统计结果
    """
    res_l = []
    
    # 遍历所有SNP位点
    for snp_pos, snp_base in snp_bases.items():
        # 遍历所有修饰位点
        for mod_pos, mod_all_reads_l in mod_all_reads.items():
            # 构建mod_id_dict的key
            mod_dict_key = f"{chrom}_{str(mod_pos)}_{strand}"
            
            # 检查该修饰位点是否在mod_id_dict中
            if mod_dict_key not in mod_id_dict:
                continue
            
            # 获取被修饰的reads
            mod_readid = set(mod_id_dict[mod_dict_key].split(";"))
            
            # 获取未修饰的reads（所有覆盖该位点的reads - 被修饰的reads）
            unmod_readid = set(mod_all_reads_l) - mod_readid
            
            # 获取该SNP-mod对的haplotype统计
            res = get_haplotypes(snp_base, mod_readid, unmod_readid, snp_file_dict[snp_pos])
            
            if res is not None:
                A1, A2, A1_unmod, A2_unmod, A1_mod, A2_mod, snpID, maf = res
                res_l.append([
                    snp_pos,    # 0
                    mod_pos,    # 1
                    A1,         # 2
                    A2,         # 3
                    A1_unmod,   # 4
                    A2_unmod,   # 5
                    A1_mod,     # 6
                    A2_mod,     # 7
                    snpID,      # 8
                    maf,        # 9
                ])
    
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes(snp_base, mod_readid, unmod_readid, snp_file_dict_res):
    """
    获取单个SNP-mod对的haplotype统计
    
    参数:
        snp_base: SNP位点的read-碱基字典 {read1: "A", read2: "T", ...}
        mod_readid: 被修饰的read集合
        unmod_readid: 未修饰的read集合
        snp_file_dict_res: SNP信息字符串 "rsid;A1;A2;MAF"
    
    返回:
        如果有效则返回 (A1, A2, A1_unmod, A2_unmod, A1_mod, A2_mod, snpID, maf)
        否则返回 None
    """
    # 获取同时覆盖SNP位点和修饰位点的reads
    intersect_mod_read_set = set(snp_base.keys()) & mod_readid      # SNP覆盖的被修饰reads
    intersect_unmod_read_set = set(snp_base.keys()) & unmod_readid  # SNP覆盖的未修饰reads
    
    # 如果没有同时覆盖的reads，或者修饰/未修饰reads中任何一个为空，则返回None
    if (len(intersect_mod_read_set) == 0) or (len(intersect_unmod_read_set) == 0):
        return None
    
    # 解析SNP信息
    snp_info_parts = snp_file_dict_res.split(";")
    snpID = snp_info_parts[0]
    A1 = snp_info_parts[1]
    A2 = snp_info_parts[2]
    maf = float(snp_info_parts[3])
    
    # 获取携带A1和A2等位基因的reads
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    
    # 确保A1和A2都有对应的reads
    if A1_read and A2_read:
        # 统计四个计数
        A1_unmod = len(A1_read & intersect_unmod_read_set)  # A1等位基因 + 未修饰
        A2_unmod = len(A2_read & intersect_unmod_read_set)  # A2等位基因 + 未修饰
        A1_mod = len(A1_read & intersect_mod_read_set)      # A1等位基因 + 修饰
        A2_mod = len(A2_read & intersect_mod_read_set)      # A2等位基因 + 修饰
        
        return A1, A2, A1_unmod, A2_unmod, A1_mod, A2_mod, snpID, maf
    else:
        # 如果A1或A2没有对应的reads，返回None
        return None



def bayesian_modqtl_analysis(df, min_coverage, n_samples):
    """
    贝叶斯modQTL分析
    
    对每个SNP-修饰位点对进行贝叶斯检验，计算：
    - log_OR_beta: log(OR)的后验均值，作为效应量
    - log_OR_se: log(OR)后验分布的标准差
    - bayes_factor: 贝叶斯因子（H1:关联 vs H0:独立）
    - posterior_prob: P(θ₁ ≠ θ₂ | data)，双侧检验的后验概率
    
    参数:
        df: DataFrame，包含列联表信息
            必须包含列: A1_mod, A1_unmod, A2_mod, A2_unmod
        min_coverage: 最小总覆盖度，列联表四个格子的总和需要 >= 此值
        n_samples: 蒙特卡洛采样数
    
    返回:
        result_df: 原始df加上4个新列
    """
    
    # 初始化结果列
    log_OR_beta_list = []
    log_OR_se_list = []
    bayes_factor_list = []
    posterior_prob_list = []
    fisher_OR_list = []
    fisher_pvalue_list = []  
    
    for idx, row in df.iterrows():
        # 提取列联表
        n11 = int(row["A1_mod"])      # A1 & Mod
        n12 = int(row["A1_unmod"])    # A1 & Unmod
        n21 = int(row["A2_mod"])      # A2 & Mod
        n22 = int(row["A2_unmod"])    # A2 & Unmod
        
        total_coverage = n11 + n12 + n21 + n22
        
        # 检查覆盖度
        if total_coverage < min_coverage:
            log_OR_beta_list.append(np.nan)
            log_OR_se_list.append(np.nan)
            bayes_factor_list.append(np.nan)
            posterior_prob_list.append(np.nan)
            fisher_OR_list.append(np.nan)
            fisher_pvalue_list.append(np.nan)
            continue
        
        # ========== 1. 计算log(OR)的后验分布 ==========
        # 使用Beta后验进行蒙特卡洛采样
        # θ₁ | data ~ Beta(1 + n11, 1 + n12)  # A1的修饰率
        # θ₂ | data ~ Beta(1 + n21, 1 + n22)  # A2的修饰率
        
        alpha1_post = 1 + n11
        beta1_post = 1 + n12
        alpha2_post = 1 + n21
        beta2_post = 1 + n22
        
        # 采样
        theta1_samples = np.random.beta(alpha1_post, beta1_post, n_samples)
        theta2_samples = np.random.beta(alpha2_post, beta2_post, n_samples)
        
        # 计算log(OR)
        # OR = (θ₁/(1-θ₁)) / (θ₂/(1-θ₂))
        # 为避免数值问题，限制theta在(0.001, 0.999)
        theta1_samples = np.clip(theta1_samples, 0.001, 0.999)
        theta2_samples = np.clip(theta2_samples, 0.001, 0.999)
        
        log_OR_samples = (np.log(theta1_samples) - np.log(1 - theta1_samples) - 
                         np.log(theta2_samples) + np.log(1 - theta2_samples))
        
        log_OR_beta = np.mean(log_OR_samples)
        log_OR_se = np.std(log_OR_samples)
        
        log_OR_beta_list.append(log_OR_beta)
        log_OR_se_list.append(log_OR_se)
        
        # ========== 2. 计算后验概率 P(θ₁ ≠ θ₂) ==========
        delta_samples = theta1_samples - theta2_samples
        
        # 双侧检验：取 max(P(δ>0), P(δ<0))
        prob_positive = np.mean(delta_samples > 0)
        prob_negative = np.mean(delta_samples < 0)
        posterior_prob = max(prob_positive, prob_negative)
        
        posterior_prob_list.append(posterior_prob)
        
        # ========== 3. 计算贝叶斯因子 ==========
        bf = compute_independence_bayes_factor(n11, n12, n21, n22)
        bayes_factor_list.append(bf)

        # ========== 4. 计算Fisher精确检验 ==========
        fisher_or, fisher_p = compute_fisher_test(n11, n12, n21, n22)
        fisher_OR_list.append(fisher_or)
        fisher_pvalue_list.append(fisher_p)
    
    # 构建结果DataFrame
    result_df = df.copy()
    result_df["log_OR_beta"] = log_OR_beta_list
    result_df["log_OR_se"] = log_OR_se_list
    result_df["bayes_factor"] = bayes_factor_list
    result_df["posterior_prob"] = posterior_prob_list
    result_df["fisher_OR"] = fisher_OR_list
    result_df["fisher_pvalue"] = fisher_pvalue_list
    
    # 过滤掉无效行
    result_df = result_df[result_df["bayes_factor"].notna()].reset_index(drop=True)
    cols = ["log_OR_beta", "log_OR_se", "bayes_factor", "posterior_prob", "fisher_OR", "fisher_pvalue"]
    result_df[cols] = result_df[cols].round(4)
    
    return result_df


def compute_independence_bayes_factor(n11, n12, n21, n22):
    """
    计算2×2列联表独立性检验的贝叶斯因子
    
    使用 Jeffreys 先验 (Beta(0.5, 0.5)):
    """
    
    # H1: 两个等位基因有各自独立的修饰率
    # θ₁ ~ Beta(0.5, 0.5), θ₂ ~ Beta(0.5, 0.5)
    # P(data|H1) = Beta(n11+0.5, n12+0.5) * Beta(n21+0.5, n22+0.5)
    log_bf_h1 = (betaln(n11 + 0.5, n12 + 0.5) + 
                 betaln(n21 + 0.5, n22 + 0.5) - 
                 2 * betaln(0.5, 0.5))
    
    # H0: 两个等位基因共享相同的修饰率
    # θ ~ Beta(0.5, 0.5)
    # P(data|H0) = Beta(n11+n21+0.5, n12+n22+0.5)
    n_mod_total = n11 + n21
    n_unmod_total = n12 + n22
    log_bf_h0 = (betaln(n_mod_total + 0.5, n_unmod_total + 0.5) - 
                 betaln(0.5, 0.5))
    
    # 计算Bayes Factor
    log_bf = log_bf_h1 - log_bf_h0
    bf = np.exp(log_bf)
    
    return bf


def compute_fisher_test(n11, n12, n21, n22):
    """
    对2×2列联表进行Fisher精确检验
    """
    from scipy.stats import fisher_exact
    
    table = [[n11, n12], [n21, n22]]
    odds_ratio, pvalue = fisher_exact(table, alternative='two-sided')
    
    return odds_ratio, pvalue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Nanopore direct RNA data call modQTL (improved version without WGS EAF prior)."
    )
    parser.add_argument("-b", "--bam", type=str, help="bam file path")
    parser.add_argument(
        "--snp_info",
        type=str,
        help="processed snp txt file path. 需要是\t间隔的文件，里面需要有CHR SNP POS A1 A2 MAF这几列",
    )
    parser.add_argument("-o", "--outdirpre", type=str, help="outdir and pre")
    parser.add_argument(
        "-csv",
        "--modification",
        type=str,
        help="*_sites_result.csv 路径，由nanornaqtl pheno得到的结果",
    )
    parser.add_argument(
        "-pkl",
        "--read_mod_dict",
        type=str,
        help="*_read_final.pkl 路径，由nanornaqtl pheno得到的结果",
    )
    parser.add_argument("-c", "--chrom", type=str, help="chromosome")
    parser.add_argument("-s", "--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument(
        "-q",
        "--min_qscore",
        type=int,
        default=10,
        help="Base Min Query Quality(default=10)",
    )
    parser.add_argument(
        "-m",
        "--modification_type",
        type=str,
        required=True,
        choices=["m6A", "m5C", "pseU", "inosine"],
        help="修饰类型",
    )
    parser.add_argument(
        "--min_coverage",
        type=int,
        default=8,
        help="最小总覆盖度(default=8)",
    )
    parser.add_argument(
        "--mcmc_samples",
        type=int,
        default=2000,
        help="MCMC采样数(default=2000)",
    )
    args = parser.parse_args()

    output_path = f"{args.outdirpre}_{args.modification_type}QTL_{args.chrom}_{args.strand}.csv"

    base_minQ = args.min_qscore - 1 if args.min_qscore != 0 else 0

    # snp_info txt2dict
    snp_info = pd.read_csv(
        args.snp_info, sep="\t", usecols=["CHR", "SNP", "POS", "A1", "A2", "MAF"]
    )
    snp_info = snp_info[snp_info["CHR"] == args.chrom]
    snp_info["pos0"] = snp_info["POS"].astype(int) - 1
    snp_info["v"] = (
        snp_info["SNP"]
        + ";"
        + snp_info["A1"]
        + ";"
        + snp_info["A2"]
        + ";"
        + snp_info["MAF"].astype(str)
    )
    snp_dict = dict(zip(snp_info["pos0"], snp_info["v"]))

    # mod_df csv2dict
    mod_df = pd.read_csv(args.modification)
    mod_df = mod_df[(mod_df["strand"] == args.strand) & (mod_df["chrom"] == args.chrom)]
    if len(mod_df) == 0:
        print(f"{args.chrom} {args.strand} 没有{args.modification_type}位点")
        exit()
    mod_df["pos_0base"] = mod_df["pos_1base"] - 1
    mod_dict = dict(
        x for x in mod_df.apply(create_mod_dict, axis=1).tolist() if x is not None
    )

    # 读取mod_id_dict
    with open(args.read_mod_dict, "rb") as file:
        mod_id_dict = pickle.load(file)

    # 读取基因组大小信息
    geno_size_df = pd.read_csv(
        args.geno_size, sep="\t", header=None, names=["chrom", "size"]
    )
    geno_chrom_end = int(
        geno_size_df[geno_size_df["chrom"] == args.chrom]["size"].values[0]
    )

    # 初始化结果DataFrame
    haplotype_df = pd.DataFrame(
        columns=[
            "chrom",
            "strand",
            "snp_pos_1base",
            "mod_pos_1base",
            "SNP",
            "A1",
            "A2",
            "MAF",
            "A1_unmod",
            "A2_unmod",
            "A1_mod",
            "A2_mod",
            "ismod?",
            "mod_rate",
        ]
    )

    # 统计haplotype信息
    print(f"开始统计 {args.chrom} {args.strand} 的haplotype信息...")
    haplotype_result = count_haplotype(
        args.chrom,
        geno_chrom_end,
        args.strand,
        args.bam,
        mod_dict,
        mod_id_dict,
        snp_dict,
        base_minQ,
    )
    
    if haplotype_result is not None:
        haplotype_df = pd.concat([haplotype_df, haplotype_result], ignore_index=True)
    
    haplotype_df = haplotype_df.dropna().reset_index(drop=True)
    
    if len(haplotype_df) == 0:
        print(f"{args.chrom} {args.strand} 没有有效的haplotype数据")
        exit()
    
    print(f"成功获取 {len(haplotype_df)} 个SNP-{args.modification_type}对的haplotype信息")
    print(haplotype_df.head())
    
    # 按染色体、SNP位置、修饰位置排序
    df = haplotype_df.sort_values(by=["chrom", "snp_pos_1base", "mod_pos_1base"])
    
    # 进行贝叶斯modQTL分析
    res_df = bayesian_modqtl_analysis(
        df, 
        min_coverage=args.min_coverage,
        n_samples=args.mcmc_samples
    )
    
    if len(res_df) == 0:
        print(f"{args.chrom} {args.strand} 中的SNP没有足够的read覆盖度，无法进行统计")
        exit()

    # ========== 计算Fisher的FDR校正 ==========
    from statsmodels.stats.multitest import multipletests
    
    # 只对有有效p值的行进行FDR校正
    valid_pvalues = res_df['fisher_pvalue'].notna()
    
    if valid_pvalues.sum() > 0:
        pvalues = res_df.loc[valid_pvalues, 'fisher_pvalue'].values
        _, fdr_pvalues, _, _ = multipletests(pvalues, method='fdr_bh')

        res_df.loc[valid_pvalues, 'fisher_FDR'] = fdr_pvalues
        res_df.loc[~valid_pvalues, 'fisher_FDR'] = np.nan
    else:
        res_df['fisher_FDR'] = np.nan
    
    # 保存结果
    res_df = res_df.reset_index(drop=True)
    res_df.to_csv(output_path, index=False)
    print(f"\n{output_path} 已保存")
