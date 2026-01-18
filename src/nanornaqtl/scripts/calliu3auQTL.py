import pandas as pd
import pysam
import scipy
import argparse
from collections import Counter
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats
from scipy.special import betaln, gammaln, loggamma
from scipy.optimize import minimize_scalar,minimize
import pymc as pm
import arviz as az
from scipy.stats import gaussian_kde
from scipy.integrate import quad
import warnings
from scipy.stats import chi2_contingency

warnings.filterwarnings('ignore')

def count_haplotype(
    chrom,
    end,
    strand,
    bamfile,
    read_dict,
    snp_file_dict,
    base_minQ,
):
    # read_dict read对应的APA read:APA
    # snp_file_dict snp位点信息 int(pos0): rsid;A1;A2;MAF
    snp_bases = {}
    # snp_bases 存放的是snp位点 [pos0]:{read1:"A",read2:"A",read3:"T"}
    with pysam.AlignmentFile(bamfile, "rb") as samfile:
        for pileupcolumn in samfile.pileup(
            chrom,
            0,
            end,
            min_base_quality=base_minQ,
            stepper="samtools",
            max_depth=5000000,
        ):
            base_pos = pileupcolumn.reference_pos
            ## 处理snp位点的情况
            # 1判断位置是否是snp
            if base_pos in snp_file_dict:
                # 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]
                if len(set([seq for seq in seqs if seq])) > 1:
                    snp_bases[base_pos] = {
                        i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                    }
    if len(snp_bases) == 0:
        print(f"{chrom} {strand}中,没有杂合的snp被read覆盖")
        return None
    haplotype_df = pd.DataFrame(
        columns=[
            "chrom",
            "strand",
            "snp_pos_1base",
            "SNP",
            "A1",
            "A2",
            "MAF",
            "all_types",
            "A1_counts",
            "A2_counts",
        ]
    )
    result = start_get_haplotypes(snp_bases, read_dict, snp_file_dict)
    if result is None:
        print(f"{chrom} {strand}中没有符合条件的snp")
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,
                strand,
                int(i[7]) + 1,
                i[5],
                i[0],
                i[1],
                i[6],
                i[2],
                i[3],
                i[4],
            )
            ith += 1
        return haplotype_df

def start_get_haplotypes(snp_bases, read_dict, snp_file_dict):
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_dict read:APA
    res_l = []
    for snp_pos, snp_base in snp_bases.items():
        # 筛选出snp覆盖的read中有APA分类的read
        snp_base_filter = {
            key: value for key, value in snp_base.items() if key in read_dict
        }
        if len(snp_base_filter) == 0:
            continue
        res = get_haplotypes(snp_base_filter, read_dict, snp_file_dict[snp_pos])
        if res is not None:
            a1, a2, all_types, A1_counts, A2_counts, snpID, maf = res
            res_l.append([a1, a2, all_types, A1_counts, A2_counts, snpID, maf, snp_pos])
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes(snp_base, read_dict, snp_file_dict_res):
    """
    返回结构:
    - A1, A2: 等位基因标识
    - all_types: 所有分子表型的列表（排序后）
    - A1_counts, A2_counts: 对应each type的计数列表
    - snpID, maf: SNP信息
    """
    snpID, A1, A2, maf = snp_file_dict_res.split(";")[:4]
    maf = float(maf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if not A1_read or not A2_read:
        return None
    # 获取这些reads对应的分子表型
    A1_phenotypes = [read_dict[read_id] for read_id in A1_read if read_id in read_dict]
    A2_phenotypes = [read_dict[read_id] for read_id in A2_read if read_id in read_dict]
    if len(A1_phenotypes) == 0 or len(A2_phenotypes) == 0:
        return None

    # 获取所有unique的分子表型，并按字典序排序
    all_types = sorted(list(set(A1_phenotypes + A2_phenotypes)))

    # 计算每个表型在A1和A2中的计数
    A1_counter = Counter(A1_phenotypes)
    A2_counter = Counter(A2_phenotypes)

    # 按照all_types的顺序构建计数列表
    A1_counts = [A1_counter.get(ptype, 0) for ptype in all_types]
    A2_counts = [A2_counter.get(ptype, 0) for ptype in all_types]

    return A1, A2, all_types, A1_counts, A2_counts, snpID, maf

def analyze_single_snp(A1_counts, A2_counts, min_coverage=10, n_samples=1000):
    """
    单个SNP的usage QTL贝叶斯分析
    
    对2×K列联表进行贝叶斯检验，计算：
    - bayes_factor: 贝叶斯因子（H1:关联 vs H0:独立）
    - posterior_prob: P(π₁ ≠ π₂ | data)，分布差异的后验概率
    - TVD_beta: Total Variation Distance，效应大小
    - TVD_se: TVD的后验标准差
    - dominant_shift: 主要变化方向描述
    
    参数:
        A1_counts: list，A1等位基因对应各isoform的计数
        A2_counts: list，A2等位基因对应各isoform的计数
        min_coverage: 最小总覆盖度
        n_samples: 蒙特卡洛采样数
    
    返回:
        dict: 包含各统计量的字典，若不满足条件返回None
    """
    
    A1_counts = np.array(A1_counts)
    A2_counts = np.array(A2_counts)
    
    n_A1 = np.sum(A1_counts)
    n_A2 = np.sum(A2_counts)
    total_coverage = n_A1 + n_A2
    
    # 检查覆盖度
    if total_coverage < min_coverage:
        return None
    
    # 检查是否至少有2种isoform
    K = len(A1_counts)
    if K < 2:
        return None
    
    # 检查A1和A2是否都有reads
    if n_A1 == 0 or n_A2 == 0:
        return None
    
    # ========== 1. 蒙特卡洛采样计算后验分布 ==========
    # 使用Dirichlet后验
    # π₁ | data ~ Dirichlet(1 + A1_counts)
    # π₂ | data ~ Dirichlet(1 + A2_counts)
    
    alpha1_post = 1 + A1_counts  # Dirichlet参数
    alpha2_post = 1 + A2_counts
    
    # 采样
    pi1_samples = np.random.dirichlet(alpha1_post, n_samples)  # [n_samples, K]
    pi2_samples = np.random.dirichlet(alpha2_post, n_samples)  # [n_samples, K]
    
    # ========== 2. 计算TVD的后验分布 ==========
    # TVD = 0.5 * Σ|π₁ₖ - π₂ₖ|
    diff_samples = pi1_samples - pi2_samples  # [n_samples, K]
    TVD_samples = 0.5 * np.sum(np.abs(diff_samples), axis=1)  # [n_samples]
    
    TVD_beta = np.mean(TVD_samples)
    TVD_se = np.std(TVD_samples)
    
    # ========== 3. 计算后验概率 ==========
    # P(π₁ ≠ π₂) 用TVD > 小阈值来近似
    # 这里用一个很小的阈值，表示"实质性差异"
    posterior_prob = np.mean(TVD_samples > 0.01)
    
    # ========== 4. 计算贝叶斯因子 ==========
    bf = compute_independence_bf_2xK(A1_counts, A2_counts)

    # ========== 5. 计算G检验 ==========
    g_pvalue = compute_g_test(A1_counts, A2_counts)
    
    # ========== 6. 计算dominant_shift ==========
    # 基于后验均值计算每个isoform的比例变化
    pi1_mean = np.mean(pi1_samples, axis=0)  # [K]
    pi2_mean = np.mean(pi2_samples, axis=0)  # [K]
    diff_mean = pi1_mean - pi2_mean  # A1相对于A2的变化
    
    dominant_shift = get_dominant_shift(diff_mean)
    
    return {
        "bayes_factor": bf,
        "posterior_prob": posterior_prob,
        "TVD_beta": TVD_beta,
        "TVD_se": TVD_se,
        "dominant_shift": dominant_shift,
        "g_pvalue": g_pvalue,
    }


def compute_independence_bf_2xK(A1_counts, A2_counts):
    """
    2×K列联表独立性检验的贝叶斯因子 (Jeffreys prior)
    BF = P(data|H1) / P(data|H0)
    """
    A1_counts = np.array(A1_counts)
    A2_counts = np.array(A2_counts)
    K = len(A1_counts)
    
    n_A1 = np.sum(A1_counts)
    n_A2 = np.sum(A2_counts)
    col_sums = A1_counts + A2_counts
    N = n_A1 + n_A2
    
    # 使用 Jeffreys prior: Dirichlet(0.5, ..., 0.5)
    alpha = 0.5
    
    # H1: 独立的 Dirichlet-Multinomial
    log_p_h1 = (
        gammaln(K * alpha) - gammaln(n_A1 + K * alpha) + 
        np.sum(gammaln(A1_counts + alpha) - gammaln(alpha)) +
        gammaln(K * alpha) - gammaln(n_A2 + K * alpha) + 
        np.sum(gammaln(A2_counts + alpha) - gammaln(alpha))
    )
    
    # H0: 共享分布，数据是从同一个multinomial中独立抽取
    log_p_h0 = (
        gammaln(K * alpha) - gammaln(N + K * alpha) +
        np.sum(gammaln(col_sums + alpha) - gammaln(alpha))
    )
    
    log_bf = log_p_h1 - log_p_h0
    
    return np.exp(np.clip(log_bf, -700, 700))


def compute_g_test(A1_counts, A2_counts):
    """G检验"""
    contingency_table = np.array([A1_counts, A2_counts])
    try:
        g_stat, g_pvalue, dof, expected = stats.chi2_contingency(
            contingency_table, lambda_="log-likelihood"
        )
        return g_pvalue
    except ValueError:
        return np.nan


def get_dominant_shift(diff_mean):
    """
    根据isoform比例差异生成dominant_shift描述
    
    参数:
        diff_mean: π₁ - π₂ 的后验均值向量，正值表示A1中该isoform比例更高
    
    返回:
        str: 描述主要变化方向的字符串
    """
    
    # 找出变化最大的isoform（绝对值最大）
    abs_diff = np.abs(diff_mean)
    max_idx = np.argmax(abs_diff)
    max_diff = diff_mean[max_idx]
    
    # 只报告变化超过一定阈值的isoform
    threshold = 0.05  # 5%的比例变化
    
    significant_changes = []
    for i, d in enumerate(diff_mean):
        if abs(d) >= threshold:
            direction = "A1↑" if d > 0 else "A1↓"
            significant_changes.append(f"type{i+1}:{direction}({d:+.2f})")
    
    if len(significant_changes) == 0:
        # 如果没有超过阈值的，报告最大的那个
        direction = "A1↑" if max_diff > 0 else "A1↓"
        return f"type{max_idx+1}:{direction}({max_diff:+.2f})"
    else:
        return "; ".join(significant_changes)


def process_all_data_simple(df, min_coverage, n_samples):
    """
    处理所有SNP的usage QTL分析
    
    参数:
        df: DataFrame，包含列: A1_counts, A2_counts, all_types等
        min_coverage: 最小总覆盖度
        n_samples: 蒙特卡洛采样数
    
    返回:
        result_df: 原始df加上新列
    """
    
    # 初始化结果列
    bf_list = []
    posterior_prob_list = []
    TVD_beta_list = []
    TVD_se_list = []
    dominant_shift_list = []
    g_pvalue_list = []
    
    for idx, row in df.iterrows():
        # 提取计数（可能是字符串格式，需要转换）
        A1_counts = row["A1_counts"]
        A2_counts = row["A2_counts"]
        
        # 如果是字符串格式，转换为列表
        if isinstance(A1_counts, str):
            A1_counts = eval(A1_counts)
        if isinstance(A2_counts, str):
            A2_counts = eval(A2_counts)
        
        # 分析
        result = analyze_single_snp(
            A1_counts, A2_counts,
            min_coverage=min_coverage,
            n_samples=n_samples
        )
        
        if result is None:
            bf_list.append(np.nan)
            posterior_prob_list.append(np.nan)
            TVD_beta_list.append(np.nan)
            TVD_se_list.append(np.nan)
            dominant_shift_list.append(np.nan)
            g_pvalue_list.append(np.nan)
        else:
            bf_list.append(result["bayes_factor"])
            posterior_prob_list.append(result["posterior_prob"])
            TVD_beta_list.append(result["TVD_beta"])
            TVD_se_list.append(result["TVD_se"])
            dominant_shift_list.append(result["dominant_shift"])
            g_pvalue_list.append(result["g_pvalue"])
    
    # 构建结果DataFrame
    result_df = df.copy()
    result_df["bayes_factor"] = bf_list
    result_df["posterior_prob"] = posterior_prob_list
    result_df["TVD_beta"] = TVD_beta_list
    result_df["TVD_se"] = TVD_se_list
    result_df["dominant_shift"] = dominant_shift_list
    result_df["g_pvalue"] = g_pvalue_list
    
    # 过滤掉无效行
    result_df = result_df[result_df["bayes_factor"].notna()].reset_index(drop=True)

    # G检验FDR校正
    valid_pvalues = result_df['g_pvalue'].notna()
    if valid_pvalues.sum() > 0:
        _, fdr_values, _, _ = multipletests(
            result_df.loc[valid_pvalues, 'g_pvalue'].values, method='fdr_bh'
        )
        result_df.loc[valid_pvalues, 'g_FDR'] = fdr_values
    else:
        result_df['g_FDR'] = np.nan
    
    return result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Nanopore direct RNA data call 3auQTL and iuQTL."
    )
    parser.add_argument("-b", "--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-o", "--outdirpre", type=str, help="outdir and pre")
    parser.add_argument(
        "-f", "--read_overlap_file", type=str, help="read2apa(*_APA_result.csv) or read2isoform(OUT.transcript_model_reads.tsv.gz) overlap file path"
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
    parser.add_argument("-m", "--molecular_type", type=str, required=True,choices=["APA", "isoform"], help="分子表型类型")
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

    output_path = (
        f"{args.outdirpre}_{args.molecular_type}_haplotype_{args.chrom}_{args.strand}_tmp.csv"
    )

    base_minQ = args.min_qscore - 1 if args.min_qscore != 0 else 0

    if args.molecular_type == "APA":
        overlap_df = pd.read_csv(args.read_overlap_file)
        df_dict = dict(zip(overlap_df["readID"], overlap_df["APA_type"]))  # readid:APA
    else:
        overlap_df = pd.read_csv(args.read_overlap_file, sep="\t", compression="gzip")
        overlap_df.columns = ["readID", "isoformID"]
        overlap_df = overlap_df[overlap_df["isoformID"] != "*"]
        df_dict = dict(zip(overlap_df["readID"], overlap_df["isoformID"]))  # readid:isoform

    # snp_info txt2dict
    snp_info = pd.read_csv(args.snp_info, sep="\t", usecols=["CHR","SNP","POS","A1","A2","MAF"])
    snp_info = snp_info[snp_info['CHR'] == args.chrom]
    snp_info["pos0"] = snp_info["POS"].astype(int) - 1
    snp_info["v"] = snp_info["SNP"] + ";" + snp_info["A1"] + ";" + snp_info["A2"] + ";" + snp_info["MAF"].astype(str)
    snp_dict = dict(zip(snp_info["pos0"], snp_info["v"])) # pos0: SNP;A1;A2;MAF

    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom","size"])
    geno_chrom_end = int(geno_size_df[geno_size_df['chrom'] == args.chrom]['size'].values[0])

    haplotype_df = pd.DataFrame(
            columns=[
                "chrom",
                "strand",
                "snp_pos_1base",
                "SNP",
                "A1",
                "A2",
                "MAF",
                "all_types",
                "A1_counts",
                "A2_counts",
            ]
        )

    haplotype_df = pd.concat(
        [
            haplotype_df,
            count_haplotype(
                args.chrom,
                geno_chrom_end,
                args.strand,
                args.bam,
                df_dict,
                snp_dict,
                base_minQ,
            ),
        ],
        ignore_index=True,
    )
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=["chrom", "snp_pos_1base"])
        df = process_all_data_simple(df, args.min_coverage, args.mcmc_samples)
        if df is None:
            print(
                f"{args.chrom} {args.strand}中的SNP没有足够的read覆盖度，无法进行统计"
            )
            exit()
        df = df.reset_index(drop=True)
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
