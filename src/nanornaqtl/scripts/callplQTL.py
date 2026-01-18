import pandas as pd
import pysam
import argparse
import numpy as np
from scipy import stats
import warnings
import ast
import numpy as np
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore")


def count_haplotype(
    chrom, end, strand, bamfile, read_len_dict, snp_file_dict, base_minQ
):
    # read_len_dict read对应的polyA长度 readid:length
    # snp_file_dict snp位点信息 # pos0: rsid;A1;A2;MAF
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
                seqs = [
                    i.upper() for i in pileupcolumn.get_query_sequences()
                ]  # 被覆盖位点处的所有碱基
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
            "A1_len",
            "A2_len",
        ]
    )
    result = start_get_haplotypes(snp_bases, read_len_dict, snp_file_dict)
    if result is None:
        print(f"{chrom} {strand}中没有符合条件的snp")
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,
                strand,
                int(i[6]) + 1,
                i[4],
                i[0],
                i[1],
                i[5],
                i[2],
                i[3],
            )
            ith += 1
        return haplotype_df


def start_get_haplotypes(snp_bases, read_len_dict, snp_file_dict):
    res_l = []
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_len_dict readid:length
    for snp_pos, snp_base in snp_bases.items():
        # 筛选出snp覆盖的read有稳定性分数的read
        snp_base_filter = {
            key: value for key, value in snp_base.items() if key in read_len_dict
        }
        if len(snp_base_filter) == 0:
            continue
        res = get_haplotypes(snp_base_filter, read_len_dict, snp_file_dict[snp_pos])
        if res is not None:
            a1, a2, A1_len, A2_len, snpID, maf = res
            res_l.append([a1, a2, A1_len, A2_len, snpID, maf, snp_pos])
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes(snp_base, read_len_dict, snp_file_dict_res):
    snpID, A1, A2, maf = snp_file_dict_res.split(";")[:4]
    maf = float(maf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if A1_read and A2_read:
        A1_len = [
            read_len_dict[read_id] for read_id in A1_read if read_id in read_len_dict
        ]
        A2_len = [
            read_len_dict[read_id] for read_id in A2_read if read_id in read_len_dict
        ]
        if len(A1_len) == 0 or len(A2_len) == 0:
            return None
        return A1, A2, A1_len, A2_len, snpID, maf
    else:
        return None


def bayesian_plqtl_single_snp(row, mcmc_samples=2000, random_seed=42):
    """
    基于分布重叠度的贝叶斯检验
    """
    np.random.seed(random_seed)
    
    # 数据清理
    try:
        A1_l = ast.literal_eval(row['A1_len']) if isinstance(row['A1_len'], str) else row['A1_len']
        A2_l = ast.literal_eval(row['A2_len']) if isinstance(row['A2_len'], str) else row['A2_len']
    except:
        A1_l, A2_l = row['A1_len'], row['A2_len']
    
    A1_l = np.array(A1_l, dtype=np.float64)
    A2_l = np.array(A2_l, dtype=np.float64)
    A1_l = A1_l[(A1_l > 0) & np.isfinite(A1_l)]
    A2_l = A2_l[(A2_l > 0) & np.isfinite(A2_l)]
    
    if len(A1_l) < 5 or len(A2_l) < 5:
        return None
    
    n1, n2 = len(A1_l), len(A2_l)
    m1, m2 = np.mean(A1_l), np.mean(A2_l)
    
    # ===== 核心：计算观测统计量 =====
    # 使用KS统计量衡量分布差异
    observed_ks, _ = stats.ks_2samp(A1_l, A2_l)
    
    # ===== 贝叶斯排列检验 =====
    combined = np.concatenate([A1_l, A2_l])
    n_perm = mcmc_samples
    perm_ks = np.zeros(n_perm)
    
    for i in range(n_perm):
        np.random.shuffle(combined)
        perm_A1 = combined[:n1]
        perm_A2 = combined[n1:]
        perm_ks[i], _ = stats.ks_2samp(perm_A1, perm_A2)
    
    # ===== 计算BF =====
    # 后验：观测到的KS值在排列分布中的位置
    # 如果观测KS远大于排列KS，说明分布确实不同
    
    # 使用经验分布计算概率
    prob_h0 = np.mean(perm_ks >= observed_ks)  # H0下看到这么大KS的概率
    prob_h0 = np.clip(prob_h0, 0.001, 0.999)   # 避免极端值
    
    # 先验：假设H0和H1各50%
    prior_odds = 1.0
    
    # 后验odds
    # P(H1|data) / P(H0|data) = P(data|H1) / P(data|H0)
    # P(data|H0) ∝ prob_h0
    # P(data|H1) ∝ 1 - prob_h0 (简化)
    posterior_odds = (1 - prob_h0) / prob_h0
    
    bf = posterior_odds / prior_odds
    
    # 效应量：log ratio
    beta_est = np.log(m1) - np.log(m2) if m1 > 0 and m2 > 0 else 0
    se_est = np.sqrt(np.var(A1_l)/(n1*m1**2) + np.var(A2_l)/(n2*m2**2)) if m1 > 0 and m2 > 0 else 0
    
    # 频率派参考
    try:
        _, p_welch = stats.ttest_ind(A1_l, A2_l, equal_var=False)
        _, p_mw = stats.mannwhitneyu(A1_l, A2_l, alternative='two-sided')
        _, p_ks = stats.ks_2samp(A1_l, A2_l)
    except:
        p_welch, p_mw, p_ks = 1.0, 1.0, 1.0

    return {
        'beta': round(beta_est, 4),
        'SE': round(se_est, 4),
        'KS_stat': round(observed_ks, 4),
        'bayes_factor': round(bf, 2),
        'posterior_prob': round(prob_h0, 4),
        'p_welch': round(p_welch, 6),
        'p_mw': round(p_mw, 6),
        'p_ks': round(p_ks, 6)
    }


def process_all_plqtl_data(haplotype_df, min_coverage, mcmc_samples, random_seed=42):
    """
    遍历所有SNP进行贝叶斯polyA尾长QTL分析
    
    参数:
        haplotype_df: 包含SNP信息和polyA长度的DataFrame
        min_coverage: 最小总覆盖度阈值
        mcmc_samples: MCMC采样数
        random_seed: 随机种子基数
        
    返回:
        DataFrame: 包含所有SNP的分析结果
    """
    print(f"开始贝叶斯polyA尾长QTL分析...")
    print(f"参数: min_coverage={min_coverage}, mcmc_samples={mcmc_samples}")
    
    # ===== 1. 数据质控 =====
    print(f"\n步骤1: 数据质控...")
    
    def count_valid_lengths(length_list):
        """计算有效尾长数量"""
        if isinstance(length_list, str):
            try:
                length_list = ast.literal_eval(length_list)
            except:
                return 0
        if not isinstance(length_list, (list, np.ndarray)):
            return 0
        lengths = np.array(length_list)
        return np.sum((lengths > 0) & np.isfinite(lengths))
    
    haplotype_df = haplotype_df.copy()
    haplotype_df['_n_A1'] = haplotype_df['A1_len'].apply(count_valid_lengths)
    haplotype_df['_n_A2'] = haplotype_df['A2_len'].apply(count_valid_lengths)
    haplotype_df['_total_cov'] = haplotype_df['_n_A1'] + haplotype_df['_n_A2']
    
    # 过滤
    filtered_df = haplotype_df[haplotype_df['_total_cov'] >= min_coverage].copy()
    filtered_df = filtered_df.reset_index(drop=True)
    
    print(f"  过滤前SNP数量: {len(haplotype_df)}")
    print(f"  过滤后SNP数量: {len(filtered_df)}")
    
    if len(filtered_df) == 0:
        print("警告: 没有SNP通过质控")
        return None
    
    # 清理临时列
    filtered_df = filtered_df.drop(columns=['_n_A1', '_n_A2', '_total_cov'])
    
    # ===== 2. 遍历分析 =====
    print(f"\n步骤2: 贝叶斯分析...")
    results = []
    n_success = 0
    n_fail = 0
    
    for idx, row in filtered_df.iterrows():
        if (idx + 1) % 100 == 0 or idx == 0:
            print(f"  进度: {idx+1}/{len(filtered_df)}")
        
        bayes_result = bayesian_plqtl_single_snp(
            row,
            mcmc_samples=mcmc_samples,
            random_seed=random_seed + idx
        )
        
        if bayes_result is not None:
            result = {
                'chrom': row['chrom'],
                'strand': row['strand'],
                'snp_pos_1base': row['snp_pos_1base'],
                'SNP': row['SNP'],
                'A1': row['A1'],
                'A2': row['A2'],
                'MAF': row['MAF'],
                'A1_len': row['A1_len'],
                'A2_len': row['A2_len'],
            }
            result.update(bayes_result)
            results.append(result)
            n_success += 1
        else:
            n_fail += 1
    
    print(f"\n分析完成: 成功 {n_success}, 失败 {n_fail}")
    
    if len(results) == 0:
        print("错误: 所有SNP分析都失败")
        return None
    
    results_df = pd.DataFrame(results)
    return results_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Nanopore direct RNA data call plQTL(QTLs related to polyA tail length)."
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
        "--polya_csv",
        type=str,
        help="*_polyAlen_result.csv路径，由nanornaqtl pheno得到的结果",
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
        f"{args.outdirpre}_polyAlen_haplotype_{args.chrom}_{args.strand}_tmp.csv"
    )

    base_minQ = args.min_qscore - 1 if args.min_qscore != 0 else 0

    # csv2dict
    read_len = pd.read_csv(args.polya_csv)
    read_len_dict = read_len.set_index("readID")["polyA_length"].to_dict()  # readid:length
    # txt2dict
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
    snp_dict = dict(zip(snp_info["pos0"], snp_info["v"]))  # pos0: SNP;A1;A2;MAF

    geno_size_df = pd.read_csv(
        args.geno_size, sep="\t", header=None, names=["chrom", "size"]
    )
    geno_chrom_end = int(
        geno_size_df[geno_size_df["chrom"] == args.chrom]["size"].values[0]
    )


    haplotype_df = pd.DataFrame(
        columns=[
            "chrom",
            "strand",
            "snp_pos_1base",
            "SNP",
            "A1",
            "A2",
            "MAF",
            "A1_len",
            "A2_len",
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
                read_len_dict,
                snp_dict,
                base_minQ,
            ),
        ],
        ignore_index=True,
    )
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=["chrom", "snp_pos_1base"])
        df = process_all_plqtl_data(
            df,
            min_coverage=args.min_coverage,
            mcmc_samples=args.mcmc_samples,
        )
        if df is None:
            print(
                f"{args.chrom} {args.strand}中的SNP没有足够的read覆盖度，无法进行统计"
            )
            exit()
        df = df.reset_index(drop=True)
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
