import pandas as pd
import pysam
import argparse
import numpy as np
from scipy import stats
import warnings
import ast
import numpy as np
import json
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore")


def count_haplotype(
    chrom, end, strand, bamfile, read_st_dict, snp_file_dict, base_minQ
):
    # read_st_dict read对应的稳定性指标 readid:rate
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
            "A1_IR",
            "A2_IR",
        ]
    )
    result = start_get_haplotypes(snp_bases, read_st_dict, snp_file_dict)
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


def start_get_haplotypes(snp_bases, read_st_dict, snp_file_dict):
    res_l = []
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_st_dict readid:rate
    for snp_pos, snp_base in snp_bases.items():
        # 筛选出snp覆盖的read有稳定性分数的read
        snp_base_filter = {
            key: value for key, value in snp_base.items() if key in read_st_dict
        }
        if len(snp_base_filter) == 0:
            continue
        res = get_haplotypes(snp_base_filter, read_st_dict, snp_file_dict[snp_pos])
        if res is not None:
            a1, a2, A1_ir, A2_ir, snpID, maf = res
            res_l.append([a1, a2, A1_ir, A2_ir, snpID, maf, snp_pos])
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes(snp_base, read_st_dict, snp_file_dict_res):
    snpID, A1, A2, maf = snp_file_dict_res.split(";")[:4]
    maf = float(maf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if A1_read and A2_read:
        A1_ir = [
            read_st_dict[read_id] for read_id in A1_read if read_id in read_st_dict
        ]
        A2_ir = [
            read_st_dict[read_id] for read_id in A2_read if read_id in read_st_dict
        ]
        if len(A1_ir) == 0 or len(A2_ir) == 0:
            return None
        return A1, A2, A1_ir, A2_ir, snpID, maf
    else:
        return None


def bayesian_irqtl_single_snp(row, mcmc_samples=2000, random_seed=42):
    """
    基于分布差异的贝叶斯内含子滞留率QTL分析
    """
    np.random.seed(random_seed)
    
    # 数据解析
    try:
        A1_ir = ast.literal_eval(row['A1_IR']) if isinstance(row['A1_IR'], str) else row['A1_IR']
        A2_ir = ast.literal_eval(row['A2_IR']) if isinstance(row['A2_IR'], str) else row['A2_IR']
    except:
        A1_ir, A2_ir = row['A1_IR'], row['A2_IR']
    
    A1_ir = np.array(A1_ir, dtype=np.float64)
    A2_ir = np.array(A2_ir, dtype=np.float64)
    A1_ir = A1_ir[np.isfinite(A1_ir) & (A1_ir >= 0) & (A1_ir <= 1)]
    A2_ir = A2_ir[np.isfinite(A2_ir) & (A2_ir >= 0) & (A2_ir <= 1)]
    
    if len(A1_ir) < 5 or len(A2_ir) < 5:
        return None
    
    # 检查：如果合并数据90%以上都是0，跳过
    combined_all = np.concatenate([A1_ir, A2_ir])
    if np.mean(combined_all == 0) > 0.9:
        return None
    
    n1, n2 = len(A1_ir), len(A2_ir)
    m1, m2 = np.mean(A1_ir), np.mean(A2_ir)
    
    # KS统计量
    observed_ks, _ = stats.ks_2samp(A1_ir, A2_ir)
    
    # 贝叶斯排列检验
    combined = np.concatenate([A1_ir, A2_ir])
    n_perm = mcmc_samples
    perm_ks = np.zeros(n_perm)
    
    for i in range(n_perm):
        np.random.shuffle(combined)
        perm_A1 = combined[:n1]
        perm_A2 = combined[n1:]
        perm_ks[i], _ = stats.ks_2samp(perm_A1, perm_A2)
    
    # 计算BF
    prob_h0 = np.mean(perm_ks >= observed_ks)
    prob_h0 = np.clip(prob_h0, 0.001, 0.999)
    
    posterior_odds = (1 - prob_h0) / prob_h0
    bf = posterior_odds  # prior_odds = 1
    bf = min(bf, 1000)
    
    # 效应量：均值差异的logit变换
    m1_clip = np.clip(m1, 1e-6, 1 - 1e-6)
    m2_clip = np.clip(m2, 1e-6, 1 - 1e-6)
    beta_est = np.log(m1_clip / (1 - m1_clip)) - np.log(m2_clip / (1 - m2_clip))
    
    # SE (delta method for logit difference)
    v1, v2 = np.var(A1_ir, ddof=1), np.var(A2_ir, ddof=1)
    se1 = np.sqrt(v1 / n1) / (m1_clip * (1 - m1_clip)) if m1_clip > 0 and m1_clip < 1 else 0
    se2 = np.sqrt(v2 / n2) / (m2_clip * (1 - m2_clip)) if m2_clip > 0 and m2_clip < 1 else 0
    se_est = np.sqrt(se1**2 + se2**2)
    
    # 频率派参考
    try:
        _, p_welch = stats.ttest_ind(A1_ir, A2_ir, equal_var=False)
        _, p_mw = stats.mannwhitneyu(A1_ir, A2_ir, alternative='two-sided')
        _, p_ks = stats.ks_2samp(A1_ir, A2_ir)
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

def process_all_irqtl_data(haplotype_df, min_coverage, mcmc_samples, random_seed=42):
    """
    遍历所有SNP进行贝叶斯内含子滞留率QTL分析
    
    参数:
        haplotype_df: 包含SNP信息和内含子滞留率的DataFrame
        min_coverage: 最小总覆盖度阈值
        mcmc_samples: MCMC采样数
        random_seed: 随机种子基数
        
    返回:
        DataFrame: 包含所有SNP的分析结果
    """
    print(f"开始贝叶斯内含子滞留率QTL分析...")
    print(f"参数: min_coverage={min_coverage}, mcmc_samples={mcmc_samples}")
    
    # ===== 1. 数据质控 =====
    print(f"\n步骤1: 数据质控...")
    
    def count_valid_ir(ir_list):
        """计算有效内含子滞留率数量"""
        if isinstance(ir_list, str):
            try:
                ir_list = ast.literal_eval(ir_list)
            except:
                return 0
        if not isinstance(ir_list, (list, np.ndarray)):
            return 0
        ir_arr = np.array(ir_list)
        # 有效值：0 <= ir < 1 且非NaN
        return np.sum(np.isfinite(ir_arr) & (ir_arr >= 0) & (ir_arr < 1))
    
    haplotype_df = haplotype_df.copy()
    haplotype_df['_n_A1'] = haplotype_df['A1_IR'].apply(count_valid_ir)
    haplotype_df['_n_A2'] = haplotype_df['A2_IR'].apply(count_valid_ir)
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
        
        bayes_result = bayesian_irqtl_single_snp(
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
                'A1_IR': row['A1_IR'],
                'A2_IR': row['A2_IR'],
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
        description="Nanopore direct RNA data call irQTL(QTLs related to intron retention rate)."
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
        "--ir_csv",
        type=str,
        help="*_intronRetention_result.csv路径，由nanornaqtl pheno得到的结果",
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
        f"{args.outdirpre}_intronretention_haplotype_{args.chrom}_{args.strand}_tmp.csv"
    )

    base_minQ = args.min_qscore - 1 if args.min_qscore != 0 else 0

    # csv2dict
    read_st = pd.read_csv(args.ir_csv, usecols=[0, 4])
    read_st_dict = read_st.set_index("readID")[
        "IntronRetentionRate"
    ].to_dict()  # readid:rate
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
            "A1_IR",
            "A2_IR",
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
                read_st_dict,
                snp_dict,
                base_minQ,
            ),
        ],
        ignore_index=True,
    )
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=["chrom", "snp_pos_1base"])
        df = process_all_irqtl_data(
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
