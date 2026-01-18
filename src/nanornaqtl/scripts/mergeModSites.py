import os
import glob
import pandas as pd
import pickle
from subprocess import call
import re
from itertools import groupby
import logomaker
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
import warnings


def fasta_iter(fasta_name):
    with open(fasta_name) as filename:
        faiter = (x[1] for x in groupby(filename, lambda line: line[0] == ">"))
        for header in faiter:
            header_str = next(header)[1:].strip()
            seq = "".join(s.strip() for s in next(faiter))
            yield header_str, seq


def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans("ATGCatgc", "TACGtacg")
    finalseq = seqreverse.translate(transtable)
    return finalseq


def get_motif(df, s, e, n, fa_file, pre_path, modification_type):
    # 创建副本以避免SettingWithCopyWarning
    df = df.copy()
    df["name"] = df["chrom"] + "_" + df["pos_1base"].astype(str) + "_" + df["strand"]
    df.loc[df["strand"] == "+", "Start"] = (
        df.loc[df["strand"] == "+", "pos_1base"] - s - 1
    )
    df.loc[df["strand"] == "+", "End"] = df.loc[df["strand"] == "+", "pos_1base"] + e
    df.loc[df["strand"] == "-", "Start"] = (
        df.loc[df["strand"] == "-", "pos_1base"] - e - 1
    )
    df.loc[df["strand"] == "-", "End"] = df.loc[df["strand"] == "-", "pos_1base"] + s
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)

    # 修复警告：明确创建新的DataFrame
    get_fasta = df[["chrom", "Start", "End", "name"]].copy()
    get_fasta["score"] = 0
    get_fasta["strand"] = df["strand"]

    get_fasta.to_csv(
        f"{pre_path}_{modification_type}_{n}.bed", header=None, sep="\t", index=False
    )
    call(
        f"bedtools getfasta -name -s -fi {fa_file} -bed {pre_path}_{modification_type}_{n}.bed -fo {pre_path}_{modification_type}_{n}.fa",
        shell=True,
    )
    return df


def check_base(pre_path, df, n, pos, expected_base, modification_type):
    # 创建副本以避免SettingWithCopyWarning
    df = df.copy()
    fasta_dict = {}
    with open(f"{pre_path}_{modification_type}_{n}.fa", "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_dict[record.id.split("::")[0]] = str(record.seq)
    df["motif_fa"] = df["name"].apply(lambda x: fasta_dict[x])
    df["base"] = df["motif_fa"].apply(lambda x: x[pos])
    df = df[(df["base"] == expected_base) | (df["base"] == expected_base.lower())]
    return df


def plot_motif_logo(df, figsize_width, plot_path):
    sequences = df["motif_fa"].str.upper().tolist()
    if len(sequences) == 0:
        print(f"Warning: No sequences found for motif plotting at {plot_path}!!!\n")
        return

    counts_df = logomaker.alignment_to_matrix(sequences, to_type="counts")
    freq_df = counts_df.div(counts_df.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(figsize_width, 4))
    logo = logomaker.Logo(freq_df, ax=ax)

    color_scheme = {"A": "red", "C": "blue", "G": "orange", "T": "green"}
    logo.style_glyphs(color_scheme=color_scheme)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_ylabel("Frequency")
    logo.ax.set_xlabel("Position")
    logo.ax.xaxis.set_ticks_position("none")
    logo.ax.yaxis.set_ticks_position("none")
    logo.ax.xaxis.set_tick_params(width=0)
    logo.ax.yaxis.set_tick_params(width=0)

    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(plot_path, format="pdf", dpi=600)
    plt.close()


def find_and_merge_files(output_prefix, modification_type):
    """查找并合并正负链文件"""
    sites_file_f = f"{output_prefix}_{modification_type}_sites_+_tmp.csv"
    sites_file_r = f"{output_prefix}_{modification_type}_sites_-_tmp.csv"
    reads_file_f = f"{output_prefix}_{modification_type}_reads_+_tmp.pkl"
    reads_file_r = f"{output_prefix}_{modification_type}_reads_-_tmp.pkl"

    # 基于 sites 文件判断存在哪些链；reads 文件缺失时仅告警并继续
    found_strands = []
    plus_sites_exists = os.path.exists(sites_file_f)
    minus_sites_exists = os.path.exists(sites_file_r)
    plus_reads_exists = os.path.exists(reads_file_f)
    minus_reads_exists = os.path.exists(reads_file_r)

    if plus_sites_exists:
        found_strands.append("+")
        if not plus_reads_exists:
            print(
                f"Warning: + strand reads file missing for {modification_type}, will proceed without reads.\n"
            )
    elif plus_reads_exists and not plus_sites_exists:
        print(
            f"Warning: + strand sites file missing for {modification_type}, + strand will be skipped for sites merge.\n"
        )

    if minus_sites_exists:
        found_strands.append("-")
        if not minus_reads_exists:
            print(
                f"Warning: - strand reads file missing for {modification_type}, will proceed without reads.\n"
            )
    elif minus_reads_exists and not minus_sites_exists:
        print(
            f"Warning: - strand sites file missing for {modification_type}, - strand will be skipped for sites merge.\n"
        )

    if not found_strands:
        detail = []
        for label, exists in (
            ("sites +", plus_sites_exists),
            ("reads +", plus_reads_exists),
            ("sites -", minus_sites_exists),
            ("reads -", minus_reads_exists),
        ):
            detail.append(f"{label}:{'Y' if exists else 'N'}")
        raise FileNotFoundError(
            f"No {modification_type} sites files found with prefix {output_prefix}. Presence -> "
            + ", ".join(detail)
        )

    if len(found_strands) == 1:
        print(
            f"Warning: Only {found_strands[0]} strand sites found for {modification_type}\n"
        )
    else:
        print(f"Found both + and - strand sites for {modification_type}")

    # 合并sites文件
    dfs_all = []
    if "+" in found_strands:
        df_f = pd.read_csv(sites_file_f)
        dfs_all.append(df_f)
    if "-" in found_strands:
        df_r = pd.read_csv(sites_file_r)
        dfs_all.append(df_r)

    df_merged = pd.concat(dfs_all, ignore_index=True)
    df_merged = df_merged.sort_values(by=["chrom", "pos_1base"])

    # 合并reads文件
    final_reads_dict = {}
    if plus_reads_exists:
        try:
            with open(reads_file_f, "rb") as f:
                dict_f = pickle.load(f)
            final_reads_dict.update(dict_f)
        except Exception as e:
            print(f"Warning: Failed to load + strand reads file: {e}")
    if minus_reads_exists:
        try:
            with open(reads_file_r, "rb") as f:
                dict_r = pickle.load(f)
            final_reads_dict.update(dict_r)
        except Exception as e:
            print(f"Warning: Failed to load - strand reads file: {e}")

    return df_merged, final_reads_dict, found_strands


def process_m6A(df, args):
    """处理m6A修饰"""
    s, e = 2, 2
    df = get_motif(df, s, e, "m6A", args.fa, args.output_prefix, args.modification_type)
    df = check_base(args.output_prefix, df, "m6A", 2, "A", args.modification_type)

    if args.motifPaint:
        motif_path = f"{args.output_prefix}_m6A_motif.pdf"
        plot_motif_logo(df, 5, motif_path)
        print(f"m6A motif plot saved to {motif_path}")

    return df[["chrom", "pos_1base", "strand", "mod_num", "cov", "mod_rate", "motif"]]


def process_m5C(df, args):
    """处理m5C修饰"""
    # 使用copy()确保我们操作的是副本
    df_CG = df[df["motif_classification"] == "CG"].copy()
    df_CHH = df[df["motif_classification"] == "CHH"].copy()
    df_CHG = df[df["motif_classification"] == "CHG"].copy()

    processed_dfs = []

    if not df_CG.empty:
        df_CG = get_motif(
            df_CG, 0, 1, "CG", args.fa, args.output_prefix, args.modification_type
        )
        df_CG = check_base(
            args.output_prefix, df_CG, "CG", 0, "C", args.modification_type
        )
        processed_dfs.append(df_CG)
    if not df_CHH.empty:
        df_CHH = get_motif(
            df_CHH, 0, 2, "CHH", args.fa, args.output_prefix, args.modification_type
        )
        df_CHH = check_base(
            args.output_prefix, df_CHH, "CHH", 0, "C", args.modification_type
        )
        processed_dfs.append(df_CHH)
    if not df_CHG.empty:
        df_CHG = get_motif(
            df_CHG, 0, 2, "CHG", args.fa, args.output_prefix, args.modification_type
        )
        df_CHG = check_base(
            args.output_prefix, df_CHG, "CHG", 0, "C", args.modification_type
        )
        processed_dfs.append(df_CHG)

    if processed_dfs:
        merged_df = pd.concat(processed_dfs, ignore_index=True)
        return merged_df[
            [
                "chrom",
                "pos_1base",
                "strand",
                "mod_num",
                "cov",
                "mod_rate",
                "motif",
                "motif_classification",
            ]
        ]
    else:
        return pd.DataFrame()


def process_pseU(df, args):
    # 使用motif_classification分组(更准确)
    if 'motif_classification' in df.columns:
        print("使用motif_classification列进行分组")
        df1 = df[df['motif_classification'] == 'pus1'].copy()
        df4 = df[df['motif_classification'] == 'pus4'].copy()
        df7 = df[df['motif_classification'] == 'pus7'].copy()
        
        print(f"  pus1: {len(df1)} 个位点")
        print(f"  pus4: {len(df4)} 个位点")
        print(f"  pus7: {len(df7)} 个位点")
    else:
        # 回退到按motif长度分组(兼容旧数据)
        print("警告: 没有motif_classification列,使用motif长度分组")
        df1 = df[df["motif"].apply(lambda x: len(x) == 3)].copy()
        df7 = df[df["motif"].apply(lambda x: len(x) == 5)].copy()
        df4 = df[df["motif"].apply(lambda x: len(x) == 6)].copy()
        
        print(f"  长度=3 (pus1): {len(df1)} 个位点")
        print(f"  长度=5 (pus7): {len(df7)} 个位点")
        print(f"  长度=6 (pus4): {len(df4)} 个位点")

    processed_dfs = []

    if not df1.empty:
        df1 = get_motif(
            df1, 2, 0, "pus1", args.fa, args.output_prefix, args.modification_type
        )
        df1 = check_base(
            args.output_prefix, df1, "pus1", 2, "T", args.modification_type
        )
        processed_dfs.append(df1)
        if args.motifPaint:
            pus1_path = f"{args.output_prefix}_pseU_motif_pus1.pdf"
            plot_motif_logo(df1, 3, pus1_path)
            print(f"PUS1 motif plot saved to {pus1_path}")
    else:
        print("No PUS1 motif found for pseU")

    if not df4.empty:
        df4 = get_motif(
            df4, 2, 3, "pus4", args.fa, args.output_prefix, args.modification_type
        )
        df4 = check_base(
            args.output_prefix, df4, "pus4", 2, "T", args.modification_type
        )
        processed_dfs.append(df4)
        if args.motifPaint:
            pus4_path = f"{args.output_prefix}_pseU_motif_pus4.pdf"
            plot_motif_logo(df4, 6, pus4_path)
            print(f"PUS4 motif plot saved to {pus4_path}")
    else:
        print("No PUS4 motif found for pseU")

    if not df7.empty:
        df7 = get_motif(
            df7, 2, 2, "pus7", args.fa, args.output_prefix, args.modification_type
        )
        df7 = check_base(
            args.output_prefix, df7, "pus7", 2, "T", args.modification_type
        )
        processed_dfs.append(df7)
        if args.motifPaint:
            pus7_path = f"{args.output_prefix}_pseU_motif_pus7.pdf"
            plot_motif_logo(df7, 5, pus7_path)
            print(f"PUS7 motif plot saved to {pus7_path}")
    else:
        print("No PUS7 motif found for pseU")

    if processed_dfs:
        merged_df = pd.concat(processed_dfs, ignore_index=True)
        return merged_df[
            ["chrom", "pos_1base", "strand", "mod_num", "cov", "mod_rate", "motif"]
        ]
    else:
        return pd.DataFrame()


def process_inosine(df, args):
    """处理inosine修饰"""
    s, e = 1, 1
    df = get_motif(
        df, s, e, "inosine", args.fa, args.output_prefix, args.modification_type
    )
    df = check_base(args.output_prefix, df, "inosine", 1, "A", args.modification_type)

    return df[["chrom", "pos_1base", "strand", "mod_num", "cov", "mod_rate", "motif"]]


def cleanup_tmp_files(output_prefix, modification_type):
    """删除所有临时文件"""
    # 定义需要删除的文件模式
    file_patterns = [
        f"*tmp*",  # 临时文件
        f"{output_prefix}_{modification_type}*.bed",  # bed文件
        f"{output_prefix}_{modification_type}*.fa",  # fasta文件
    ]

    all_tmp_files = []
    for pattern in file_patterns:
        all_tmp_files.extend(glob.glob(pattern))

    if all_tmp_files:
        for file in all_tmp_files:
            try:
                os.remove(file)
            except OSError:
                pass
        print(f"Cleaned up {len(all_tmp_files)} temporary files")


if __name__ == "__main__":
    # 抑制pandas警告（可选）
    warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

    parser = argparse.ArgumentParser(
        description="Unified script to merge modification results and optionally plot motif logos",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python mergeModSites.py -o output_prefix -m m6A -f reference.fa
  python mergeModSites.py -o output_prefix -m m5C -f reference.fa --motifPaint --metaPlotR
        """,
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        type=str,
        required=True,
        help="Output directory and prefix (same as used in previous steps)",
    )
    parser.add_argument(
        "-m",
        "--modification_type",
        type=str,
        required=True,
        choices=["m6A", "m5C", "pseU", "inosine"],
        help="Type of modification to process",
    )
    parser.add_argument(
        "-f", "--fa", type=str, required=True, help="Reference fasta file path"
    )
    parser.add_argument(
        "--motifPaint", action="store_true", help="Generate motif logo plots"
    )
    parser.add_argument(
        "--metaPlotR", action="store_true", help="Generate metaPlotR bed file"
    )

    args = parser.parse_args()

    try:
        # 查找并合并文件
        print(f"Processing {args.modification_type} modification...\n")
        df_merged, reads_dict, found_strands = find_and_merge_files(
            args.output_prefix, args.modification_type
        )

        # 根据修饰类型处理数据
        if args.modification_type == "m6A":
            df_final = process_m6A(df_merged, args)
        elif args.modification_type == "m5C":
            df_final = process_m5C(df_merged, args)
        elif args.modification_type == "pseU":
            df_final = process_pseU(df_merged, args)
        elif args.modification_type == "inosine":
            df_final = process_inosine(df_merged, args)

        # 保存结果
        if not df_final.empty:
            res_path = f"{args.output_prefix}_{args.modification_type}_sites_result.csv"
            df_final.to_csv(res_path, index=False)
            print(f"Results saved to {res_path}")

            # 过滤并保存reads字典
            df_final = df_final.copy()
            df_final["pos_0base"] = df_final["pos_1base"] - 1
            df_final["k"] = (
                df_final["chrom"]
                + "_"
                + df_final["pos_0base"].astype(str)
                + "_"
                + df_final["strand"]
            )
            readid_set = set(df_final["k"])
            filtered_reads_dict = {
                key: value for key, value in reads_dict.items() if key in readid_set
            }

            pkl_path = f"{args.output_prefix}_{args.modification_type}_read_final.pkl"
            with open(pkl_path, "wb") as file:
                pickle.dump(filtered_reads_dict, file)
            print(f"Filtered reads dictionary saved to {pkl_path}")

            # 清理临时文件
            cleanup_tmp_files(args.output_prefix, args.modification_type)
            print(
                f"\nProcessing completed for {args.output_prefix}_{args.modification_type}"
            )

            # 生成metaPlotR文件
            if args.metaPlotR:
                metaPlotR_path = (
                    f"{args.output_prefix}_{args.modification_type}_metaPlotR.bed"
                )
                df_meta = df_final.copy()

                df_meta["index"] = pd.Series(range(1, len(df_meta) + 1))
                df_meta = df_meta[
                    ["chrom", "pos_0base", "pos_1base", "index", "mod_rate", "strand"]
                ]
                df_meta.to_csv(metaPlotR_path, index=False, sep="\t", header=False)
                print(f"MetaPlotR file saved to {metaPlotR_path}")
        else:
            print("Warning: No valid sites found after filtering")

    except Exception as e:
        print(f"Error: {e}")
    call(f'rm {args.output_prefix}*tmp*', shell=True)
