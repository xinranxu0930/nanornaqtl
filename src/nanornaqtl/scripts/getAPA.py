import pysam
from subprocess import call
import argparse
import pandas as pd
import sys

def get_apa(bamfile, strand):
    res_dict = set()
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for read in samfile:
        chrom = read.reference_name
        if strand == '+':
            s = read.reference_end -1
            e = read.reference_end +1
        else:  # strand == '-'
            s = read.reference_start
            e = read.reference_start+2
        res_dict.add((chrom, s, e, read.query_name, 0, strand))
    return res_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Get APA for every read from Mapped BAM file",
    formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Data Format Description:
  PAS Sites File from APAdb is a BED format file with the following format and content:
  chr1	16442	16450	WASH7P.1:16442-16450	43	-	Intron	16443.0	16443	-
  chr1	134934	134953	LOC729737.1:134934-134953	26	-	UTR3	134944.0	134944	-
  """
  )
    parser.add_argument("-b","--bam", type=str, required=True, help="Mapped BAM file")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-a","--apa_file", type=str, help="APA File from APAdb")
    args = parser.parse_args()
    try:
        bam0 = args.bam.replace(".bam", "0.bam")
        bam16 = args.bam.replace(".bam", "16.bam")

        res_dict_0 = get_apa(bam0, "+")
        res_dict_16 = get_apa(bam16, "-")

        merged_set = res_dict_0.union(res_dict_16)

        # 将集合内容写入BED文件
        with open(f'{args.output_prefix}_readend.bed', 'w') as file:
            for item in merged_set:
                file.write('\t'.join(map(str, item)) + '\n')

        call(f'bedtools intersect -s -a {args.output_prefix}_readend.bed -b {args.apa_file} -wa -wb> {args.output_prefix}_overlap.bed',shell=True)
        apa_sites = pd.read_csv(f'{args.output_prefix}_overlap.bed', sep='\t', header=None, usecols=[3,9])
        apa_sites.columns = ["readID","APA_type"]
        apa_sites.to_csv(f'{args.output_prefix}_APA_result.csv', index=False)
        call(f'rm {args.output_prefix}_readend.bed {args.output_prefix}_overlap.bed',shell=True)
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)