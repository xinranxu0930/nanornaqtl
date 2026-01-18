import argparse
import pysam

def res_0(bam_file,threads,map_bam_file0,min_mapq):
    with pysam.AlignmentFile(bam_file, "rb",threads=threads) as samfile:
        with pysam.AlignmentFile(map_bam_file0, "wb", header=samfile.header) as out_file:
            for read in samfile:
                if (read.flag == 0) and (read.mapping_quality >= min_mapq):
                    out_file.write(read)
def res_16(bam_file,threads,map_bam_file16,min_mapq):
    with pysam.AlignmentFile(bam_file, "rb",threads=threads) as samfile:
        with pysam.AlignmentFile(map_bam_file16, "wb", header=samfile.header) as out_file:
            for read in samfile:
                if (read.flag == 16) and (read.mapping_quality >= min_mapq):
                    out_file.write(read)
def res_all(bam_file,threads,map_bam_file,min_mapq):
    with pysam.AlignmentFile(bam_file, "rb",threads=threads) as samfile:
        with pysam.AlignmentFile(map_bam_file, "wb", header=samfile.header) as out_file:
            for read in samfile:
                if (read.flag == 0 or read.flag == 16) and (read.mapping_quality >= min_mapq):
                    out_file.write(read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--basecall_bam", type=str, help="Dorado basecall bam file")
    parser.add_argument("-p", "--dir_pre", type=str, help="Result dir pre\nUltimately, three files will be generated : *_calls_sorted_map0.bam *_calls_sorted_map16.bam *_calls_sorted_map.bam")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Threads (default 4)")
    parser.add_argument("-q", "--min_mapq", type=int, default=0, help="Minimum Read MAPQ score (default 0)")
    args = parser.parse_args()

    map_bam_file0 = f"{args.dir_pre}_calls_sorted_map0.bam"
    map_bam_file16 = f"{args.dir_pre}_calls_sorted_map16.bam"
    map_bam_file = f"{args.dir_pre}_calls_sorted_map.bam"

    res_16(args.basecall_bam,args.threads,map_bam_file16,args.min_mapq)
    res_0(args.basecall_bam,args.threads,map_bam_file0,args.min_mapq)
    res_all(args.basecall_bam,args.threads,map_bam_file,args.min_mapq)


