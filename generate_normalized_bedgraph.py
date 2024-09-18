import subprocess as sp
import sys

merged_bam_filename = sys.argv[1]
output_filename = sys.argv[2]

# count number of lines
result = sp.run(f"samtools idxstats {merged_bam_filename} | awk '{{sum+=$3}}END{{print sum}}'", shell=True, capture_output=True)
num_reads = int(result.stdout)
scale = 1e6/num_reads
print(f'Counted {num_reads} reads in input bam. Normalizing output bedgraph by scale factor {scale:.3}.')

# generate normalized bedgraph
sp.run(f'bedtools genomecov -ibam  {merged_bam_filename} -bga -scale {scale:.3} > {output_filename}', shell=True)