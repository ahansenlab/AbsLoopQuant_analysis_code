import argparse
import subprocess as sp
import os

parser = argparse.ArgumentParser(description = 'align Micro-C data with bwa-mem2, recover multimapped reads in user-defined regions of interest, and generate deduplicated pairs files')
parser.add_argument('--files_1', '-1', help='comma-separated list of fastq files containing left mates of paired end reads')
parser.add_argument('--files_2', '-2', help='comma-separated list of fastq files containing left mates of paired end reads')
parser.add_argument('--genome', '-g', help='path to genome file to align to (also need chromsizes file with the suffix .sorted.chrom.sizes in the same directory)')
parser.add_argument('--assembly', '-a', help='genome assembly name')
parser.add_argument("--rois", "-r", help="regions of interest that could contain the primary hits of the multimapping reads, formatted as a comma-separated list of chr,start,end, e.g., 'chr15,11722732,11732141,chr15,11567241,11572268,chr1,34257542,34257653,chr4,132978078,132978190,chr8,13511989,13512092'")
parser.add_argument('--threads', '-t', help='number of threads to use', default='1')
parser.add_argument('--name', '-n', help='sample name prefix to use for output files')
parser.add_argument('--outdir', '-o', help='output directory')

args = parser.parse_args()

files1 = args.files_1.split(',')
files2 = args.files_2.split(',')
genome_path = args.genome
assembly = args.assembly
ROI_string = args.rois
nproc = args.threads
sample_name = args.name
outdir = args.outdir

scripts_path = os.path.dirname(os.path.realpath(__file__))  # path to directory of current file (where all the scripts are contained)

# create temporary directory for sorting and merging pairs files
tempdir = os.path.join(outdir, 'tempdir')
sp.run(f'mkdir -p {tempdir}', shell=True)

# parse the ROI string
regions_of_interest = [(ROI_string.split(',')[i], int(ROI_string.split(',')[i+1]), int(ROI_string.split(',')[i+2])) for i in range(0,len(ROI_string.split(',')),3)]
region_conditions_1 = [f'((mapq1<30) and (chrom1=="{r[0]}") and (pos1>={r[1]}) and (pos1<={r[2]}))' for r in regions_of_interest]
region_conditions_2 = [f'((mapq2<30) and (chrom2=="{r[0]}") and (pos2>={r[1]}) and (pos2<={r[2]}))' for r in regions_of_interest]
region_conditions_all = region_conditions_1 + region_conditions_2
condition = ' or '.join(region_conditions_all)  # to be passed to pairtools select

# alignment
for i in range(len(files1)):
    file1 = files1[i]
    file2 = files2[i]
    sp.run(f'''bwa-mem2 mem -t {nproc} -SP {genome_path} {file1} {file2} | pairtools parse --add-columns mapq --walks-policy all -c {genome_path}.sorted.chrom.sizes --assembly {assembly} --min-mapq 30 --drop-sam --drop-seq --drop-readid --nproc-in {nproc} | pairtools sort --tmpdir {tempdir} --nproc {nproc} -o {os.path.join(outdir, sample_name)}_{i+1}.pairs.gz''', shell=True)
    sp.run(f"""bwa-mem2 mem -t {nproc} -SP {genome_path} {file1} {file2} | pairtools parse --add-columns mapq,algn_ref_span --walks-policy all -c {genome_path}.sorted.chrom.sizes --assembly {assembly} --min-mapq 0 --drop-seq --drop-readid --nproc-in {nproc} | pairtools select -o {os.path.join(outdir, sample_name)}_in_selected_regions_{i+1}.pairs '{condition}'""", shell=True)

# redistribution of multimapped reads
sp.run(f'python {os.path.join(scripts_path, "redistribute_multimapped_reads.py")} --name {sample_name} --filenames {",".join([f"{sample_name}_in_selected_regions_{i+1}.pairs" for i in range(len(files1))])} --rois {ROI_string} --genome {genome_path} --min_mapq 30', shell=True, cwd=outdir)

# fix header and sort recovered reads
sp.run(f'pairtools header transfer -r {sample_name}_1.pairs.gz -o {sample_name}_in_selected_regions_redist_header.pairs {sample_name}_in_selected_regions_redist.pairs', shell=True, cwd=outdir)
sp.run(f'pairtools sort {sample_name}_in_selected_regions_redist_header.pairs -o {sample_name}_in_selected_regions_redist_sorted.pairs', shell=True, cwd=outdir)

# merge and deduplicate pairs files
sp.run(f'''pairtools merge --tmpdir tempdir --nproc {nproc} {' '.join([f"{os.path.join(outdir, sample_name)}_{i+1}.pairs.gz" for i in range(len(files1))])} {sample_name}_in_selected_regions_redist_sorted.pairs | pairtools dedup --max-mismatch 1 --mark-dups --output {sample_name}_nodups.pairs.gz --output-unmapped {sample_name}_unmapped.pairs.gz --output-dups {sample_name}_dups.pairs.gz --output-stats {sample_name}_dedup.stats''', shell=True, cwd=outdir)
