import argparse
import subprocess as sp
import warnings
import os

parser = argparse.ArgumentParser(description = 'from deduplicated pairs files, downsample trans reads (if necessary), double synthetic regions (if necessary) and generate mcool file')
parser.add_argument('--name', '-n', help='sample name prefix (input file should be <name>_nodups.pairs.gz in <outdir>)')
parser.add_argument('--genome', '-g', help='path to genome file that data was aligned to (need chromsizes file with the suffix .sorted.chrom.sizes in the same directory)')
parser.add_argument('--assembly', '-a', help='genome assembly name')
parser.add_argument('--threads', '-t', help='number of threads to use', default='1')
parser.add_argument('--transfraction', '-p', help='target fraction of trans reads for trans read downsampling, if necessary')
parser.add_argument('--double', '-d', help='double synthetic regions corresponding to TetO+3xCTCF and 3xCTCF+LacO regions in TetO+LacO+3xCTCF cell line from Mach et al. 2022', action='store_true')
parser.add_argument('--ignoredist', '-i', help='distance in bp to ignore from the diagonal when balancing contact map (only used at resolutions where it corresponds to more than 2 diagonals)')
parser.add_argument('--outdir', '-o', help='output directory')
parser.add_argument('--binsize', '-b', help='smallest cooler binsize to make', default='250')

args = parser.parse_args()
sample_name = args.name
genome_path = args.genome
assembly = args.assembly
nproc = args.threads
target_trans_fraction = float(args.transfraction)
double_syn_regions = args.double
ignore_dist = args.ignoredist
outdir = args.outdir
binsize = int(args.binsize)

scripts_path = os.path.dirname(os.path.realpath(__file__))  # path to directory of current file (where all the scripts are contained)

if target_trans_fraction is None:
    skip_trans_downsampling = True
else:
    skip_trans_downsampling = False
    
    # calculate the current fraction of trans reads
    dedup_stats_file = open(os.path.join(outdir, f'{sample_name}_dedup.stats'))
    dedup_stats_str = dedup_stats_file.read()
    dedup_stats_lines = dedup_stats_str.split('\n')
    for line in dedup_stats_lines:
        if line.split('\t')[0]=='cis':
            num_cis_reads = int(line.split('\t')[1])
        if line.split('\t')[0]=='trans':
            num_trans_reads = int(line.split('\t')[1])
    current_trans_fraction = num_trans_reads/(num_cis_reads+num_trans_reads)

    # calculate the necessary scaling factor by which to downsample the trans reads
    scaling_factor = target_trans_fraction/current_trans_fraction * (1-current_trans_fraction)/(1-target_trans_fraction)
    if scaling_factor > 1:
        warnings.warn(f'The target trans fraction ({target_trans_fraction:.4f}) is greater than the current trans fraction ({current_trans_fraction:.4f}); skipping trans downsampling.', stacklevel=2)
        skip_trans_downsampling = True

if skip_trans_downsampling:
    pairs_filename_for_cool = f'{sample_name}_nodups.pairs.gz'
else:
    # perform downsampling of trans reads
    sp.run(f'python {os.path.join(scripts_path, "downsample_trans_reads.py")} {sample_name}_nodups.pairs.gz {scaling_factor:.2f} {sample_name}_nodups_trans_downsampled_{scaling_factor:.2f}.pairs.gz', shell=True, cwd=outdir)

    pairs_filename_for_cool = f'{sample_name}_nodups_trans_downsampled_{scaling_factor:.2f}.pairs.gz'

# index pairs file
sp.run(f'pairix {pairs_filename_for_cool}', shell=True, cwd=outdir)

# generate cooler
sp.run(f'bgzip -cd -@ {nproc} {pairs_filename_for_cool} | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {assembly} {genome_path}.sorted.chrom.sizes:{binsize} - {sample_name}_{binsize}.cool', shell=True, cwd=outdir)

if double_syn_regions:
    # double synthetic regions
    sp.run(f'python {os.path.join(scripts_path, "double_read_counts_in_TetO_LacO_bins.py")} {sample_name}_{binsize}.cool {sample_name}_doubled_syn_regions.cool', shell=True, cwd=outdir)
    cool_filename_for_mcool = f'{sample_name}_doubled_syn_regions.cool'
else:
    cool_filename_for_mcool = f'{sample_name}_{binsize}.cool'

# calculate bin sizes
all_bin_sizes_list = ','.split('10000000,5000000,2000000,1000000,500000,200000,100000,50000,20000,10000,5000,2000,1000,500,250')
all_bin_sizes_list = [i for i in all_bin_sizes_list if i >= binsize]
all_bin_sizes = ','.join(all_bin_sizes_list)

# make mcool
sp.run(f'''cooler zoomify --nproc {nproc} --balance --balance-args '--nproc {nproc} --ignore-diags 2 --ignore-dist {ignore_dist}' --out {sample_name}.mcool --resolutions {allPlease move or remove them before you merge._bin_sizes} {cool_filename_for_mcool}''', shell=True, cwd=outdir)

# remove temporary directory
sp.run('rm -r tempdir', shell=True, cwd=outdir)
    
    
    
    
    
    
    
    