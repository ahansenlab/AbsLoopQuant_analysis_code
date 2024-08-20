# Absolute Loop Quantification analysis code

This repository contains source code for James' absolute loop quantification project. (Eventually add the name of the manuscript in this introduction section.)

Code is provided as shell scripts, Python scripts, or Jupyter notebooks to be run in conda environments containing the required packages.


### Main scripts for Micro-C data processing
__________________

#### Alignment of Micro-C reads (microc_bwamem_with_recovery.py)

Given paired-end reads from a Micro-C experiment in .fastq format, this script aligns the reads, recovers multiply mapped reads in user-defined regions of interest, parses the reads into pairs representing genomic contacts, and removes optical/PCR duplicates. The output is given in .pairs format for downstream processing.

Example usage:
```
python microc_bwamem_with_recovery.py --files_1 sample_name_S1_L001_R1_001.fastq.gz,sample_name_S1_L002_R1_001.fastq.gz --files_2 sample_name_S1_L001_R2_001.fastq.gz,sample_name_S2_L002_R2_001.fastq.gz --genome mm39.fasta --assembly mm39 --rois chr15,11722732,11732141,chr15,11567241,11572268,chr1,34257542,34257653,chr4,132978078,132978190,chr8,13511989,13512092 --threads 8 --name sample_name --outdir output_directory
```

#### Process pairs to mcool (process_pairs_to_mcool.py)

This script takes deduplicated pairs in .pairs format, downsamples trans pairs to achieve a desired cis/trans ratio (if necessary), doubles read counts in heterozygous regions (if necessary), and produces an output .mcool file containing the balanced Micro-C contact matrix at various resolutions.

Example usage:
```
python process_pairs_to_mcool.py --name sample_name --genome mm39.fasta --assembly mm39 --threads 8 --transfraction 0.05 --ignoredist 20000 --outdir output_directory
```

### Helper scripts for Micro-C data processing
__________________

The following scripts are called by the two main scripts for Micro-C data processing described above.

#### Redistribution of multimapping reads (redistribute_multimapped_reads.py)

This script recovers reads within a user-defined region of interest that fail to map uniquely to the genome, and redistributes each one randomly among its possible positions. The input is a .pairs file (or a list of .pairs files) containing pairs where one or both sides are multimapping and possibly map to the region of interest. The output is a new .pairs file with new positions assigned to multimapping reads.

Example usage:
```
python redistribute_multimapped_reads.py --name sample_name --filenames sample_name_in_ROIs_1.pairs,sample_name_in_ROIs_2.pairs --rois chr15,11722732,11732141,chr15,11567241,11572268,chr1,34257542,34257653,chr4,132978078,132978190,chr8,13511989,13512092 --genome mm39_modified.fasta --min_mapq 30
```

#### Downsample trans pairs (downsample_trans_reads.py)

This script removes trans pairs randomly (independently with a fixed probability) from a .pairs file to reduce the fraction of pairs that are trans to a desired number. The second argument is the scaling factor, which is the fraction by which to reduce the number of trans reads.

Example usage:
```python downsample_trans_reads.py input.pairs 0.4 output.pairs```

#### Double read counts in TetO_LacO (double_read_counts_in_TetO_LacO_bins.py)

This script doubles the raw read counts in the heterozygous synthetic insertions ("TetO" & "LacO") of the TetO-LacO+3xCTCF (1B1) cell line in order to ensure fair comparison with homozygous regions. The code was written specifically to perform the operation on a raw .cool file with a bin size of 250 bp.

Example usage:
```python double_read_counts_in_TetO_LacO_bins.py TetO_LacO_rep1_250bp.cool TetO_LacO_rep1_syn_regions_doubled_250bp.cool```

### 3D polymer simulations
__________________

### Loop calling, quantification, and classification
__________________

#### Calculate P(s) curves (calculate_P_s_curves.py)

This script calculates chromosome-averaged P(s) curves at 1 kb resolution for the four Micro-C samples generated in this project, plus the ultra-deep merged Micro-C dataset which includes data from past studies ("all_merged").

#### Process and combine Mustache loops (process_combine_mustache_loops.ipynb)

This Jupyter notebook contains code to process the outputs from Mustache and merge the loops called at different resolutions (1kb, 2kb, and 5kb) into a single list of loops without duplicates. This code calls the helper script `merge_1kb_2kb_5kb_loops.sh` to perform the merging.

#### Filtering of quantifiable chromatin loops (filter_loops.py)

This Python script filters the chromatin loops based on various criteria for quantifiability.

#### Quantify loops

#### Process ChIP-seq data for loop classification

#### Classify loops by mechanism


