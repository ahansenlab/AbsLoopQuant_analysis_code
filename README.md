# Absolute Loop Quantification analysis code

This repository contains source code for James' absolute loop quantification project. (Eventually add the name of the manuscript in this introduction section.)

Code is provided as shell scripts, Python scripts, or Jupyter notebooks to be run in conda environments containing the required packages.


### Main scripts for Micro-C data processing
__________________

#### Alignment of Micro-C reads (microc_bwamem_with_recovery.py)

Given paired-end reads from a Micro-C experiment in .fastq format, this script aligns the reads, recovers multiply mapped reads in user-defined regions of interest, parses the reads into pairs representing genomic contacts, and removes optical/PCR duplicates. The output is given in .pairs format for downstream processing.

Example usage:
```python process_pairs_to_mcool.py --name sample_name --genome mm39.fasta --assembly mm39 --threads 8 --transfraction 0.05 --ignoredist 20000 --outdir output_directory```

#### Process pairs to mcool (process_pairs_to_mcool.py)

This script takes deduplicated pairs in .pairs format, downsamples trans reads to achieve a desired cis/trans ratio (if necessary), doubles read counts in heterozygous regions (if necessary), and produces an output .mcool file containing the balanced Micro-C contact matrix at various resolutions.

### Helper scripts for Micro-C data processing
__________________

The following scripts are called by the two main scripts for Micro-C data processing described above.

#### Redistribution of multimapping reads (redistribute_multimapped_reads.py)

This script recovers reads within a user-defined region of interest that fail to map uniquely to the genome, and redistributes each one randomly among its possible positions.

#### Double read counts in TetO_LacO (double_read_counts_in_TetO_LacO_bins.py)

This script doubles the read counts in 