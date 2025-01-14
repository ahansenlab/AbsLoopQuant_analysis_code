# Absolute Loop Quantification analysis code

This repository contains source code used for analysis in the preprint _Genome-wide absolute quantification of chromatin looping_ by James M. Jusuf, Simon Grosse-Holz, Michele Gabriele, Pia Mach, Ilya M. Flyamer, Christoph Zechner, Luca Giorgetti, Leonid A. Mirny, and Anders S. Hansen (2025).

Most of the data needed to run this code is in the folder `AbsLoopQuant_data` and will be shared to Zenodo upon final publication. The raw and processed Micro-C data (`.mcool` files) are available on NCBI GEO (accession number GSE286495).

Code is provided as shell scripts, Python scripts, or Jupyter notebooks to be run in conda environments.

## Loop quantification demo

The Jupyter notebook `loop_quantification_demo.ipynb` contains a short demo of quantifying a loop using AbLE (Absolute Loop Estimator) in one of our Micro-C maps. The code uses the _Npr3_ replicate 1 and merged mESC Micro-C maps, which must be downloaded from GEO (GSE286495); approximately 36 GB total.

The `looptools.py` module should be placed in the working directory and imported as shown in the first cell of the demo; there are no additional installation steps. The following standard packages are required:
* numpy (1.24.4)
* scipy (1.11.2)
* pandas (1.5.3)
* matplotlib (3.8.0)
* opencv-python (4.10.0.82)
* cooler (0.9.3)
* cooltools (0.5.4)

The software has been tested on the version numbers indicated in the parentheses on a computer running Ubuntu 20.04. The demo should take a few seconds to run on a typical desktop computer.

## Main scripts for Micro-C data processing

#### Alignment of Micro-C reads (microc_bwamem_with_recovery.py)

Given paired-end reads from a Micro-C experiment in .fastq format, this Python script aligns the reads, recovers multiply mapped reads in user-defined regions of interest, parses the reads into pairs representing genomic contacts, and removes optical/PCR duplicates. The output is given in .pairs format.

Example usage:
```
python microc_bwamem_with_recovery.py --files_1 sample_name_S1_L001_R1_001.fastq.gz,sample_name_S1_L002_R1_001.fastq.gz --files_2 sample_name_S1_L001_R2_001.fastq.gz,sample_name_S2_L002_R2_001.fastq.gz --genome mm39.fasta --assembly mm39 --rois chr15,11722732,11732141,chr15,11567241,11572268,chr1,34257542,34257653,chr4,132978078,132978190,chr8,13511989,13512092 --threads 8 --name sample_name --outdir output_directory
```

#### Process pairs to mcool (process_pairs_to_mcool.py)

This Python script takes deduplicated pairs in .pairs format, downsamples trans pairs to achieve a desired cis/trans ratio (if necessary), doubles read counts in heterozygous regions (if necessary), and produces an output .mcool file containing the balanced Micro-C contact matrix at various resolutions.

Example usage:
```
python process_pairs_to_mcool.py --name sample_name --genome mm39.fasta --assembly mm39 --threads 8 --transfraction 0.05 --ignoredist 20000 --outdir output_directory
```

## Helper scripts for Micro-C data processing

The following scripts are called by the two main scripts for Micro-C data processing described above.

#### Redistribution of multimapping reads (redistribute_multimapped_reads.py)

This Python script recovers reads within a user-defined region of interest that fail to map uniquely to the genome, and redistributes each one randomly among its possible positions. The input is a .pairs file (or a list of .pairs files) containing pairs where one or both sides are multimapping and possibly map to the region of interest. The output is a new .pairs file with new positions assigned to multimapping reads.

Example usage:
```
python redistribute_multimapped_reads.py --name sample_name --filenames sample_name_in_ROIs_1.pairs,sample_name_in_ROIs_2.pairs --rois chr15,11722732,11732141,chr15,11567241,11572268,chr1,34257542,34257653,chr4,132978078,132978190,chr8,13511989,13512092 --genome mm39_modified.fasta --min_mapq 30
```

#### Downsample trans pairs (downsample_trans_reads.py)

This Python script removes trans pairs randomly (independently with a fixed probability) from a .pairs file to reduce the fraction of pairs that are trans to a desired number. The second argument is the scaling factor, which is the fraction by which to reduce the number of trans reads.

Example usage:
```python downsample_trans_reads.py input.pairs 0.4 output.pairs```

#### Double read counts in TetO_LacO (double_read_counts_in_TetO_LacO_bins.py)

This Python script doubles the raw read counts in the heterozygous synthetic insertions ("TetO" & "LacO") of the TetO-LacO+3xCTCF (1B1) cell line in order to ensure fair comparison with homozygous regions. The code was written specifically to perform the operation on a raw .cool file with a bin size of 250 bp.

Example usage:
```python double_read_counts_in_TetO_LacO_bins.py TetO_LacO_rep1_250bp.cool TetO_LacO_rep1_syn_regions_doubled_250bp.cool```


## Loop calling, quantification, and classification

#### Calculate P(s) curves (calculate_P_s_curves.py)

This Python script calculates chromosome-averaged P(s) curves at 1 kb resolution for the four Micro-C samples generated in this project, plus the ultra-deep merged Micro-C dataset which includes data from past studies (`all_merged`).

#### Process and combine Mustache loops (process_combine_mustache_loops.ipynb)

This Jupyter notebook contains code to process the outputs from Mustache and merge the loops called at different resolutions (1kb, 2kb, and 5kb) into a single list of loops without duplicates. This code calls the helper script `merge_1kb_2kb_5kb_loops.sh` to perform the merging.

#### Helper scripts for combining loops (merge_1kb_2kb_5kb_loops.sh and loopComparer.py)

These scripts are called by `process_combine_mustache_loops.ipynb` to merge lists of loops called at different resolutions.

#### Filtering of quantifiable chromatin loops (filter_loops.py)

This Python script filters the chromatin loops based on various criteria for quantifiability, saving the output as a .csv file. This script uses multiprocessing to speed up the operations.

#### Quantify loops (quantify_loops.py)

This Python script calculates the AbLE scores of the filtered loops, saving the output as a .csv file. This script uses multiprocessing to speed up the operations.

#### Absolute Looping Estimator (AbLE) module (looptools.py)

This is the main Python module for loop quantification by AbLE. It contains the functions to calculate the average P(s) curves (used in `calculate_P_s_curves.py`) and to calculate the AbLE scores of loops (used in `quantify_loops.py`).

#### Process ChIP-seq data for loop classification (process_loop_classification_epigenomics_data.sh)

This shell script processes the ChIP peak data (H3K4me1, H3K27ac, CTCF, and SMC1A), CTCF motif data, and TSS data, which is subsequently used to identify the locations of promoters, enhancers, and CTCF/cohesin-bound sites across the genome.

#### Classify loops by mechanism (classify_loops.py)

This Python script identifies whether the left and right anchor of each loop is a promoter, enhancer, or CTCF/cohesin-bound site, and classifies the loops as enhancer-promoter, promoter-promoter, cis-regulatory, CTCF-CTCF, mixed, etc. A final .csv table containing the loop positions, absolute looping probabilities, and loop classification details is generated.


## Analysis of epigenomic features at loop anchors

#### Alignment scripts (spikeinChIP_SE_alignment.py/spikeinChIP_PE_alignment.py)

These Python scripts were used to align raw reads (`SE` for single-end, `PE` for paired-end) from publicly available epigenomics experiments, mainly ChIP-seq experiments.

#### Generate bedGraph of pileups (generate_normalized_bedgraph.py)

Given a merged .bam file containing the alignments of an epigenomics experiment, this Python script generates a bedgraph file of the pileups across all genomic positions, normalized to reads per million.

Example usage:
```
python generate_normalized_bedgraph.py CTCF_GSE90994_merged_across_reps.bam CTCF_GSE90994.bedgraph
```

The resulting bedGraph can then be converted to a bigWig using the UCSC bedGraphToBigWig executable.

#### Calculate epigenomic feature signal at loop anchors (calculate_bigwig_signal_at_anchors.py)

This Python script calculates the signal of a given epigenomic feature (saved as a bigWig file) at the left and right anchors of every loop in the list of filtered loops. The output is a numpy array containing two columns (one for the left anchor signal and one for the right anchor signal), saved as a .txt file.

Example usage:
```python calculate_bigwig_signal_at_anchors.py CTCF_GSE90994```


## 3D polymer simulations

#### Running polymer simulations

The folder `3D_polysim_code` contains the code used to perform the polymer simulations themselves. There are three different conditions:

* `3D_polysim_with_loopextr_with_EP_3kBT.py`: with loop extrusion and enhancer-promoter attraction (generates main simulation data)
* `3D_polysim_with_loopextr_no_EP.py`: with loop extrusion, no enhancer-promoter attraction (data is used to estimate the ground-truth background interaction probability of enhancer-promoter contacts)
* `3D_polysim_no_loopextr_no_EP.py`: no loop extrusion, no enhancer-promoter attraction (data is used for calibration of length and time scales)
* `DSB_smcTranslocator_v2.pyx`: a Cython file that performs the 1D simulations of loop extrusion used to generate the bonds in the main 3D simulation

#### Processing simulated Micro-C map (simulation_data_processing.ipynb)

This Jupyter notebook contains the code to generate a simulated Micro-C contact map of the entire chromosome, store it in .cool format, and balance the contact map in the same fashion as for real experimental data.

#### Calculate looping probabilities and Micro-C dot strengths from simulations (calculate_simulation_calibration_curve.ipynb)

This Jupyter notebook is used to calculate the ground-truth looping probabilities (y-values) from the simulated trajectories, as well as the Micro-C dot strengths (x-values) from the simulated Micro-C contact map, to be used in the simulated calibration curve.


## How to cite

This work is shared under an MIT license. If you make use of analysis scripts or data from this work, please cite as follows:

Jusuf, J.M., Grosse-Holz, S., Gabriele, M., Mach, P., Flyamer, I.M., Zechner, C., Giorgetti, L., Mirny, Leonid, L.A., Hansen, A.S., Genome-wide absolute quantification of chromatin looping.