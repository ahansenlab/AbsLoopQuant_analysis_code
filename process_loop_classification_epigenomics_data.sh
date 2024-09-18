# copy promoters file to data folder
cp /mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/loop_classification/mm39_TSS_pad_2kb.bed /mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data

# identify enhancers
cd /mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/
bedtools intersect -a H3K4me1_GSE90893/H3K4me1_GSE90893_peaks.broadPeak -b H3K27ac_GSE90893/H3K27ac_GSE90893_broad_peaks.broadPeak > loop_classification/H3K4me1_H3K27ac.bed

cd /mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/loop_classification
bedtools subtract -a H3K4me1_H3K27ac.bed -b mm39_TSS_pad_2kb.bed > /mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/mm39_enhancers.bed

# identify CTCF sites
cd /mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/
bedtools intersect -a /mnt/coldstorage/shares/Miles/Analysis/Motifs/mm39/CTCFmm39loose/CTCFmm39.1e3.sites.bed -b CTCF_GSE90994/CTCF_GSE90994_peaks.narrowPeak -wb -f 1 > loop_classification/CTCF_sites_CTCF.peaks

cd /mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/loop_classification
Rscript /mnt/md0/Tools/scripts/GetUniqueMotifsFromFIMObedWithPeaks.R -b CTCF_sites_CTCF.peaks -o CTCF_sites_CTCF_unique_peaks.bed -t  # get the single best peak within each CTCF 

# get CTCF and cohesin bound sites, put result in data folder
cd /mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/
bedtools intersect -a loop_classification/CTCF_sites_CTCF_unique_peaks.bed -b SMC1A_GSE123636/SMC1A_GSE123636_peaks.narrowPeak -u > /mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/mm39_CTCF_cohesin.bed