import numpy as np
import pandas as pd
import pyBigWig
import sys

loops_df = pd.read_csv('../AbsLoopQuant_data/filtered_loops.csv', index_col=0)

dataset_name = sys.argv[1]
bw = pyBigWig.open(f'/mnt/coldstorage/jjusuf/Past_Epigenomics_Experiments/{dataset_name}/{dataset_name}_normalized.bw')

pad = 2500

def get_left_anchor_signal(i):
    chrom = loops_df.iloc[i]['chr']
    coord = loops_df.iloc[i]['left']
    signal_arr = bw.values(chrom, coord-pad, coord+pad)
    return np.nansum(signal_arr)

def get_right_anchor_signal(i):
    chrom = loops_df.iloc[i]['chr']
    coord = loops_df.iloc[i]['right']
    signal_arr = bw.values(chrom, coord-pad, coord+pad)
    return np.nansum(signal_arr)

anchor_signal_df = np.zeros((len(loops_df), 2))

for i in range(len(loops_df)):
    anchor_signal_df[i, 0] = get_left_anchor_signal(i)
    anchor_signal_df[i, 1] = get_right_anchor_signal(i)
    
    if i % 1000 == 0:
        np.savetxt(f'../AbsLoopQuant_data/bigwig_signal_at_anchors/{dataset_name}_normalized_signal_at_anchors_pad_2.5kb.txt', anchor_signal_df, fmt='%.4e')
np.savetxt(f'../AbsLoopQuant_data/bigwig_signal_at_anchors/{dataset_name}_normalized_signal_at_anchors_pad_2.5kb.txt', anchor_signal_df, fmt='%.4e')