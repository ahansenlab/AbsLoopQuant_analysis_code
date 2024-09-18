import numpy as np
import pandas as pd
import cooler
import time
from multiprocess import Pool

## setup ##

# import the loop quantification module
import sys
sys.path.insert(1, '/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_analysis_code')
import looptools

# load the list of filtered loops
filtered_loops_df = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/filtered_loops.csv', index_col=0)

# load the coolers & P(s) curves
coolers = {}
res = 1000
coolers['T1'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/T1_final/T1.mcool::/resolutions/{res}')
coolers['T2'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/T2_final/T2.mcool::/resolutions/{res}')
coolers['C3'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/C3_final/C3.mcool::/resolutions/{res}')
coolers['C4'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/C4_final/C4.mcool::/resolutions/{res}')
coolers['all_merged'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_all_merged.mcool::/resolutions/{res}')

P_s_curves = {}
for sample_name in coolers.keys():
    P_s_curves[sample_name] = np.loadtxt(f'/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/P_s_{sample_name}_1000bp.txt')

# functions to transform coordinates from WT genomes to synthetic genomes
def to_C_coords(chrom, pos):
    '''Convert WT mm39 coordinates to the proper coordinates for the C36 modified genome.'''
    if chrom=='chr18':
        if pos >= 58619130:
            return pos+14525
        elif pos >= 58104202:
            return pos+12081
        else:
            return pos
    else:
        return pos
    
def to_T_coords(chrom, pos):
    '''Convert WT mm39 coordinates to the proper coordinates for the TetO-LacO modified genome.'''
    if chrom=='chr15':
        if pos >= 11717705:
            return pos+14438
        elif pos >= 11567241:
            return pos+5029
        else:
            return pos
    else:
        return pos
    
## define function to perform loop quantification ##
    
def quantify_loop_and_save_result(i):
    
    chr_name, start, end, size = filtered_loops_df.loc[i]
    local_region_size = int(np.round(np.sqrt(size/1000/(32/25**2)))) * 1000  # in bp
    if convert_coords:
        start = coords_convert_function(chr_name, start)
        end = coords_convert_function(chr_name, end)
    score = lq.quantify_loop(chr_name, start, end,
                                 local_region_size=local_region_size,
                                 quant_region_size=10000,
                                 k_min=2,
                                 clr_for_outlier_detection=coolers['all_merged'], P_s_values_for_outlier_detection=P_s_curves['all_merged'])
    return score

## set up the multiprocessing ##

nproc = 40  # number of threads to use for multiprocessing
chunk_size = 40

num_chunks = int(np.ceil(len(filtered_loops_df)/chunk_size))
chunk_starts = np.arange(num_chunks)*chunk_size
chunk_ends = (np.arange(num_chunks)+1)*chunk_size
chunk_ends[-1] = len(filtered_loops_df)

## perform the loop quantification ##

# initialize dataframe to store results
loop_quant_df = pd.DataFrame(index=filtered_loops_df.index, columns=['T1','T2','C3','C4'])

for sample_name in ['T1','T2','C3','C4']:

    clr = coolers[sample_name]
    P_s_values = P_s_curves[sample_name]
    lq = looptools.LoopQuantifier(clr, P_s_values)
    if sample_name.startswith('T'):
        convert_coords = True
        coords_convert_function = to_T_coords
    elif sample_name.startswith('C'):
        convert_coords = True
        coords_convert_function = to_C_coords
    else:
        convert_coords = False
    
    for chunk_index in np.arange(num_chunks):
        start = time.time()
        with Pool(nproc) as p:
            indices_in_chunk = filtered_loops_df.index[np.arange(chunk_starts[chunk_index],chunk_ends[chunk_index])]
            scores_in_chunk = p.map(quantify_loop_and_save_result, indices_in_chunk)
            loop_quant_df.loc[indices_in_chunk,sample_name] = np.array(scores_in_chunk)
        end = time.time()
        print(f'{sample_name}, chunk {chunk_index+1} of {num_chunks} ({end-start:.2f} s)')
        if chunk_index % 50 == 0:
            loop_quant_df.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/MicroC_loop_quantification_scores.csv', float_format='%.6f')
    loop_quant_df.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/MicroC_loop_quantification_scores.csv', float_format='%.6f')