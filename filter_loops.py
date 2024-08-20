import numpy as np
import pandas as pd
import subprocess as sp
import time
from multiprocessing import Pool
import cv2
import cooler

# load coolers
coolers = {}
res = 1000
coolers['T1'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/T1_final/T1.mcool::/resolutions/{res}')
coolers['T2'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/T2_final/T2.mcool::/resolutions/{res}')
coolers['C3'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/C3_final/C3.mcool::/resolutions/{res}')
coolers['C4'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/C4_final/C4.mcool::/resolutions/{res}')
coolers['all_merged'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_all_merged.mcool::/resolutions/{res}')

# load P(s) curves for each sample
P_s_curves = {}
for sample_name in coolers.keys():
    P_s_curves[sample_name] = np.loadtxt(f'/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/P_s_{sample_name}_1000bp.txt')
    
# load table of loops
loops_bedpe = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/mESC_all_merged_loops_Mustache_final/mESC_all_merged_loops_Mustache_final.bedpe', sep='\t', header=None)

loops_df = pd.DataFrame(columns=['chr','left','right','size'])
loops_df['chr'] = loops_bedpe[0]
loops_df['left'] = (loops_bedpe[1]+loops_bedpe[2])//2
loops_df['right'] = (loops_bedpe[4]+loops_bedpe[5])//2
loops_df['size'] = loops_df['right']-loops_df['left']

# parameters
min_read_counts_per_pixel = 0.4
local_region_size = 10000

# calculate matrices that depend on local_region_size
a = local_region_size//res
y_px, x_px = np.meshgrid(np.arange(-a, a+1),np.arange(-a, a+1))
    
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
    
def acceptable_size_and_location(clr, chrom, left, right, min_size=32000, convert_coords_function=None):
    '''
    Check that loop is an appropriate size (>=min_size) and is not too close to end of chromosome.
    Requires global variable local_region_size.
    Returns True if passed, False if failed.
    '''
    # convert coordinates if necessary
    if convert_coords_function is not None:
        left = convert_coords_function(chrom, left)
        right = convert_coords_function(chrom, right)
        
    if right-left < min_size:
        return False
    if left-local_region_size < 0 or right+local_region_size > clr.chromsizes[chrom]:
        return False
    return True

def no_NaNs_near_center(clr, chrom, left, right, na_stripe_dist_to_center_px_cutoff=5, convert_coords_function=None):
    '''
    Check that there are no NaN stripes too close to the center (<=na_stripe_dist_to_center_px_cutoff away).
    Requires global variable local_region_size.
    Returns True if passed, False if failed.
    '''
    # convert coordinates if necessary
    if convert_coords_function is not None:
        left = convert_coords_function(chrom, left)
        right = convert_coords_function(chrom, right)
        
    # get the image
    img = clr.matrix().fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')
    
    # find NA stripes, if any
    length_of_img = img.shape[0]  # height/width of square image
    ver_na_stripe_indices = np.where(np.sum(np.isnan(img),0)==length_of_img)[0]  # get indices of NA stripes
    hor_na_stripe_indices = np.where(np.sum(np.isnan(img),1)==length_of_img)[0]
    any_na_stripes = len(ver_na_stripe_indices)>0 or len(hor_na_stripe_indices)>0

    if any_na_stripes:
        middle_index = length_of_img//2
        ver_na_stripe_indices_from_middle = np.abs(ver_na_stripe_indices - middle_index)
        hor_na_stripe_indices_from_middle = np.abs(hor_na_stripe_indices - middle_index)
        if np.any(ver_na_stripe_indices_from_middle<=na_stripe_dist_to_center_px_cutoff) or np.any(hor_na_stripe_indices_from_middle<=na_stripe_dist_to_center_px_cutoff):
            return False  # NA values are too close to center; can't be resolved
        
    return True

def read_counts_per_pixel(clr, chrom, left, right, convert_coords_function=None):
    '''
    Calculate the number of reads divided by the number of pixels in the local region.
    Requires global variable local_region_size.
    '''
    # convert coordinates if necessary
    if convert_coords_function is not None:
        left = convert_coords_function(chrom, left)
        right = convert_coords_function(chrom, right)
        
    img_unbalanced = clr.matrix(balance=False).fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')
    read_count = np.sum(img_unbalanced)
    num_pixels = np.size(img_unbalanced)
    return read_count/num_pixels

def global_maximum_dist_to_center(clr, chrom, left, right, P_s_data, s_px_matrix, convert_coords_function=None, gaussian_blur_sigma_px=2.5, ignore_diag_cutoff_px=5):
    '''
    Calculate the Euclidean distance (in pixels) of the global maximum to the center of the image (the location of the loop). The global maximum is calculated on the observed/expected matrix.
    Requires global variable local_region_size.
    '''
    # convert coordinates if necessary
    if convert_coords_function is not None:
        left = convert_coords_function(chrom, left)
        right = convert_coords_function(chrom, right)
        
    # get the image
    img = clr.matrix().fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')

    # get the expected global background image
    bg_img = P_s_data[s_px_matrix]

    # in the image and background image, make all pixels near diagonal NA
    img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan
    bg_img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan

    # divide the image by the expected global background
    img_over_bg = img/bg_img

    # resolve any NA values (only do this if no NaNs near center)
    img_over_bg_NAs_removed = np.nan_to_num(img_over_bg, nan=np.nanmedian(img))  # replace NA values with median value in the image
    
    # blur image
    ksize = int(np.ceil(3*gaussian_blur_sigma_px)//2*2+1)  # round up to next odd integer >= 3 sigma
    img_over_bg_blurred = cv2.GaussianBlur(img_over_bg_NAs_removed,ksize=(ksize,ksize),sigmaX=gaussian_blur_sigma_px)
    
    # find global maximum
    center_pixel_indices = np.array([i[0] for i in np.where(np.logical_and(x_px==0,y_px==0))])
    brightest_pixel_indices = np.array(np.unravel_index(np.nanargmax(img_over_bg_blurred), img_over_bg_blurred.shape))
    dist_to_brightest_pixel = np.linalg.norm(brightest_pixel_indices-center_pixel_indices)

    return dist_to_brightest_pixel

def get_image(clr, chrom, left, right, P_s_data=None, over_background=False, convert_coords_function=None, ignore_diag_cutoff_px=5):
    '''
    Get the image (for diagnostic purposes).
    Requires global variable local_region_size (and s_px_matrix if over_background==True).
    '''
    # convert coordinates if necessary
    if convert_coords_function is not None:
        left = convert_coords_function(chrom, left)
        right = convert_coords_function(chrom, right)
        
    # get the image
    img = clr.matrix().fetch(f'{chrom}:{left-local_region_size}-{left+local_region_size}',f'{chrom}:{right-local_region_size}-{right+local_region_size}').astype('float')
    
    # make all pixels near diagonal NA
    img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan
        
    if over_background:
        # get the expected global background image
        bg_img = P_s_data[s_px_matrix]
        
        # make all pixels near diagonal NA
        bg_img[s_px_matrix<=ignore_diag_cutoff_px] = np.nan

        # divide the image by the expected global background
        img_over_bg = img/bg_img
        
        return img_over_bg
    else:
        return img
    
def run_loop_filtering(k):

    # get loop
    loop = loops_df.loc[k]
    chrom = loop['chr']
    left = loop['left']
    right = loop['right']
    size = loop['size']
    # adjust left and right to be in bin centers
    left = left//res*res + res//2
    right = right//res*res + res//2
    
    # calculate s_px_matrix
    loop_size_px = right//res-left//res
    s_px_matrix = loop_size_px+y_px-x_px  # genomic separation in units of res
    s_px_matrix[s_px_matrix<0] = 0  # don't allow negative values of s

    size_loc_pass = acceptable_size_and_location(coolers['all_merged'], chrom, left, right)
    NaN_pass_T1 = no_NaNs_near_center(coolers['T1'], chrom, left, right, convert_coords_function=to_T_coords)
    NaN_pass_T2 = no_NaNs_near_center(coolers['T2'], chrom, left, right, convert_coords_function=to_T_coords)
    NaN_pass_C3 = no_NaNs_near_center(coolers['C3'], chrom, left, right, convert_coords_function=to_C_coords)
    NaN_pass_C4 = no_NaNs_near_center(coolers['C4'], chrom, left, right, convert_coords_function=to_C_coords)
    NaN_pass_all_merged = no_NaNs_near_center(coolers['all_merged'], chrom, left, right)

    passed_step_1 = np.all([size_loc_pass, NaN_pass_T1, NaN_pass_T2, NaN_pass_C3, NaN_pass_C4, NaN_pass_all_merged])

    if not passed_step_1:
        return size_loc_pass, NaN_pass_T1, NaN_pass_T2, NaN_pass_C3, NaN_pass_C4, NaN_pass_all_merged, None, None, None, None, None

    read_counts_per_pixel_T1 = read_counts_per_pixel(coolers['T1'], chrom, left, right, convert_coords_function=to_T_coords)
    read_counts_per_pixel_T2 = read_counts_per_pixel(coolers['T2'], chrom, left, right, convert_coords_function=to_T_coords)
    read_counts_per_pixel_C3 = read_counts_per_pixel(coolers['C3'], chrom, left, right, convert_coords_function=to_C_coords)
    read_counts_per_pixel_C4 = read_counts_per_pixel(coolers['C4'], chrom, left, right, convert_coords_function=to_C_coords)

    passed_step_2 = np.all([read_counts_per_pixel_T1>min_read_counts_per_pixel,
                           read_counts_per_pixel_T2>min_read_counts_per_pixel,
                           read_counts_per_pixel_C3>min_read_counts_per_pixel,
                           read_counts_per_pixel_C4>min_read_counts_per_pixel])

    if not passed_step_2:
        return size_loc_pass, NaN_pass_T1, NaN_pass_T2, NaN_pass_C3, NaN_pass_C4, NaN_pass_all_merged, read_counts_per_pixel_T1, read_counts_per_pixel_T2, read_counts_per_pixel_C3, read_counts_per_pixel_C4, None

    global_max_dist_all_merged = global_maximum_dist_to_center(coolers['all_merged'], chrom, left, right, P_s_curves['all_merged'], s_px_matrix)
    
    return size_loc_pass, NaN_pass_T1, NaN_pass_T2, NaN_pass_C3, NaN_pass_C4, NaN_pass_all_merged, read_counts_per_pixel_T1, read_counts_per_pixel_T2, read_counts_per_pixel_C3, read_counts_per_pixel_C4, global_max_dist_all_merged

# initialize dataframe for storing results
loop_filtering_criteria_df = pd.DataFrame(index=loops_df.index, columns=['size_loc_pass','NaN_pass_T1','NaN_pass_T2','NaN_pass_C3','NaN_pass_C4','NaN_pass_all_merged','read_counts_per_pixel_T1','read_counts_per_pixel_T2','read_counts_per_pixel_C3','read_counts_per_pixel_C4','global_max_dist_all_merged'])

# run parameters
chunk_size = 40
nproc = 40

# set up multiprocessing
num_chunks = int(np.ceil(len(loops_df)/chunk_size))
chunk_starts = np.arange(num_chunks)*chunk_size
chunk_ends = (np.arange(num_chunks)+1)*chunk_size
chunk_ends[-1] = len(loops_df)

for chunk_index in np.arange(num_chunks):
    start = time.time()
    with Pool(nproc) as p:
        indices_in_chunk = np.arange(chunk_starts[chunk_index],chunk_ends[chunk_index])
        results_in_chunk = p.map(run_loop_filtering, indices_in_chunk)
        loop_filtering_criteria_df.loc[indices_in_chunk] = results_in_chunk
    end = time.time()

    # print progress
    print(f'chunk {chunk_index+1} of {num_chunks} ({end-start:.2f} s)')

    # save periodically
    if chunk_index % 100 == 0:
        loop_filtering_criteria_df.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/loop_filtering_criteria.csv')
loop_filtering_criteria_df.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/loop_filtering_criteria.csv')  # save everything when finished

filtered_loops_df = loops_df.loc[loop_filtering_criteria_df['global_max_dist_all_merged']<=2.5]
filtered_loops_df.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/filtered_loops.csv')