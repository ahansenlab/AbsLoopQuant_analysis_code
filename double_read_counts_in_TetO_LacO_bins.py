# For a contact matrix containing raw read counts (.cool file) for the TetO-LacO+3xCTCF cell line from Mach et al. 2022, this Python script doubles the number of read counts in the bins corresponding to the TetO+3xCTCF and 3xCTCF+LacO sequences (which are heterozygous) to make the read counts as would be expected from homozygous insertions. This step is performed before zoomifying/balancing the cooler.

# arguments:
# (1) filename for input cooler (assumes it is a .cool file with 250 bp bin sizes)
# (2) desired filename for output cooler

import numpy as np
import cooler
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

clr = cooler.Cooler(input_filename)

# coordinates of synthetic regions rounded to nearest 250 bp
LacO_start_coord_250 = 11567250
LacO_end_coord_250 = 11572250
TetO_start_coord_250 = 11722750
TetO_end_coord_250 = 11732250

LacO_start_bin_250, LacO_end_bin_250 = clr.extent(f'chr15:{LacO_start_coord_250}-{LacO_end_coord_250}')
TetO_start_bin_250, TetO_end_bin_250 = clr.extent(f'chr15:{TetO_start_coord_250}-{TetO_end_coord_250}')

chunksize = 10000000

spans = cooler.util.partition(0, len(clr.pixels()), chunksize)

def get_modified_pixels(lo, hi):
    pixels = clr.pixels()[lo:hi]
    bin1_in_LacO = np.logical_and(pixels['bin1_id']>=LacO_start_bin_250, pixels['bin1_id']<LacO_end_bin_250)
    bin2_in_LacO = np.logical_and(pixels['bin2_id']>=LacO_start_bin_250, pixels['bin2_id']<LacO_end_bin_250)
    bin1_in_TetO = np.logical_and(pixels['bin1_id']>=TetO_start_bin_250, pixels['bin1_id']<TetO_end_bin_250)
    bin2_in_TetO = np.logical_and(pixels['bin2_id']>=TetO_start_bin_250, pixels['bin2_id']<TetO_end_bin_250)
    pixels.loc[np.any((bin1_in_LacO,bin2_in_LacO,bin1_in_TetO,bin2_in_TetO),0),'count'] *= 2
    return pixels

chunk_generator = (get_modified_pixels(lo, hi) for lo, hi in spans)
cooler.create_cooler(output_filename, clr.bins()[:], chunk_generator)