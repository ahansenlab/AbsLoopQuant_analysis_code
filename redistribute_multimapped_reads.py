# script to redistribute multiply mapped reads in Micro-C, version 3 (updated June 2024)
# use --help option for more info

import numpy as np
import pandas as pd
from io import StringIO
from pyfaidx import Fasta
from Bio.Seq import Seq
import warnings
import re
import sys
import argparse

parser = argparse.ArgumentParser(description = "distribute multiply mapped reads among possible positions randomly")
parser.add_argument("--name", "-n", help="sample name, e.g., T1")
parser.add_argument("--filenames", "-f", help="comma-separated list of input .pairs files with pairs to be redistributed, e.g., 'T1_in_selected_regions_1.pairs,T1_in_selected_regions_2.pairs,...,T1_in_selected_regions_8.pairs'")
parser.add_argument("--rois", "-r", help="regions of interest that could contain the primary hits of the multimapping reads, formatted as a comma-separated list of chr,start,end, e.g., 'chr15,11722732,11732141,chr15,11567241,11572268,chr1,34257542,34257653,chr4,132978078,132978190,chr8,13511989,13512092'")
parser.add_argument("--genome", "-g", help="path to genome file, e.g., /mnt/md0/jjusuf/genomes/mm39/mm39_TetO_LacO_3xCTCF.fasta")
parser.add_argument("--min_mapq", "-m", help="minimum mapq score that was used to classify reads as unique; use the same value as in main alignment (usually 30)", type=int)
parser.add_argument("--verbose", "-v", help="run in verbose mode", action='store_true')
args = parser.parse_args()

sample_name = args.name
filenames_string = args.filenames
ROI_string = args.rois
genome_path = args.genome
min_mapq = args.min_mapq
verbose = args.verbose

# other setup
genome = Fasta(genome_path)
filenames = filenames_string.split(',')
regions_of_interest = [(ROI_string.split(',')[i], int(ROI_string.split(',')[i+1]), int(ROI_string.split(',')[i+2])) for i in range(0,len(ROI_string.split(',')),3)]
list_of_region_seqs = [genome[c[0]][c[1]:c[2]].seq.upper() for c in regions_of_interest]

# define functions
def open_pairs_file(filename, dropsam=False):
    """
    Open a .pairs file and return a pandas dataframe containing the data.
    """
    
    f = open(filename)
    s = f.read()
    lines = s.split('\n')
    
    # remove comment lines
    for i in range(len(lines)):
        if not lines[i].startswith('#'):
            break
    last_comment_line = lines[i-1]
    column_names = last_comment_line.split()[1:]
    lines = lines[i:]

    # remove last line if empty
    if lines[-1] == '':
        lines = lines[:-1]
        
    df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t', header=None, names=column_names)
    
    return df

def add_pos0L_columns(df):
    """
    Add columns for the 0-based left position of the reads in a pandas dataframe of pairs data. Returns the modified dataframe.
    """
    
    for i in df.index:
        for j in (1, 2):
            if df.loc[i, f'strand{j}']=='+':
                df.loc[i, f'pos0L{j}'] = df.loc[i, f'pos{j}'] - 1
            if df.loc[i, f'strand{j}']=='-':
                df.loc[i, f'pos0L{j}'] = df.loc[i, f'pos{j}'] - df.loc[i, f'algn_ref_span{j}']
    
    df['pos0L1'] = df['pos0L1'].astype('int')
    df['pos0L2'] = df['pos0L2'].astype('int')
    
    return df

def add_pos_columns(df):
    """
    Add columns for the 1-based 5' position of the reads in a pandas dataframe, based on existing columns for 0-based left positions.
    """
    
    for i in df.index:
        for j in (1, 2):
            if df.loc[i, f'strand{j}']=='+':
                df.loc[i, f'pos{j}'] = df.loc[i, f'pos0L{j}'] + 1
            if df.loc[i, f'strand{j}']=='-':
                df.loc[i, f'pos{j}'] = df.loc[i, f'pos0L{j}'] + df.loc[i, f'algn_ref_span{j}']
                
    df['pos1'] = df['pos1'].astype('int')
    df['pos2'] = df['pos2'].astype('int')
    
    return df

def convert_pos_0L(h, algn_ref_span):
    """
    Convert the 0-based left position of a hit into the 1-based 5' position. The hit should be given as a tuple of the form (chrom, pos0L, strand).
    """
    chrom, pos0L, strand = h
    if strand=='+':
        return pos0L + 1
    if strand=='-':
        return pos0L + algn_ref_span

def get_seq(chrom, pos0L, refspan, genome):
    """
    Return the sequence of a read given the chromosome, 0-based left position, span in the reference genome, and the reference genome itself as a Fasta object
    """
    
    if chrom == '!':  # null mapped
        return None
    seq = genome[chrom][pos0L:pos0L+refspan].seq.upper()
    return seq

def rc(seq):
    """
    Return the reverse complement of a sequence.
    """
    
    return str(Seq(seq).reverse_complement())

def find_all_hits(seq, strand, genome, regions_of_interest):
    """
    Find all "hits" of a sequence within the regions of interest.
    
    Arguments:
        - sequence (str)
        - the strand of the sequence ('+' or '-')
        - the reference genome as a Fasta object
        - a list of all regions_of_interest to search within, given as a list of tuples of the form (chrom, pos_start, pos_end).
          pos_start and pos_end should be 0-based
          
    Returns:
        - hits given as a list of tuples of the form (chrom, pos0L, strand), where pos0L is the 0-based left position of the hit
    """
    
    list_of_region_seqs = [genome[c[0]][c[1]:c[2]].seq.upper() for c in regions_of_interest]
    list_of_region_chroms = [c[0] for c in regions_of_interest]
    list_of_region_start_pos = [c[1] for c in regions_of_interest]
    hits = []
    
    if strand == '+':
        opp_strand = '-'
    else:  # strand == '-'
        opp_strand = '+'
        
    if seq != None:
        for i in range(len(list_of_region_seqs)):
            region_seq = list_of_region_seqs[i]
            region_chrom = list_of_region_chroms[i]
            region_start_pos = list_of_region_start_pos[i]
            hits += [(region_chrom, match.start(1)+region_start_pos, strand) for match in re.finditer(rf'(?=({seq}))', region_seq)]
            hits += [(region_chrom, match.start(1)+region_start_pos, opp_strand) for match in re.finditer(rf'(?=({rc(seq)}))', region_seq)]
    return hits

def P_s(dist):
    """
    Returns the relative probability of contact of two loci on the same chromosome given their genomic separation (dist).
    The absolute scaling of this function does not matter because the probabilities will be normalized later.
    """
    
    if dist < 1000:
        return 1/1000
        
    return 1/dist

# initialize dataframe of pairs
pairs_data = pd.DataFrame(columns=('readID','chrom1','pos1','chrom2','pos2','strand1','strand2','pair_type','mapq1','mapq2','algn_ref_span1','algn_ref_span2','sam1','sam2','pos0L1','pos0L2'))
for i in range(len(filenames)):
    new_pairs_data = open_pairs_file(filenames[i], dropsam=False)
    new_pairs_data = add_pos0L_columns(new_pairs_data)
    pairs_data = pd.concat((pairs_data, new_pairs_data), ignore_index=True)
    
# initialize dataframe to store newly distributed pairs
pairs_data_redist = pd.DataFrame(index=pairs_data.index, columns=('readID','chrom1','pos1','chrom2','pos2','strand1','strand2','pair_type','mapq1','mapq2','algn_ref_span1','algn_ref_span2','sam1','sam2','pos0L1','pos0L2'))

# redistribute pairs
for i in pairs_data.index:
    pair = pairs_data.loc[i]
    if verbose:
        print(f"pair {i}")
        print(f"side1 raw: chrom={pair['chrom1']}, pos={pair['pos1']}, strand={pair['strand1']}, mapq={pair['mapq1']}, algn_ref_span={pair['algn_ref_span1']}")
        print(f"side2 raw: chrom={pair['chrom2']}, pos={pair['pos2']}, strand={pair['strand2']}, mapq={pair['mapq2']}, algn_ref_span={pair['algn_ref_span2']}")
        print(f"side1 primary hit: ({pair['chrom1']},{pair['pos0L1']},{pair['strand1']})")
        print(f"side2 primary hit: ({pair['chrom2']},{pair['pos0L2']},{pair['strand2']})")
        
    # skip pair if either side is null
    if 'N' in pair['pair_type'].upper():
        if verbose:
            print('one or both sides is null, skipping...\n')
        continue
        
    # find hits
    seq1 = get_seq(pair['chrom1'], pair['pos0L1'], pair['algn_ref_span1'], genome)
    seq2 = get_seq(pair['chrom2'], pair['pos0L2'], pair['algn_ref_span2'], genome)
    hits1 = find_all_hits(seq1, pair['strand1'], genome, regions_of_interest)
    hits2 = find_all_hits(seq2, pair['strand2'], genome, regions_of_interest)
    primary_hit1 = (pair['chrom1'], pair['pos0L1'], pair['strand1'])
    if len(hits1)==0:
        hits1.append(primary_hit1)
    else:
        if not primary_hit1 in hits1:
            warnings.warn(f'The primary alignment of side1 ({primary_hit1}) was not detected as a hit by the script!', stacklevel=2)
    primary_hit2 = (pair['chrom2'], pair['pos0L2'], pair['strand2'])
    if len(hits2)==0:
        hits2.append(primary_hit2)
    else:
        if not primary_hit2 in hits2:
            warnings.warn(f'The primary alignment of side2 ({primary_hit2}) was not detected as a hit by the script!')
    if verbose:
        print(f"side1 all hits: {' '.join([f'({h[0]},{h[1]},{h[2]})' for h in hits1])}")
        print(f"side2 all hits: {' '.join([f'({h[0]},{h[1]},{h[2]})' for h in hits2])}")

    # skip pair if either side is multiply mapped outside ROI (unresolvable)
    if (pair['mapq1'] < min_mapq and len(hits1)==1) or (pair['mapq2'] < min_mapq and len(hits2)==1):
        if verbose:
            print('one or both sides is multiply mapped outside the ROIs, skipping...\n')
        continue

    # redistribute read
    chroms_in_hits1 = np.unique([h[0] for h in hits1])
    chroms_in_hits2 = np.unique([h[0] for h in hits2])
    chroms_in_both = np.intersect1d(chroms_in_hits1, chroms_in_hits2)
    if len(chroms_in_both)==0:  # read must be trans
        if verbose:
            print('proceeding as trans read; assigning both sides randomly with uniform distribution')
        chosen_chrom1 = np.random.choice(chroms_in_hits1)
        chosen_chrom2 = np.random.choice(chroms_in_hits2)
        hits1 = [h for h in hits1 if h[0]==chosen_chrom1]
        hits2 = [h for h in hits2 if h[0]==chosen_chrom2]
        chosen_hit1 = hits1[np.random.choice(range(len(hits1)))]
        chosen_hit2 = hits2[np.random.choice(range(len(hits2)))]
    else:  # read is cis (or assume it is cis with high probability)
        if verbose:
            print('proceeding as cis read\nall possible pairs and probabilities:')
        hits1 = [h for h in hits1 if h[0] in chroms_in_both]
        hits2 = [h for h in hits2 if h[0] in chroms_in_both]

        num_possible_pairs = len(hits1) * len(hits2)
        possible_pairs_side1 = [None] * num_possible_pairs
        possible_pairs_side2 = [None] * num_possible_pairs
        probabilities = [0] * num_possible_pairs
        c = 0
        for j in range(len(hits1)):
            for k in range(len(hits2)):
                possible_pairs_side1[c] = hits1[j]
                possible_pairs_side2[c] = hits2[k]
                c += 1
        for c in range(num_possible_pairs):
            side1 = possible_pairs_side1[c]
            side2 = possible_pairs_side2[c]
            if side1[0]==side2[0]:  # cis contact
                side1_pos0M = side1[1]+(pair['algn_ref_span1']-1)/2 if side1[2]=='+' else side1[1]-(pair['algn_ref_span1']-1)/2
                side2_pos0M = side2[1]+(pair['algn_ref_span2']-1)/2 if side2[2]=='+' else side2[1]-(pair['algn_ref_span2']-1)/2
                dist = np.abs(side1_pos0M-side2_pos0M)
                probabilities[c] = P_s(dist)
            else:  # trans contact
                probabilities[c] = 0
        probabilities = probabilities / np.sum(probabilities)
        if verbose:
            for c in range(num_possible_pairs):
                side1 = possible_pairs_side1[c]
                side2 = possible_pairs_side2[c]
                print(f"{c} ({side1[0]},{side1[1]},{side1[2]}) ({side2[0]},{side2[1]},{side2[2]}) {probabilities[c]:.2e}")

        chosen_index = np.random.choice(range(num_possible_pairs), p=probabilities)
        chosen_hit1 = possible_pairs_side1[chosen_index]
        chosen_hit2 = possible_pairs_side2[chosen_index]
    
    pairs_data_redist.loc[i, 'readID'] = pairs_data.loc[i, 'readID']
    pairs_data_redist.loc[i, 'chrom1'] = chosen_hit1[0]
    pairs_data_redist.loc[i, 'pos1'] = np.nan
    pairs_data_redist.loc[i, 'chrom2'] = chosen_hit2[0]
    pairs_data_redist.loc[i, 'pos2'] = np.nan
    pairs_data_redist.loc[i, 'strand1'] = chosen_hit1[2]
    pairs_data_redist.loc[i, 'strand2'] = chosen_hit2[2]
    pairs_data_redist.loc[i, 'pair_type'] = pairs_data.loc[i, 'pair_type']
    pairs_data_redist.loc[i, 'mapq1'] = pairs_data.loc[i, 'mapq1']
    pairs_data_redist.loc[i, 'mapq2'] = pairs_data.loc[i, 'mapq2']
    pairs_data_redist.loc[i, 'algn_ref_span1'] = pairs_data.loc[i, 'algn_ref_span1']
    pairs_data_redist.loc[i, 'algn_ref_span2'] = pairs_data.loc[i, 'algn_ref_span2']
    pairs_data_redist.loc[i, 'pos0L1'] = chosen_hit1[1]
    pairs_data_redist.loc[i, 'pos0L2'] = chosen_hit2[1]
    if verbose:
        print(f"side1 chosen hit: ({chosen_hit1[0]},{chosen_hit1[1]},{chosen_hit1[2]})")
        print(f"side2 chosen hit: ({chosen_hit2[0]},{chosen_hit2[1]},{chosen_hit2[2]})")
        print(f"side1 raw new: chrom={chosen_hit1[0]}, pos={convert_pos_0L(chosen_hit1,pairs_data.loc[i, 'algn_ref_span1'])}, strand={chosen_hit1[2]}, mapq={pairs_data.loc[i, 'mapq1']}, algn_ref_span={pairs_data.loc[i, 'algn_ref_span1']}")
        print(f"side2 raw new: chrom={chosen_hit2[0]}, pos={convert_pos_0L(chosen_hit2,pairs_data.loc[i, 'algn_ref_span2'])}, strand={chosen_hit2[2]}, mapq={pairs_data.loc[i, 'mapq2']}, algn_ref_span={pairs_data.loc[i, 'algn_ref_span2']}\n")

# remove any rows that were not redistributed
pairs_data_redist = pairs_data_redist.loc[~np.all(pd.isna(pairs_data_redist), axis=1)]

# populate columns pos1 & pos2 (1-based coordinates of 5' end)
pairs_data_redist = add_pos_columns(pairs_data_redist)

# format for export
pairs_data_redist_final = pairs_data_redist.copy()
pairs_data_redist_final['mapq1']=60
pairs_data_redist_final['mapq2']=60
pairs_data_redist_final.drop(['algn_ref_span1', 'algn_ref_span2', 'pos0L1', 'pos0L2', 'sam1', 'sam2'], axis=1, inplace=True)
pairs_data_redist_final_str = pairs_data_redist_final.to_csv(sep='\t', index=False, header=False)
output_filename = f'{sample_name}_in_selected_regions_redist.pairs'
with open(output_filename, 'w') as f:
    f.write(pairs_data_redist_final_str)