import numpy as np
import pandas as pd
import cooler
import pybedtools

# import the loop quantification module
import sys
sys.path.insert(1, '/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_analysis_code')
import looptools  # all the back end code is here

loops = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/filtered_loops.csv', index_col=0)

promoters_df = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/mm39_TSS_pad_2kb.bed',sep='\t', header=None, names=('chr','start','end'))

enhancers_df = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/mm39_enhancers.bed', sep='\t', header=None, names=('chr','start','end'), index_col=None, usecols=[0,1,2])

CTCF_SMC1A_df = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/mm39_CTCF_cohesin.bed',sep='\t', header=None, names=('chr','start','end','strand'), usecols=[0,1,2,5])
CTCF_SMC1A_pos_df = CTCF_SMC1A_df.loc[CTCF_SMC1A_df['strand']=='+']
CTCF_SMC1A_neg_df = CTCF_SMC1A_df.loc[CTCF_SMC1A_df['strand']=='-']

enhancers_bed = pybedtools.BedTool(enhancers_df.to_string(header=False, index=False), from_string=True)
promoters_bed = pybedtools.BedTool(promoters_df.to_string(header=False, index=False), from_string=True)
CTCF_SMC1A_bed = pybedtools.BedTool(CTCF_SMC1A_df.to_string(header=False, index=False), from_string=True)
CTCF_SMC1A_pos_bed = pybedtools.BedTool(CTCF_SMC1A_pos_df.to_string(header=False, index=False), from_string=True)
CTCF_SMC1A_neg_bed = pybedtools.BedTool(CTCF_SMC1A_neg_df.to_string(header=False, index=False), from_string=True)

loops_bed_df = loops[np.array(['chr','left','right'])]
loops_bed = pybedtools.BedTool(loops_bed_df.to_string(header=False, index=False), from_string=True)

left_anchors_bed = pybedtools.BedTool(pd.concat([loops['chr'], loops['left']-2500, loops['left']+2500, pd.Series(loops.index, index=pd.Series(loops.index))], axis=1).to_string(header=False, index=False), from_string=True)
right_anchors_bed = pybedtools.BedTool(pd.concat([loops['chr'], loops['right']-2500, loops['right']+2500, pd.Series(loops.index, index=pd.Series(loops.index))], axis=1).to_string(header=False, index=False), from_string=True)

# promoter at left anchors
intersection = left_anchors_bed.intersect(promoters_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'P_L']=True
loops.loc[pd.isna(loops['P_L']),'P_L']=False

# promoter at right anchors
intersection = right_anchors_bed.intersect(promoters_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'P_R']=True
loops.loc[pd.isna(loops['P_R']),'P_R']=False

# enhancer at left anchors
intersection = left_anchors_bed.intersect(enhancers_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'E_L']=True
loops.loc[pd.isna(loops['E_L']),'E_L']=False

# enhancer at right anchors
intersection = right_anchors_bed.intersect(enhancers_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'E_R']=True
loops.loc[pd.isna(loops['E_R']),'E_R']=False

# positive CTCF at left anchors (convergent)
intersection = left_anchors_bed.intersect(CTCF_SMC1A_pos_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'C_pos_L']=True
loops.loc[pd.isna(loops['C_pos_L']),'C_pos_L']=False

# negative CTCF at left anchors (divergent)
intersection = left_anchors_bed.intersect(CTCF_SMC1A_neg_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'C_neg_L']=True
loops.loc[pd.isna(loops['C_neg_L']),'C_neg_L']=False

# any CTCF at left anchors
loops.loc[:,'C_L']=np.logical_or(loops['C_pos_L'], loops['C_neg_L'])

# positive CTCF at right anchors (divergent)
intersection = right_anchors_bed.intersect(CTCF_SMC1A_pos_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'C_pos_R']=True
loops.loc[pd.isna(loops['C_pos_R']),'C_pos_R']=False

# negative CTCF at right anchors (convergent)
intersection = right_anchors_bed.intersect(CTCF_SMC1A_neg_bed, u=True).to_dataframe()
loops.loc[intersection['name'].values,'C_neg_R']=True
loops.loc[pd.isna(loops['C_neg_R']),'C_neg_R']=False

# any CTCF at right anchors
loops.loc[:,'C_R']=np.logical_or(loops['C_pos_R'], loops['C_neg_R'])

# set datatype for logical columns to bool
for i in loops.columns[4:]:
    loops[i] = loops[i].astype('bool')

loops['pure_PP'] = np.all((loops.P_L,  # promoter L, nothing else
                           loops.E_L==False,
                           loops.C_L==False,
                           loops.P_R,  # promoter R, nothing else
                           loops.E_R==False,
                           loops.C_R==False), 0)
print(f'pure PP: {np.sum(loops["pure_PP"])}')

loops['pure_EE'] = np.all((loops.E_L,  # enhancer L, nothing else
                           loops.P_L==False,
                           loops.C_L==False,
                           loops.E_R,  # enhancer R, nothing else
                           loops.P_R==False,
                           loops.C_R==False), 0)
print(f'pure EE: {np.sum(loops["pure_EE"])}')

loops['pure_EP'] = np.logical_or(np.all((loops.P_L,  # promoter L, nothing else
                                         loops.E_L==False,
                                         loops.C_L==False,
                                         loops.E_R,  # enhancer R, nothing else
                                         loops.P_R==False,
                                         loops.C_R==False), 0),
                                 np.all((loops.E_L,  # enhancer L, nothing else
                                         loops.P_L==False,
                                         loops.C_L==False,
                                         loops.P_R,  # promoter R, nothing else
                                         loops.E_R==False,
                                         loops.C_R==False), 0))
print(f'pure EP: {np.sum(loops["pure_EP"])}')

loops['mixed_CRE'] = np.all((np.logical_or(loops.P_L, loops.E_L),  # promoter or enhancer L, not CTCF
                            loops.C_L==False,
                            np.logical_or(loops.P_R, loops.E_R),  # promoter or enhancer R, not CTCF
                            loops.C_R==False,
                            loops.pure_PP==False,  # must not be pure PP/pure EP/pure EE
                            loops.pure_EP==False,
                            loops.pure_EE==False), 0)
print(f'mixed CRE: {np.sum(loops["mixed_CRE"])}')

loops['pure_CRE'] = np.any((loops.pure_PP, loops.pure_EP, loops.pure_EE, loops.mixed_CRE), 0)
print(f'pure CRE: {np.sum(loops["pure_CRE"])}')
print()

loops['pure_CC_convergent'] = np.all((loops.C_pos_L,  # presence of CTCF pos at L, no enhancer/promoter
                                      loops.P_L==False,
                                      loops.E_L==False,
                                      loops.C_neg_R,  # presence of CTCF neg at R, no enhancer/promoter
                                      loops.P_R==False,
                                      loops.E_R==False), 0)
print(f'pure CC convergent: {np.sum(loops["pure_CC_convergent"])}')

loops['pure_CC_nonconvergent'] = np.all((loops.C_L,  # CTCF L, no enhancer/promoter
                                         loops.P_L==False,
                                         loops.E_L==False,
                                         loops.C_R,  # CTCF R, no enhancer/promoter
                                         loops.P_R==False,
                                         loops.E_R==False,
                                         loops.pure_CC_convergent==False  # must not have convergent pair
                                        ), 0)
print(f'pure CC nonconvergent: {np.sum(loops["pure_CC_nonconvergent"])}')

loops['pure_CC'] = np.logical_or(loops['pure_CC_convergent'], loops['pure_CC_nonconvergent'])
print(f'pure CC: {np.sum(loops["pure_CC"])}')
print()

two_sided = np.all((np.any((loops.C_L,loops.P_L,loops.E_L), 0),  # feature L
                    np.any((loops.C_R,loops.P_R,loops.E_R), 0)), 0)  # feature R
has_CRE = np.any((loops.P_L, loops.P_R, loops.E_L, loops.E_R), 0)
has_C = np.logical_or(loops.C_L, loops.C_R)
loops['mixed'] = np.all((two_sided, has_C, has_CRE), 0)
print(f'mixed: {np.sum(loops["mixed"])}')
print()

loops['one_sided'] = np.logical_xor(np.any((loops.C_L,loops.P_L,loops.E_L), 0),  # feature L
                                    # XOR
                                    np.any((loops.C_R,loops.P_R,loops.E_R), 0),  # feature R
                                   )
print(f'one-sided: {np.sum(loops["one_sided"])}')

loops['other_other'] = np.all((loops.C_L==False,  # no features L
                               loops.P_L==False,
                               loops.E_L==False,
                               loops.C_R==False,  # no features R
                               loops.P_R==False,
                               loops.E_R==False), 0)
print(f'other-other: {np.sum(loops["other_other"])}')

# append columns for loop scores
scores = pd.read_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/MicroC_loop_quantification_scores.csv', index_col=0)
scores['abs_loop_percent'] = scores.mean(1)*18.6  # slope of calibration curve from Figure 1

# save final dataframe containing all info
final_df = pd.concat((loops, scores), axis=1)
final_df.to_csv('/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/final_loop_data.csv')
