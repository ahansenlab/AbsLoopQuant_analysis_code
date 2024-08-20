# load coolers
coolers = {}
res = 1000
coolers['T1'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/T1_final/T1.mcool::/resolutions/{res}')
coolers['T2'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/T2_final/T2.mcool::/resolutions/{res}')
coolers['C3'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/C3_final/C3.mcool::/resolutions/{res}')
coolers['C4'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/20230912_MicroC/C4_final/C4.mcool::/resolutions/{res}')
coolers['all_merged'] = cooler.Cooler(f'/mnt/coldstorage/jjusuf/Past_MicroC_Experiments/all_WTgenome/mESC_all_merged.mcool::/resolutions/{res}')

# calculate P(s) curves
for sample_name, clr in coolers.items():
    looptools.calculate_and_save_avg_Ps_curve(clr, nproc=50, output_filename=f'/mnt/md0/jjusuf/absloopquant/AbsLoopQuant_data/P_s_{sample_name}_1000bp.txt')