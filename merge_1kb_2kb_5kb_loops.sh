mkdir -p mESC_all_merged_loops_Mustache_1kb_5kb_combined
mkdir -p mESC_all_merged_loops_Mustache_1kb_2kb_combined
mkdir -p mESC_all_merged_loops_Mustache_2kb_5kb_combined
mkdir -p mESC_all_merged_loops_Mustache_final

# get pairwise intersections using Miles' script
python /mnt/md0/Tools/scripts/loopComparer.py --bedpe1 mESC_all_merged_loops_Mustache_1kb.bedpe --bedpe2 mESC_all_merged_loops_Mustache_5kb.bedpe --names 1kb 5kb --slop 5000 --outdir mESC_all_merged_loops_Mustache_1kb_5kb_combined
python /mnt/md0/Tools/scripts/loopComparer.py --bedpe1 mESC_all_merged_loops_Mustache_1kb.bedpe --bedpe2 mESC_all_merged_loops_Mustache_2kb.bedpe --names 1kb 2kb --slop 5000 --outdir mESC_all_merged_loops_Mustache_1kb_2kb_combined
python /mnt/md0/Tools/scripts/loopComparer.py --bedpe1 mESC_all_merged_loops_Mustache_2kb.bedpe --bedpe2 mESC_all_merged_loops_Mustache_5kb.bedpe --names 2kb 5kb --slop 5000 --outdir mESC_all_merged_loops_Mustache_2kb_5kb_combined

# 1 kb only
FILE1=mESC_all_merged_loops_Mustache_1kb_2kb_combined/1kbnot2kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_1kb_5kb_combined/1kbnot5kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/1kbonly.bedpe

# 2 kb only
FILE1=mESC_all_merged_loops_Mustache_1kb_2kb_combined/2kbnot1kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_2kb_5kb_combined/2kbnot5kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/2kbonly.bedpe

# 5 kb only
FILE1=mESC_all_merged_loops_Mustache_1kb_5kb_combined/5kbnot1kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_2kb_5kb_combined/5kbnot2kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/5kbonly.bedpe

# 1 and 2 kb only (not 5 kb)
FILE1=mESC_all_merged_loops_Mustache_1kb_2kb_combined/1kband2kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_1kb_5kb_combined/1kbnot5kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/1kband2kbonly.bedpe

# 1 and 5 kb only (not 2 kb)
FILE1=mESC_all_merged_loops_Mustache_1kb_5kb_combined/1kband5kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_1kb_2kb_combined/1kbnot2kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/1kband5kbonly.bedpe

# 2 and 5 kb only (not 1 kb)
FILE1=mESC_all_merged_loops_Mustache_2kb_5kb_combined/2kband5kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_1kb_2kb_combined/2kbnot1kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/2kband5kbonly.bedpe

# all
FILE1=mESC_all_merged_loops_Mustache_1kb_2kb_combined/1kband2kb.bedpe
FILE2=mESC_all_merged_loops_Mustache_1kb_5kb_combined/1kband5kb.bedpe
pairToPair -a $FILE1 -b $FILE2 | cut -f 1-6 > mESC_all_merged_loops_Mustache_final/1kband2kband5kb.bedpe