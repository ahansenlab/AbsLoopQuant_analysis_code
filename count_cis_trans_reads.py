# Count the number of cis reads and trans reads in a .pairs file. If either side is unmapped (chr=! in the pairs file), the pair will be skipped.
# Returns # of cis reads, # of trans reads, # of total reads, % of cis reads, and % of trans reads

import sys
import subprocess as sp

filename = sys.argv[1]

if filename.endswith('.gz'):
    result_cis = sp.run(f"""gzip -cd {filename} | awk 'BEGIN {{ FS="\t"; count=0 }} $1 !~ /^#/ && $2 != "!" && $4 != "!" && $2 == $4 {{ count++ }} END {{ print count }}'""", shell=True, capture_output=True)
    result_trans = sp.run(f"""gzip -cd {filename} | awk 'BEGIN {{ FS="\t"; count=0 }} $1 !~ /^#/ && $2 != "!" && $4 != "!" && $2 != $4 {{ count++ }} END {{ print count }}'""", shell=True, capture_output=True)
else:
    result_cis = sp.run(f"""awk 'BEGIN {{ FS="\t"; count=0 }} $1 !~ /^#/ && $2 != "!" && $4 != "!" && $2 == $4 {{ count++ }} END {{ print count }}' {filename}""", shell=True, capture_output=True)
    result_trans = sp.run(f"""awk 'BEGIN {{ FS="\t"; count=0 }} $1 !~ /^#/ && $2 != "!" && $4 != "!" && $2 != $4 {{ count++ }} END {{ print count }}' {filename}""", shell=True, capture_output=True)

num_cis_reads = int(result_cis.stdout)
num_trans_reads = int(result_trans.stdout)
num_total_reads = num_cis_reads + num_trans_reads

print(f'#cis\t{num_cis_reads}')
print(f'#trans\t{num_trans_reads}')
print(f'#total\t{num_total_reads}')
print(f'%cis\t{num_cis_reads/num_total_reads*100:.2f}')
print(f'%trans\t{num_trans_reads/num_total_reads*100:.2f}')