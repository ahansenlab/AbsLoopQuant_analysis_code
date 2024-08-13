# Randomly remove trans reads in a .pairs file by a user-defined probability (the scaling factor) in order to make the cis/trans ratio constant across samples. Reads where one or both sides are unmapped are not affected.

# arguments:
# (1) input filename
# (2) scaling factor (must be between 0 and 1) by which to downsample the trans reads
# (3) output filename

import sys
import subprocess as sp

input_filename = sys.argv[1]
scaling_factor = sys.argv[2]
output_filename = sys.argv[3]

if input_filename.endswith('.gz'):
    sp.run(f'''gzip -cd {input_filename} | awk 'BEGIN {{ srand() }} {{ if ($1 ~ /^#/ || $2 == "!" || $4 == "!" || $2 == $4 || rand() < {scaling_factor}) print $0 }}' FS='\t' | bgzip > {output_filename}''', shell=True)
else:
    sp.run(f'''awk 'BEGIN {{ srand() }} {{ if ($1 ~ /^#/ || $2 == "!" || $4 == "!" || $2 == $4 || rand() < {scaling_factor}) print $0 }}' FS='\t' {input_filename} | bgzip > {output_filename}''', shell=True)