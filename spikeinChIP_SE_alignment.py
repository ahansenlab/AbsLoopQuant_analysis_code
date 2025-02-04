#Adapted by James from Miles' original code (spikeinChIP_PE_alignment.py)
#Aligning and processing single end fastq files for spike-in ChIP-seq using bowtie2

from sys import exit
import subprocess as sp
import argparse
import multiprocessing
import uuid
import pandas as pd

parser = argparse.ArgumentParser(description = "run bowtie2 on paired end fastq files with spikein to produce aligned bam files")
parser.add_argument("--filename", "-f", help = "a tab-separated file containing one each line the path to the fastq with the reads, and the desired output name for the aligned file")
parser.add_argument("--genome", "-g", help = "genome build to align to - mouse or human - required")
parser.add_argument("--spikegenome", "-s", help = "spikein genome build to align to - mouse or human")
parser.add_argument("--threads", "-t", help = "number of threads to use for bowtie2 - default is 1", default = "1")
parser.add_argument("--outname", "-o", help = "name for the table to store counts", default = "alignmentcounts")
parser.add_argument("--outdir", help = "a directory to store output files - default is current directory", default = "./")
args = parser.parse_args()

genome = args.genome
spikegenome = args.spikegenome
threads = args.threads
outname = args.outname
outdir = args.outdir

#Read in file to determine what to process
files = pd.read_csv(args.filename, sep='\t', header=None, names = ['inputs', 'outputs'])
file = files.inputs.values.tolist()
outnames = files.outputs.values.tolist()

#Check requirements are fulfilled:
condapacks = sp.run("conda list".split(), capture_output=True)
condapacksstr = str(condapacks.stdout)

if "bowtie2" not in condapacksstr or "samtools" not in condapacksstr or "sambamba" not in condapacksstr:
    print("Please make sure bowtie2, samtools and sambamba are installed in your current conda environment (check conda list)")
    exit()

#Check that outdir ends with a /, add one if it doesn't
if args.outdir is not None and not outdir.endswith("/"):
    outdir = outdir + "/"
    
#Make outdir if it doesn't already exist
sp.run(f'mkdir {outdir}', shell=True)

#Check a genome was specified
if args.genome is None:
    print("Genome not specified - check help for formatting")
    parser.print_usage()
    exit()

#Check if file and outnames are present and the same lengths
if file is None or outnames is None:
    print("Input files not specified - check help for formatting")
    parser.print_usage()
    exit()
elif type(file) == list and type(file) == list and type(outnames) == list and len(file) > 1:
    if len(file) != len(outnames):
        print("Mismatch in number of input files or output names, check arguments")
        exit()

#Check that a sensible number of threads has been requested - more protections here are possible - at the moment users are trusted to be sensible
cpucount = multiprocessing.cpu_count()
if args.threads is None:
    print("Defaulting to one thread")
    threads = 1
elif int(args.threads) >= cpucount:
    print("Too many threads requested, resetting to default")
    threads = 1
else:
    print(f"Running alignment with {threads} threads...")

#Check that the user has entered a valid genome to align to
if genome == "mm10" or genome == "mm39":
    print(f"Primary genome is mouse genome {genome}")
elif genome == "hg19" or genome == "hg38":
    print(f"Primary genome is human genome {genome}")
elif genome == "sacCer3":
    print(f"Primary genome is yeast genome {genome}")
else:
    print("Genome option not recognised or not entered. Please use mm10/39 or hg19/38 or ask Miles to change the script to accommodate your new organism/genome")
    exit()

#Process ID (used to make unique sorttemp, so these are not overlapping for multiple processes in the same outdir)
uniqueid = str(uuid.uuid4())

#Create a place to store counts
allcountslist = []

if spikegenome is None:
    spikecheck = 0
    print(f"No spikein")
else:
    spikecheck = 1
    print(f"Spike genome is {spikegenome}")

#Process the files
if spikecheck:
    for i, fastq in enumerate(file):
        name = outnames[i]
        print(f"Aligning {name} to {genome} and {spikegenome}")
        #Explanation: align to genome and spikein genome, then remove multiply aligned reads with grep (XS: indicates the score of the next best aligning read, if it exists), then make a bam of all mapped reads (-F4 removes reads with a SAM flag of 4, which means unmapped)
        sp.run(f"bowtie2 -p {threads} --no-mixed --no-discordant -U {fastq} -x /mnt/md0/DataRepository/genomes/{genome}.{spikegenome}/{genome}.{spikegenome} 2> {outdir}{name}_alignment_stats.txt | grep -v XS: - | samtools view -bh -F4 - > {outdir}{name}_UniqMapped.bam", shell=True, executable="/bin/bash")
        #Sort aligned reads
        sp.run(f"sambamba sort --tmpdir {outdir}{uniqueid}/ -t {threads} -m 30G -o {outdir}{name}_UniqMapped_sorted.bam {outdir}{name}_UniqMapped.bam", shell=True, executable="/bin/bash")
        #Remove duplicates
        sp.run(f"sambamba markdup --tmpdir {outdir}{uniqueid}/ -r -t {threads} {outdir}{name}_UniqMapped_sorted.bam {outdir}{name}_UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")

        #Next, need to separate out reads from each genome:
        print(f"Extracting reads aligning uniquely to {genome}")
        #Use samtools view to open the file, grep to remove those with the spikein genome name in them, then change the chromosomes with the genome name to just the chromosome numbers as normally used, then put the file back into bam.
        sp.run(f"samtools view -h {outdir}{name}_UniqMapped_sorted_rmdup.bam | grep -v {spikegenome} | sed s/{genome}_chr/chr/g | samtools view -bhS - > {outdir}{name}_{genome}.UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
        #Do the same thing for the spikein
        print(f"Extracting reads aligning uniquely to {spikegenome}.")
        sp.run(f"samtools view -h {outdir}{name}_UniqMapped_sorted_rmdup.bam | grep -v {genome} | sed s/{spikegenome}_chr/chr/g | samtools view -bhS - > {outdir}{name}_{spikegenome}.UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
        #Index outputs
        sp.run(f"sambamba index -t {threads} {outdir}{name}_UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
        sp.run(f"sambamba index -t {threads} {outdir}{name}_{genome}.UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
        sp.run(f"sambamba index -t {threads} {outdir}{name}_{spikegenome}.UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
        #Now clean up
        sp.run(f"rm {outdir}{name}_UniqMapped.bam", shell=True, executable="/bin/bash")
        sp.run(f"rm {outdir}{name}_UniqMapped_sorted.bam", shell=True, executable="/bin/bash")

        #Finally, count reads in each file:
        totalCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_UniqMapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")
        genomeCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_{genome}.UniqMapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")
        spikegenomeCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_{spikegenome}.UniqMapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")
        countsList = [totalCount.stdout.decode('ascii').strip(), genomeCount.stdout.decode('ascii').strip(), spikegenomeCount.stdout.decode('ascii').strip()]
        allcountslist.append(countsList)

    sp.run(f"rm {outdir}{uniqueid}/ -r", shell=True, executable="/bin/bash")
    countstable = pd.DataFrame(allcountslist, columns = ['allcounts', 'genomecounts', 'spikecounts'])
    outtable = pd.concat([files, countstable], axis=1)
    outtable.to_csv(outdir + outname + '.tsv', sep = '\t', index = False, header = True)


else:
    #No spikein version
    for i, fastq in enumerate(file):
        name = outnames[i]
        print(f"Aligning {name} to {genome}")
        #Explanation: align to genome, then remove multiply aligned reads with grep (XS: indicates the score of the next best aligning read, if it exists), then make a bam of all mapped reads (-F4 removes reads with a SAM flag of 4, which means unmapped)
        sp.run(f"bowtie2 -p {threads} --no-mixed --no-discordant -U {fastq} -x /mnt/md0/DataRepository/genomes/{genome}/{genome} 2> {outdir}{name}_alignment_stats.txt | grep -v XS: - | samtools view -bh -F4 - > {outdir}{name}_UniqMapped.bam", shell=True, executable="/bin/bash")
        #Sort aligned reads
        sp.run(f"sambamba sort --tmpdir {outdir}{uniqueid}/ -t {threads} -m 30G -o {outdir}{name}_UniqMapped_sorted.bam {outdir}{name}_UniqMapped.bam", shell=True, executable="/bin/bash")
        #Remove duplicates
        sp.run(f"sambamba markdup --tmpdir {outdir}{uniqueid}/ -r -t {threads} {outdir}{name}_UniqMapped_sorted.bam {outdir}{name}_UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")

        #Index outputs
        sp.run(f"sambamba index -t {threads} {outdir}{name}_UniqMapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")

        #Now clean up
        sp.run(f"rm {outdir}{name}_UniqMapped.bam", shell=True, executable="/bin/bash")
        sp.run(f"rm {outdir}{name}_UniqMapped_sorted.bam", shell=True, executable="/bin/bash")

        #Finally, count reads in each file:
        totalCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_UniqMapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")

        countsList = [totalCount.stdout.decode('ascii').strip()]
        allcountslist.append(countsList)

    sp.run(f"rm {outdir}{uniqueid}/ -r", shell=True, executable="/bin/bash")
    countstable = pd.DataFrame(allcountslist, columns = ['allcounts'])
    outtable = pd.concat([files, countstable], axis=1)
    outtable.to_csv(outdir + outname + '.tsv', sep = '\t', index = False, header = True)
