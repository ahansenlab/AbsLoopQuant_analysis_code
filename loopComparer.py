#Script to take loops from two files and produce 6 output files: loops common to both conditions, loops unique to each condition, and anchors common to both conditions, and anchors unique to each condition.
#This uses shell calls to bedtools pairToPair and intersectBed to carry out all the sequential steps needed to get these files, making it slightly more streamlined and less accidental typo prone.

from sys import exit
import subprocess as sp
import argparse

parser = argparse.ArgumentParser(description = "compare loops and optionally loop anchors from two sets of loops - NOTE: make sure bedpe files do not have columns beyond the required 6, and do not have headers. Remove these if they are present.")
parser.add_argument("--bedpe1", "-1", help = "first bedpe file to compare - required")
parser.add_argument("--bedpe2", "-2", help = "second bedpe file to compare - required")
parser.add_argument("--names", "-n", help = "names of conditions to use - will be used to make the output filenames", nargs = 2)
parser.add_argument("--outdir", "-o", help = "a directory to store output files - default is current directory", default = "./")
parser.add_argument("--anchors", "-a", help = "make comparisons between anchors in each file - default is to skip. Note that anchors are merged in resulting files so each anchor only appears once at most", action = "store_true")
parser.add_argument("--slop", "-s", help = "basepairs of slop to add to loop anchors when comparing - default is none")
parser.add_argument("--genome", "-g", help = "if using slop and comparing anchors, genome is required to avoid going off the ends of the chromosomes - note the the comparison of whole loops does not take this into account, so slop sensibly")
args = parser.parse_args()

bedpe1 = args.bedpe1
bedpe2 = args.bedpe2
name1 = args.names[0]
name2 = args.names[1]
slop = args.slop
anchors = args.anchors
outdir = args.outdir
genome = args.genome

if args.bedpe1 is None or args.bedpe2 is None:
    print("Input files not specified - check help for formatting")
    parser.print_usage()
    exit()

#Check that outdir ends with a /, add one if it doesn't
if args.outdir is not None and not outdir.endswith("/"):
    outdir = outdir + "/"

#Make sure the chromsizes file can be accessed if slop is needed for the anchors
if anchors and slop is not None:
    if genome is None:
        print("Genome not specified - check help for formatting")
        parser.print_usage()
        exit()
    elif genome == "mm10" or genome == "mm39":
        print("Using mouse chromosomes...")
    elif genome == "hg19" or genome == "hg38":
        print("Using human chromosomes...")
    else:
        print("Genome option not recognised or not entered. Please use mm10/39 or hg19/38 or ask Miles to change the script to accommodate your new organism/genome, and make sure your chromosome sizes file is placed in the corresponding location and named correctly (see existing files for guidance).")
        exit()
elif anchors:
    genome = ""

#Make strings necessary to slop files if requested
if slop is not None:
    pairslopst = " -slop " + str(slop)
else:
    pairslopst = ""

#Commands for pair comparisons
bothline = "pairToPair{4} -a {0} -b {1} | cut -f 1-6 > {5}{2}and{3}.bedpe && wc -l {0} && wc -l {1} && wc -l {5}{2}and{3}.bedpe"
firstonlyline = "pairToPair{4} -a {0} -b {1} | cut -f 1-6 | grep -F -v -x -f /dev/stdin {0} | cut -f 1-6 > {5}{2}not{3}.bedpe && wc -l {5}{2}not{3}.bedpe"
secondonlyline = "pairToPair{4} -a {1} -b {0} | cut -f 1-6 | grep -F -v -x -f /dev/stdin {1} | cut -f 1-6 > {5}{3}not{2}.bedpe && wc -l {5}{3}not{2}.bedpe"

bedpelines = [bothline, firstonlyline, secondonlyline]

#run it
for line in bedpelines:
    tokenized_line = line.format(bedpe1, bedpe2, name1, name2, pairslopst, outdir)
    #print(tokenized_line)
    # run
    sp.run(tokenized_line, shell=True, executable="/bin/bash")

#Set up comparison for anchors as well.
#Only do if anchor option set:
if anchors:
    #First, need to make each bedpe into a bed format - take columns 4-6 and cat them onto 1-3 (or 3-5 onto 0-2, if you index from 0), then remove the intermediates - this should basically be instant
    cutbedpe1line = "cut -f 1-3 {0} > {1}{2}.firstcut.bedpe && cut -f 4-6 {0} > {1}{2}.secondcut.bedpe && cat {1}{2}.firstcut.bedpe {1}{2}.secondcut.bedpe > {1}{2}.anchors.bed && rm {1}{2}.firstcut.bedpe {1}{2}.secondcut.bedpe"
    if slop is not None:
        bed1slopline = "slopBed -i {1}{2}.anchors.bed -g /mnt/md0/DataRepository/chromsizes/{5}/{5}.sorted.chrom.sizes -b " + str(slop) + " > {1}{2}.anchors.slopped.bed && rm {1}{2}.anchors.bed && mv {1}{2}.anchors.slopped.bed {1}{2}.anchors.bed"
    else:
        bed1slopline = ""
    cutbedpe2line = "cut -f 1-3 {3} > {1}{4}.firstcut.bedpe && cut -f 4-6 {3} > {1}{4}.secondcut.bedpe && cat {1}{4}.firstcut.bedpe {1}{4}.secondcut.bedpe > {1}{4}.anchors.bed && rm {1}{4}.firstcut.bedpe {1}{4}.secondcut.bedpe"
    if slop is not None:
        bed2slopline = "slopBed -i {1}{4}.anchors.bed -g /mnt/md0/DataRepository/chromsizes/{5}/{5}.sorted.chrom.sizes -b " + str(slop) + " > {1}{4}.anchors.slopped.bed && rm {1}{4}.anchors.bed && mv {1}{4}.anchors.slopped.bed {1}{4}.anchors.bed"
    else:
        bed2slopline = ""
    makebed = [cutbedpe1line, bed1slopline, cutbedpe2line, bed2slopline]
    for line in makebed:
        tokenized_line = line.format(bedpe1, outdir, name1, bedpe2, name2, genome)
        #print(tokenized_line)
        # run
        sp.run(tokenized_line, shell=True, executable="/bin/bash")
    #Use intersectBed from bedtools to look for matching/not matching, and merge after this so each anchor shows up only once.
    bothbedline = "intersectBed -u -f 0.1 -a {1}{2}.anchors.bed -b {1}{4}.anchors.bed > {1}{2}and{4}.anchors.bed && sortBed -i {1}{2}and{4}.anchors.bed > {1}{2}and{4}.anchors.sorted.bed && mergeBed -i {1}{2}and{4}.anchors.sorted.bed -d -1 > {1}{2}and{4}.anchors.merged.bed && wc -l {1}{2}and{4}.anchors.merged.bed && rm {1}{2}and{4}.anchors.bed {1}{2}and{4}.anchors.sorted.bed && mv {1}{2}and{4}.anchors.merged.bed {1}{2}and{4}.anchors.bed"
    bed1onlyline = "intersectBed -a {1}{2}.anchors.bed -b {1}{4}.anchors.bed -v > {1}{2}not{4}.anchors.bed && sortBed -i {1}{2}not{4}.anchors.bed > {1}{2}not{4}.anchors.sorted.bed && mergeBed -i {1}{2}not{4}.anchors.sorted.bed -d -1 > {1}{2}not{4}.anchors.merged.bed && wc -l {1}{2}not{4}.anchors.merged.bed && rm {1}{2}not{4}.anchors.bed {1}{2}not{4}.anchors.sorted.bed && mv {1}{2}not{4}.anchors.merged.bed {1}{2}not{4}.anchors.bed"
    bed2onlyline = "intersectBed -a {1}{4}.anchors.bed -b {1}{2}.anchors.bed -v > {1}{4}not{2}.anchors.bed && sortBed -i {1}{4}not{2}.anchors.bed > {1}{4}not{2}.anchors.sorted.bed && mergeBed -i {1}{4}not{2}.anchors.sorted.bed -d -1 > {1}{4}not{2}.anchors.merged.bed && wc -l {1}{4}not{2}.anchors.merged.bed && rm {1}{4}not{2}.anchors.bed {1}{4}not{2}.anchors.sorted.bed && mv {1}{4}not{2}.anchors.merged.bed {1}{4}not{2}.anchors.bed"
    bedcompare = [bothbedline, bed1onlyline, bed2onlyline]
    for line in bedcompare:
        tokenized_line = line.format(bedpe1, outdir, name1, bedpe2, name2, genome)
        #print(tokenized_line)
        # run
        sp.run(tokenized_line, shell=True, executable="/bin/bash")
