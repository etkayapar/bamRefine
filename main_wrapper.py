#! /usr/bin/env python3

import sys
import os
import pysam
import time
from subprocess import Popen
import getopt
import bamRefine_cy
import pickle


# Processing command-line options for the program ------------------

## default args
chrmFname = 'chrms.txt'
genomeF = None 
thread = 2
inName = None
ouName = 'out.bam'
lookup = 10
verbose = False
snpF = 'chr_pos_ref_alt_1240K.all.snp'


def usage():
    msg = '''
Usage: ./bamrefine [options]

OPTIONS:
        -i, --input
                      Input BAM file
        -o, --output
                      Output BAM file [out.bam]
        -s, --snps
                      BED formatted file for snps [chr_pos_ref_alt_1240K.all.snp]
        -p, --threads
                      # of threads to use [2]
        -g, --ref-genome
                      Path to ref. genome to fetch chr/contig names or, a .txt
                      file containing chr names
        -l, --pmd-length-threshold
                      pmd length threshold [10]
        -v, --verbose
                      verbose output of progress
    '''

    print(msg)
    sys.exit()

try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:],
                                           'i:s:o:p:g:l:v',
                                           ['input=',
                                            'snps=',
                                           'output=',
                                           'threads=',
                                            'ref-genome',
                                            'pmd-length-threshold'
                                           'verbose'])
except getopt.GetoptError as err:
    print(str(err))
    usage()

if len(options) < 1:
    usage()

for opt, arg in options:
    if opt in ('-i', '--input'):
        inName = arg
    elif opt in ('-o', '--output'):
        ouName = arg
    elif opt in ('-p', '--threads'):
        thread = int(arg)
    elif opt in ('-s', '--snps'):
        snpF = arg
    elif opt in ('-l', '--pmd-length-threshold'):
        lookup = arg
    elif opt in ('-g', '--ref-genome'):
        genomeF = arg
    elif opt in ('-v', '--verbose'):
        verbose = True
    else:
        usage()

## ----------------------------------------------------------------

def parallelParse(jobL, n, lookup):
    activeJobs = []
    jobN = []
    global snpF
    for i in range(n):
        c = jobL.pop()
        cmdList = ["python3 main.py",inName,c, lookup, snpF]
        cmd = " ".join(cmdList)
        p = Popen([cmd], shell = True)
        activeJobs.append(p)
        jobN.append(c)
        if verbose:
            print('Started job for chr%s' % c)

    while True:
        time.sleep(3)
        if len(jobL) == 0:
            break
        # get list of all finished processes
        finished = [i for i in range(n) if activeJobs[i].poll() != None]
        n_f = len(finished)
        if n_f > 0:
            for i in finished:
                if len(jobL) == 0:
                    continue
                c = jobL.pop()
                cmdList = ["python3 main.py",inName,c, lookup, snpF]
                cmd = " ".join(cmdList)
                p = Popen([cmd], shell = True)
                if verbose:
                    print("Finished job for chr%s" % jobN[i])
                    print("Started job for chr%s" % c)
                activeJobs[i] = p
                jobN[i] = c



    while len([i for i in range(n) if activeJobs[i].poll() != None]) < len(activeJobs):
        continue

if genomeF == None:
    with open(chrmFname) as chrmF:
        chrms = [x.strip() for x in chrmF.readlines()]
else:
    get_chr_cmdList = ['./getchrms.sh', genomeF, '> ./tmp.chr'] 
    get_chr_cmd = " ".join(get_chr_cmdList)
    fetching = Popen([get_chr_cmd], shell = True)
    if verbose:
        print("Fetching Chromosomes from fasta...")
    while fetching.poll() == None:
        continue
    chrmFname = 'tmp.chr'
    with open(chrmFname) as chrmF:
        chrms = [x.strip() for x in chrmF.readlines()]
    if verbose:
        print("Done.")


jobs = chrms.copy()

print('Started bam filtering\n')
parallelParse(jobs, thread, lookup)

print("Please wait...")
time.sleep(10)

print("Finished BAM filtering\nMerging BAM files...")

with open('toMerge_bamlist.txt', 'w') as f:
    for c in chrms:
        f.write(c+'.bam\n')

toMergeF = 'toMerge_bamlist.txt'


merge_cmd = "samtools merge -b " + toMergeF + " -O BAM -@ "
merge_cmd += str(thread) + " " +  ouName

merging = Popen([merge_cmd], shell = True)

while merging.poll() == None:
    continue

rehead_cmd = "samtools view -H " + inName + " | samtools reheader -P  - " + ouName
rehead_cmd += " > rehead.bam" + " ; mv rehead.bam " + ouName

reheading = Popen([rehead_cmd], shell = True)

while reheading.poll() == None:
    continue

print("Finished merging.")

cleanup_cmd = "for i in `cat toMerge_bamlist.txt`; do rm $i ;done" 
cleanup_cmd += " ; rm toMerge_bamlist.txt"

print("Cleaning up temp. files...")

cleanup = Popen([cleanup_cmd], shell = True)

while cleanup.poll() == None:
    continue

print("Program finished successfully.")

