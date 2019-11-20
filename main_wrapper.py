#! /usr/bin/env python3

import sys
import os
import pysam
import time
from subprocess import Popen
import getopt
import bamRefine_cy
import pickle
import shutil


dirN = os.path.dirname(os.path.realpath(__file__))#dirname of the script
# Processing command-line options for the program ------------------

## default args
chrmFname = None
genomeF = None
thread = None
inName = None
ouName = None
lookup = None
verbose = False
snpF = None


def usage():
    msg = '''
Usage: ./bamrefine [options] <in.bam> <out.bam>

OPTIONS:
        -s, --snps                    BED  or SNP formatted file for snps
        -p, --threads                 # of threads to use
        -g, --ref-genome              Path to ref. genome to fetch chr/contig names
        -l, --pmd-length-threshold    pmd length threshold
FLAGS:
        -h, --help                    display this message end exit
        -v, --verbose                 verbose output of progress
    '''

    print(msg)
    sys.exit()

try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:],
                                           's:p:g:l:vh',
                                           ['snps=',
                                           'threads=',
                                            'ref-genome',
                                            'pmd-length-threshold'
                                           'verbose',
                                            'help'])
except getopt.GetoptError as err:
    print(str(err))
    usage()


for opt, arg in options:
    if opt in ('-p', '--threads'):
        thread = int(arg)
    elif opt in ('-s', '--snps'):
        snpF = arg
    elif opt in ('-l', '--pmd-length-threshold'):
        lookup = arg
    elif opt in ('-g', '--ref-genome'):
        genomeF = arg
    elif opt in ('-v', '--verbose'):
        verbose = True
    elif opt in ('-h', '--help'):
        usage()
    else:
        usage()


if len(options) < 4:
    print("All options need arguments")
    usage()

if remainder[0].startswith('/'):
    inName = remainder[0]
else:
    inName = './'+remainder[0]

if remainder[1].startswith('/'):
    ouName = remainder[1]
else:
    ouName = './'+remainder[1]

if snpF.startswith('/'):
    pass
else:
    snpF = './'+snpF

if genomeF.startswith('/'):
    pass
else:
    genomeF = './'+genomeF

inName  = os.path.abspath(inName)
ouName  = os.path.abspath(ouName)
snpF    = os.path.abspath(snpF)
genomeF = os.path.abspath(genomeF)



## ----------------------------------------------------------------

def parallelParse(jobL, n, lookup):
    activeJobs = []
    jobN = []
    global snpF

    ## careful with this part ------------------------------------

    '''
    Run the first contig independent of the parallelism
    for the SNP pickle file not being written to and reading
    from at the same time, it causes errors.
    '''

    firstC = jobL.pop()
    cmdList = ["python3", dirN+ "/main.py",inName,firstC, lookup, snpF]
    cmd = " ".join(cmdList)
    firstCjob = Popen([cmd], shell = True)
    ## activeJobs.append(p)
    ## jobN.append(c)

    while firstCjob.poll() != None: # wait until it finishes
        continue

    del([cmdList, cmd, firstCjob, firstC])

    ### ----------------------------------------------------------

    for i in range(n):
        c = jobL.pop()
        cmdList = ["python3", dirN+ "/main.py",inName,c, lookup, snpF]
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
                cmdList = ["python3", dirN+ "/main.py",inName,c, lookup, snpF]
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
elif os.path.isfile(genomeF+".fai"):
    print("Fetching Chromosomes from fasta...")
    chrmFname = genomeF+".fai"
    with open(chrmFname) as chrmF:
        chrms = [x.strip().split()[0] for x in chrmF.readlines()]
    print("Done.")
elif os.path.isfile(genomeF):
    msg = '''
Genome fasta is not indexed. Please index your genome with
samtools faidx
    '''
else:
    msg = '''
Can't find genome fasta in the specified path.
    '''
    print(msg)
    exit()

jobs = chrms.copy()
## Create a tmp directory
try:
    os.mkdir('.tmp_bamrefine')
except FileExistsError:
    shutil.rmtree('.tmp_bamrefine')
    os.mkdir('.tmp_bamrefine')


print('Started bam filtering\n')
os.chdir('.tmp_bamrefine')
parallelParse(jobs, thread, lookup)

print("Please wait...")
time.sleep(10)

print("Finished BAM filtering\nMerging BAM files...")



## samtools merge commands -------------
with open('toMerge_bamlist.txt', 'w') as f:
    for c in chrms:
        f.write(c+'.bam\n')

toMergeF = 'toMerge_bamlist.txt'
merge_cmd = "samtools merge -b " + toMergeF + " -O BAM -@ "
merge_cmd += " ".join([str(thread),ouName])

## debugging stuff --------
print(os.getcwd())
print(merge_cmd)

### ----------------------

merging = Popen([merge_cmd], shell = True)

while merging.poll() == None:
    continue

rehead_cmd = "samtools view -H " + inName + " | samtools reheader -P  - " + ouName
rehead_cmd += " > rehead.bam" + " ; mv rehead.bam " + ouName

print(rehead_cmd)

reheading = Popen([rehead_cmd], shell = True)

while reheading.poll() == None:
    continue

print("Finished merging.")

### -----------------------------------
## Cleaning up temp files -------------
os.chdir('../')

shutil.rmtree('.tmp_bamrefine')

### ----------------------------------

print("Program finished successfully.")

