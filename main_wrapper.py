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
import glob
from datetime import datetime


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
addTags = False
keeptmp = False

nOptions=0

def usage(help=False):
    msg = '''

Usage: ./bamrefine [options] <in.bam> <out.bam>

OPTIONS:
        -s, --snps                    BED  or SNP formatted file for snps
        -p, --threads                 # of threads to use
        -l, --pmd-length-threshold    pmd length threshold
FLAGS:
        -t, --add-tags                add maskings stats as optional SAM fields to the alignments
        -v, --verbose                 verbose output of progress
        -k, --keep-tmp                don't remove the temporary directory (./.tmp_bamrefine)
        -h, --help                    display this message end exit


    '''

    helpmsg = '''
Consider reading the man page. You can do so by navigating to the
program's directory and run 'man ./bamrefine.1'
'''

    if help:
        msg = helpmsg + msg

    print(msg)
    sys.exit()

try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:],
                                           's:p:l:tvhk',
                                           ['snps=',
                                           'threads=',
                                            'pmd-length-threshold=',
                                            'add-tags',
                                           'verbose',
                                            'help',
                                           'keep-tmp'])
except getopt.GetoptError as err:
    print(str(err))
    usage()


for opt, arg in options:
    if opt in ('-p', '--threads'):
        thread = int(arg)
        nOptions += 1
    elif opt in ('-s', '--snps'):
        snpF = arg
        nOptions += 1
    elif opt in ('-l', '--pmd-length-threshold'):
        lookup = arg
        nOptions += 1
    elif opt in ('-t', '--add-tags'):
        addTags = True
    elif opt in ('-v', '--verbose'):
        verbose = True
    elif opt in ('-k', '--keep-tmp'):
        keeptmp = True
    elif opt in ('-h', '--help'):
        usage(help=True)
    else:
        usage()


if nOptions != 3:
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

inName  = os.path.abspath(inName)
ouName  = os.path.abspath(ouName)
snpF    = os.path.abspath(snpF)
ouDir   = os.path.dirname(ouName)


## ----------------------------------------------------------------
def waitJobs(runningJobs):
    while len([i for i in runningJobs if i.poll() != None]) < len(runningJobs):
        continue

def parallelParse(jobL, n, lookup):
    activeJobs = []
    jobN = []
    global snpF

    for i in range(n):
        try:
            c = jobL.pop()
        except IndexError:
            waitJobs(activeJobs)
            return None

        cmdList = ["python3", dirN+ "/main.py",inName,c, lookup, snpF, str(int(addTags))]
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
                if verbose:
                    print("Finished job for chr%s" % jobN[i])
                if len(jobL) == 0:
                    continue
                c = jobL.pop()
                cmdList = ["python3", dirN+ "/main.py",inName,c, lookup, snpF,str(int(addTags))]
                cmd = " ".join(cmdList)
                p = Popen([cmd], shell = True)
                if verbose:
                    print("Started job for chr%s" % c)
                activeJobs[i] = p
                jobN[i] = c

    waitJobs(activeJobs)


if os.path.isfile(inName+".bai"):
    print("Fetching Chromosomes")
    chrms = bamRefine_cy.fetchChromosomes(inName)
    print("Done.")
elif os.path.isfile(inName):
    msg = '''
    input BAM is not currently indexed. Indexing...
    '''
    print(msg)
    pysam.index("-@ " + str(thread), inName)
    print("Done.")
    chrms = bamRefine_cy.fetchChromosomes(inName)
    print("Fetching Chromosomes")
    print("Done.")
else:
    msg = '''
Can't find input BAM in the specified path.
    '''
    print(msg)
    exit()

## Create a tmp directory
try:
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    oubase = os.path.basename(ouName)
    tmpname = "".join([ouDir, "/.", now,"_", oubase, "_tmp_bamrefine"])
    os.mkdir(tmpname)
except FileExistsError:
    shutil.rmtree(tmpname)
    os.mkdir(tmpname)


print('\nStarted bam filtering\n')
os.chdir(tmpname)

jobs, bypass = bamRefine_cy.createBypassBED(inName, chrms, snpF)
jobs_c = jobs.copy()

parallelParse(jobs, thread, lookup)

print("Finished BAM filtering\n\nMerging BAM files...")



## samtools merge commands -------------
with open('toMerge_bamlist.txt', 'w') as f:
    for c in jobs_c:
        f.write(c+'.bam\n')
    if bypass:
        f.write('bypassed.bam\n')

toMergeF = 'toMerge_bamlist.txt'
### ----------------------

if bypass:
    pysam.view("-@", str(thread), "--no-PG", "-L", "bypass.bed", "-b", "-o" "bypassed.bam", inName, catch_stdout=False) ##pysam bug

pysam.merge("--no-PG", "-c", "-p", "-b" , toMergeF, "-O", "BAM", "-@", str(thread), ouName)

# pysam.view("-H", "-o", "header.sam", inName, catch_stdout=False)
# pysam.reheader("-P",  "-i", "header.sam", ouName)

# while reheading.poll() == None:
#     continue

print("Finished merging.")

### -----------------------------------

## Wrapping up stats -----------------
stats = []
for statFile in glob.glob('./*_stats.txt'):
    with open(statFile) as f:
        stats.append(f.readline().strip().split('\t'))

stats_5 = sum([int(x[0]) for x in stats])
stats_3 = sum([int(x[1]) for x in stats])

stats_msg = "\nmasked %d 5' and %d 3' positions in total"

print(stats_msg % (stats_5, stats_3))

stats_cmd = "cat *_stats.txt > all_stats.txt ; "
stats_cmd += "ls *_stats.txt | grep -v all | sed 's/_stats.txt//g' > all_names.txt ; "
stats_cmd += "paste all_names.txt all_stats.txt | "
stats_cmd += "sort -n -k1 > " + ouName + "_bamrefine_stats.txt"

os.system(stats_cmd)

### ----------------------------------

## Cleaning up temp files -------------
os.chdir('../')

if not keeptmp:
    shutil.rmtree(tmpname)
### ----------------------------------

print("Program finished successfully.")

