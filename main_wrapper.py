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

def usage():
    msg = '''
Usage: ./bamrefine [options]

OPTIONS:
        -i, --input:        Input BAM file
        -o, --output:       Output BAM file [out.bam]
        -s, --snps:         Snp collection file [chr_pos_ref_alt_1240K.all.snp]
        -c, --chromosomes:  List of chromosomes.
        -p, --threads:      # of threads to use [2]
    '''

    print(msg)
    sys.exit()

options, remainder = getopt.gnu_getopt(sys.argv[1:],
                                       'i:s:c:o:p:k',
                                       ['input=',
                                        'snps=',
                                        'chromosomes=',
                                       'output=',
                                       'threads=',
                                       'keeptmp'])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inName = arg
    elif opt in ('-o', '--output'):
        ouName = arg
    elif opt in ('-p', '--threads'):
        thread = int(arg)
    elif opt in ('-s', '--snps'):
        snpF = arg
    elif opt in ('-c', '--chromosomes'):
        chrmF = arg
    else:
        usage()

## default args
chrmF = 'chrms.txt'
thread = 2
inName = None
ouName = 'out.bam'
## ----------------------------------------------------------------

def parallelParse(jobL, n):
    activeJobs = []
    jobN = []
    for i in range(n):
        c = jobL.pop()
        cmd = "python3 main.py " + inName + " " + c
        p = Popen([cmd], shell = True)
        activeJobs.append(p)
        jobN.append(c)
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
                cmd = "python3 main.py " + inName + " " + c
                p = Popen([cmd], shell = True)
                print("Finished job for chr%s" % jobN[i])
                print("Started job for chr%s" % c)
                activeJobs[i] = p
                jobN[i] = c



    while len([i for i in range(n) if activeJobs[i].poll() != None]) < len(activeJobs):
        continue

def handleSNPs(fName):
    pickleName = fName + ".brf"
    f_exists = os.path.isfile(pickleName)
    if f_exists:
        f = open(pickleName, 'rb')
        snps = pickle.load(f)
        f.close()
    else:
        snps = bamRefine_cy.parseSNPs('chr_pos_ref_alt_1240K.all.snp')
        f = open(pickleName, 'wb')
        pickle.dump(snps, f)
        f.close()

    return snps


snps = handleSNPs(snpF)

del(bamRefine_cy)


jobs = chrms.copy()
parallelParse(jobs, thread)

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

print("Finished merging.")

cleanup_cmd = "for i in `cat toMerge_bamlist.txt`; do rm $i ;done"

print("Cleaning up temp. files...")

cleanup = Popen([cleanup_cmd], shell = True)

while cleanup.poll() == None:
    continue

print("Program finished successfully.")

