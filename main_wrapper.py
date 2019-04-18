import sys
import pysam
import time
from subprocess import Popen

chrms = [str(x) for x in range(1,23)]
chrms += ['X', 'Y']

inName = sys.argv[1]
thread = int(sys.argv[2])

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
            for i in range(n_f):
                if len(jobL) == 0:
                    continue
                c = jobL.pop() 
                cmd = "python3 main.py " + inName + " " + c
                p = Popen([cmd], shell = True)
                print("Finished job for chr%s" % jobN[i])
                activeJobs[i] = p
                jobN[i] = c


jobs = chrms.copy()
parallelParse(jobs, thread)

time.sleep(3)


print("Finished BAM filtering\nMerging BAM files...")

with open('toMerge_bamlist.txt', 'w') as f:
    for c in chrms:
        f.write(c+'.bam\n')

toMergeF = 'toMerge_bamlist.txt'


merge_cmd = "samtools merge -b " + toMergeF + " -O BAM -@ " + str(thread) +" merged_out.bam"

merging = Popen([merge_cmd], shell = True)

while merging.poll() == None:
    continue

print("Finished merging.\nProgram executed successfully")
