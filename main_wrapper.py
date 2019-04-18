import pysam
from subprocess import Popen

chrms = [str(x) for x in range(1,23)]
chrms += ['X', 'Y']


procs = []
for c in chrms:
    cmd = "python3 main.py " + c
    p = Popen([cmd], shell = True)
    procs.append(p)

while len([x.poll() for x in procs if x.poll() != None]) != len(procs):
    continue

print("Merging BAM files...")

with open('toMerge_bamlist.txt', 'w') as f:
    for c in chrms:
        f.write(c+'.bam\n')

toMergeF = 'toMerge_bamlist.txt'


merge_cmd = "samtools merge -b " + toMergeF + " -O BAM -@ 8 merged_out.bam"

merging = Popen([merge_cmd], shell = True)

while merging.poll() == None:
    continue

print("Finished merging.\nProgram executed successfully")
