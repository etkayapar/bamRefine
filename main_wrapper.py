import pysam
from subprocess import Popen

chrms = range(1,25)

procs = []
for c in chrms:
    cmd = "python3 main.py " + c
    p = Popen([cmd], shell = True)
    procs.append(p)

while len([x.poll() for x in procs if x.poll() != None]) != len(procs):
    continue

print("Merging BAM files...")



merge_cmd = "samtools merge -b " + toMerge + " -O BAM merged_out.bam"

merging = Popen([merge_cmd], shell = True)

while merging.poll() == None:
    continue

print("Finished merging.\nProgram executed successfully")
