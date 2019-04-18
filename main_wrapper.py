import pysam

chrms = range(1,25)

for c in chrms:
    cmd = "python3 main.py " + c
    os.system(cmd)


merge_cmd = "samtools merge -b " + toMerge + " -O BAM merged_out.bam"
