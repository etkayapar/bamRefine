import sys
from bamRefine_cy import *
import os
# import pysam


inName = sys.argv[1]
contig = sys.argv[2]
lookup = sys.argv[3]
snps = sys.argv[4]


os.chdir('tmp_bamrefine')

snps = handleSNPs(snps)
# snps = list(snps.values())
# snps = snps[:3000]

ouName = contig + ".bam"

inBAM = pysam.AlignmentFile(inName, 'rb')
ouBAM = pysam.AlignmentFile(ouName, 'wb', template=inBAM)

processBAM(inBAM, ouBAM, snps, contig, lookup)
